/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center, 
Division of Medical and Biological Informatics.
All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without 
even the implied warranty of MERCHANTABILITY or FITNESS FOR 
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/

#ifndef mitkDicomSeriesReader_h
#define mitkDicomSeriesReader_h

#include "mitkDataNode.h"
#include "mitkConfig.h"

#include <itkGDCMImageIO.h>

#include <itkImageSeriesReader.h>
#include <itkCommand.h>

#ifdef NOMINMAX
#  define DEF_NOMINMAX
#  undef NOMINMAX
#endif

#include <gdcmConfigure.h>

#ifdef DEF_NOMINMAX
#  ifndef NOMINMAX
#    define NOMINMAX
#  endif
#  undef DEF_NOMINMAX
#endif

#include <gdcmDataSet.h>
#include <gdcmRAWCodec.h>
#include <gdcmSorter.h>
#include <gdcmScanner.h>
#include <gdcmPixmapReader.h>
#include <gdcmStringFilter.h>

namespace mitk
{

/** 
 \brief Loading DICOM images as MITK images.

 - \ref DicomSeriesReader_purpose
 - \ref DicomSeriesReader_limitations
 - \ref DicomSeriesReader_usage
 - \ref DicomSeriesReader_sorting
   - \ref DicomSeriesReader_sorting1
   - \ref DicomSeriesReader_sorting2
   - \ref DicomSeriesReader_sorting3
   - \ref DicomSeriesReader_sorting4
 - \ref DicomSeriesReader_tests

 \section DicomSeriesReader_purpose Purpose

 DicomSeriesReader serves as a central class for loading DICOM images as mitk::Image.
 As the term "DICOM image" covers a huge variety of possible modalities and 
 implementations, and since MITK assumes that 3D images are made up of continuous blocks
 of slices without any gaps or changes in orientation, the loading mechanism must
 implement a number of decisions and compromises.

 <b>The main intention of this implementation is not efficiency but correcness of generated slice positions!</b>

 \section DicomSeriesReader_limitations Assumptions and limitations

 The class is working only with GDCM 2.0.14 (or possibly newer). This version is the
 default of an MITK super-build. Support for other versions or ITK's DicomIO was dropped
 because of the associated complexity of DicomSeriesReader.

 \b Assumptions
  - expected to work for SOP Classes CT Image Storage and MR Image Storage (NOT for the "Enhanced" variants containing multi-frame images)
  - special treatment for a certain type of Philips 3D ultrasound (recogized by tag 3001,0010 set to "Philips3D")
  - loader will always attempt to read multiple single slices as a single 3D image volume (i.e. mitk::Image)
    - slices will be grouped by basic properties such as orientation, rows, columns, spacing and grouped into as large blocks as possible
  
 \b Options
  - images that cover the same piece of space (i.e. position, orientation, and dimensions are equal)
    can be interpreted as time-steps of the same image, i.e. a series will be loaded as 3D+t

 \b Limitations
  - the 3D+t assumption only works if all time-steps have an equal number of slices and if all
    have the Acquisition Time attribute set to meaningful values
  - Images from tilted CT gantries CAN ONLY be loaded as a series of single-slice images, since
    mitk::Image or the accompanying mapper are not (yet?) capable of representing such geometries
  - Secondary Capture images are expected to have the (0018,2010) tag describing the pixel spacing.
    If only the (0028,0030) tag is set, the spacing will be misinterpreted as (1,1)

 \section DicomSeriesReader_usage Usage

 The starting point for an application is a set of DICOM files that should be loaded.
 For convenience, DicomSeriesReader can also parse a whole directory for DICOM files,
 but an application should better know exactly what to load.

 Loading is then done in two steps:

  1. <b>Group the files into spatial blocks</b> by calling GetSeries(). 
     This method will sort all passed files into meaningful blocks that
     could fit into an mitk::Image. Sorting for 3D+t loading is optional but default.
     The \b return value of this function is a list of identifiers similar to
     DICOM UIDs, each associated to a sorted list of file names.

  2. <b>Load a sorted set of files</b> by calling LoadDicomSeries().
     This method expects go receive the sorting output of GetSeries().
     The method will then invoke ITK methods to actually load the
     files into memory and put them into mitk::Images. Again, loading
     as 3D+t is optional.
 
  Example:

\code

 // only a directory is known at this point: /home/who/dicom

 DicomSeriesReader::UidFileNamesMap allImageBlocks = DicomSeriesReader::GetSeries("/home/who/dicom/");

 // file now divided into groups of identical image size, orientation, spacing, etc.
 // each of these lists should be loadable as an mitk::Image.

 DicomSeriesReader::StringContainer seriesToLoad = allImageBlocks[...]; // decide what to load

 // final step: load into DataNode (can result in 3D+t image)
 DataNode::Pointer node = DicomSeriesReader::LoadDicomSeries( oneBlockSorted );

 Image::Pointer image = dynamic_cast<mitk::Image*>( node->GetData() );
\endcode

 \section DicomSeriesReader_sorting Logic for sorting 2D slices from DICOM images into 3D+t blocks for mitk::Image

 The general sorting mechanism (implemented in GetSeries) groups and sorts a set of DICOM files, each assumed to contain a single CT/MR slice.
 In the following we refer to those file groups as "blocks", since this is what they are meant to become when loaded into an mitk::Image.
 
 \subsection DicomSeriesReader_sorting1 Step 1: Avoiding pure non-sense

 A first pass separates slices that cannot possibly be loaded together because of restrictions of mitk::Image.
 After this steps, each block contains only slices that match in all of the following DICOM tags:

   - (0020,0037) Image Orientation
   - (0028,0030) Pixel Spacing
   - (0018,0050) Slice Thickness
   - (0028,0010) Number Of Rows
   - (0028,0011) Number Of Columns
   - (0020,000e) Series Instance UID : could be argued about, might be dropped in the future (optionally)
 
 \subsection DicomSeriesReader_sorting2 Step 2: Sort slices spatially

 Before slices are further analyzed, they are sorted spatially. As implemented by GdcmSortFunction(),
 slices are sorted by
   1. distance from origin (calculated using (0020,0032) Image Position Patient and (0020,0037) Image Orientation)
   2. when distance is equal, (0008,0032) Acquisition Time is used as a backup criterion (necessary for meaningful 3D+t sorting)

 \subsection DicomSeriesReader_sorting3 Step 3: Ensure equal z spacing

 Since inter-slice distance is not recorded in DICOM tags, we must ensure that blocks are made up of
 slices that have equal distances between neighboring slices. This is especially necessary because itk::ImageSeriesReader
 is later used for the actual loading, and this class expects (and does nocht verify) equal inter-slice distance.

 To achieve such grouping, the inter-slice distance is calculated from the first two different slice positions of a block. 
 Following slices are added to a block as long as they can be added by adding the calculated inter-slice distance to the
 last slice of the block. Slices that do not fit into the expected distance pattern, are set aside for further analysis.
 This grouping is done until each file has been assigned to a group.

 Slices that share a position in space are also sorted into separate blocks during this step. 
 So the result of this step is a set of blocks that contain only slices with equal z spacing 
 and uniqe slices at each position.

 \subsection DicomSeriesReader_sorting4 Step 4 (optional): group 3D blocks as 3D+t when possible

 This last step depends on an option of GetSeries(). When requested, image blocks from the previous step are merged again
 whenever two blocks occupy the same portion of space (i.e. same origin, number of slices and z-spacing).
 
 \section DicomSeriesReader_gantrytilt Handling of gantry tilt

 When CT gantry tilt is used, the gantry plane (= X-Ray source and detector ring) and the vertical plane do not align
 anymore. This scanner feature is used for example to reduce metal artifacs (e.g. <i>Lee C , Evaluation of Using CT
 Gantry Tilt Scan on Head and Neck Cancer Patients with Dental Structure: Scans Show Less Metal Artifacts. Presented
 at: Radiological Society of North America 2011 Scientific Assembly and Annual Meeting; November 27- December 2,
 2011 Chicago IL.</i>).

 The acquired planes of such CT series do not match the expectations of a orthogonal geometry in mitk::Image: if you
 stack the slices, they show a small shift along the Y axis:
\verbatim

  without tilt       with tilt
  
    ||||||             //////
    ||||||            //////
--  |||||| --------- ////// -------- table orientation
    ||||||          //////
    ||||||         //////

Stacked slices:

  without tilt       with tilt

 --------------    --------------
 --------------     --------------
 --------------      --------------
 --------------       --------------
 --------------        --------------

\endverbatim


 As such gemetries do not in conjunction with mitk::Image, DicomSeriesReader performs a correction for such series
 if the groupImagesWithGantryTilt or correctGantryTilt flag in GetSeries and LoadDicomSeries is set (default = on).

 The correction algorithms undoes two errors introduced by ITK's ImageSeriesReader:
  - the plane shift that is ignored by ITK's reader is recreated by applying a shearing transformation using itk::ResampleFilter.
  - the spacing is corrected (it is calculated by ITK's reader from the distance between two origins, which is NOT the slice distance in this special case)

 Both errors are introduced in 
 itkImageSeriesReader.txx (ImageSeriesReader<TOutputImage>::GenerateOutputInformation(void)), lines 176 to 245 (as of ITK 3.20)

 For the correction, we examine two slices of a series, both described as a pair (origin/orientation):
  - we calculate if the second origin is on a line along the normal of the first slice
    - if this is not the case, the geometry will not fit a normal mitk::Image/mitk::Geometry3D
    - we then project the second origin into the first slice's coordinate system to quantify the shift
    - both is done in class GantryTiltInformation with quite some comments.
 
 \section DicomSeriesReader_whynotinitk Why is this not in ITK?
 
  Some of this code would probably be better located in ITK. It is just a matter of resources that this is not the
  case yet. Any attempts into this direction are welcome and can be supported. At least the gantry tilt correction
  should be a simple addition to itk::ImageSeriesReader.
 
 \section DicomSeriesReader_tests Tests regarding DICOM loading

 A number of tests have been implemented to check our assumptions regarding DICOM loading. Please see \ref DICOMTesting

 \todo refactor all the protected helper objects/methods into a separate header so we compile faster
*/
class MITK_CORE_EXPORT DicomSeriesReader
{
public:

  /** 
    \brief Lists of filenames.
  */
  typedef std::vector<std::string> StringContainer;

  /** 
    \brief For grouped lists of filenames, assigned an ID each.
  */
  typedef std::map<std::string, StringContainer> UidFileNamesMap;

  /** 
    \brief Interface for the progress callback.
  */
  typedef void (*UpdateCallBackMethod)(float);

  /**
    \brief Provide combination of preprocessor defines that was active during compilation.

    Since this class is a combination of several possible implementations, separated only
    by ifdef's, calling instances might want to know which flags were active at compile time.
  */
  static std::string GetConfigurationString();

  /**
   \brief Checks if a specific file contains DICOM data.
  */
  static 
  bool 
  IsDicom(const std::string &filename);

  /**
   \brief see other GetSeries().
   
   Find all series (and sub-series -- see details) in a particular directory.
  */
  static UidFileNamesMap GetSeries(const std::string &dir, 
                                   bool groupImagesWithGantryTilt,
                                   const StringContainer &restrictions = StringContainer());

  /**
   \brief see other GetSeries().

   \warning Untested, could or could not work.

   This differs only by having an additional restriction to a single known DICOM series.
   Internally, it uses the other GetSeries() method. 
  */
  static StringContainer GetSeries(const std::string &dir, 
                                   const std::string &series_uid,
                                   bool groupImagesWithGantryTilt,
                                   const StringContainer &restrictions = StringContainer());

  /**
   \brief PREFERRED version of this method - scan and sort DICOM files.

   Parse a list of files for images of DICOM series.
   For each series, an enumeration of the files contained in it is created.
  
   \return The resulting maps UID-like keys (based on Series Instance UID and slice properties) to sorted lists of file names.
  
   SeriesInstanceUID will be enhanced to be unique for each set of file names
   that is later loadable as a single mitk::Image. This implies that
   Image orientation, slice thickness, pixel spacing, rows, and columns
   must be the same for each file (i.e. the image slice contained in the file).
  
   If this separation logic requires that a SeriesInstanceUID must be made more specialized,
   it will follow the same logic as itk::GDCMSeriesFileNames to enhance the UID with
   more digits and dots.
   
   Optionally, more tags can be used to separate files into different logical series by setting
   the restrictions parameter.

   \warning Adding restrictions is not yet implemented!
   */
  static
  UidFileNamesMap 
  GetSeries(const StringContainer& files, 
            bool sortTo3DPlust, 
            bool groupImagesWithGantryTilt,
            const StringContainer &restrictions = StringContainer());
  
  /**
    \brief See other GetSeries().

    Use GetSeries(const StringContainer& files, bool sortTo3DPlust, const StringContainer &restrictions) instead.
  */
  static
  UidFileNamesMap 
  GetSeries(const StringContainer& files, 
            bool groupImagesWithGantryTilt,
            const StringContainer &restrictions = StringContainer());

  /**
   Loads a DICOM series composed by the file names enumerated in the file names container.
   If a callback method is supplied, it will be called after every progress update with a progress value in [0,1].

   \param filenames The filenames to load.
   \param sort Whether files should be sorted spatially (true) or not (false - maybe useful if presorted)
   \param load4D Whether to load the files as 3D+t (if possible)
  */
  static DataNode::Pointer LoadDicomSeries(const StringContainer &filenames, 
                                           bool sort = true, 
                                           bool load4D = true, 
                                           bool correctGantryTilt = true,
                                           UpdateCallBackMethod callback = 0);

  /**
    \brief See LoadDicomSeries! Just a slightly different interface.
  */
  static bool LoadDicomSeries(const StringContainer &filenames, 
                              DataNode &node, 
                              bool sort = true, 
                              bool load4D = true, 
                              bool correctGantryTilt = true,
                              UpdateCallBackMethod callback = 0);

protected:

  /**
    \brief Return type of DicomSeriesReader::AnalyzeFileForITKImageSeriesReaderSpacingAssumption.

    Class contains the grouping result of method DicomSeriesReader::AnalyzeFileForITKImageSeriesReaderSpacingAssumption,
    which takes as input a number of images, which are all equally oriented and spatially sorted along their normal direction.

    The result contains of two blocks: a first one is the grouping result, all of those images can be loaded
    into one image block because they have an equal origin-to-origin distance without any gaps in-between. 
  */
  class SliceGroupingAnalysisResult
  {
    public:

      SliceGroupingAnalysisResult();

      /**
        \brief Grouping result, all same origin-to-origin distance w/o gaps.
      */
      StringContainer GetBlockFilenames();
      
      /**
        \brief Remaining files, which could not be grouped.
      */
      StringContainer GetUnsortedFilenames();
  
      /**
        \brief Wheter or not the grouped result contain a gantry tilt.
      */
      bool ContainsGantryTilt();

      /**
        \brief Meant for internal use by AnalyzeFileForITKImageSeriesReaderSpacingAssumption only.
      */
      void AddFileToSortedBlock(const std::string& filename);
     
      /**
        \brief Meant for internal use by AnalyzeFileForITKImageSeriesReaderSpacingAssumption only.
      */
      void AddFileToUnsortedBlock(const std::string& filename);
      
      /**
        \brief Meant for internal use by AnalyzeFileForITKImageSeriesReaderSpacingAssumption only.
        \todo Could make sense to enhance this with an instance of GantryTiltInformation to store the whole result!
      */
      void FlagGantryTilt();

      /**
        \brief Only meaningful for use by AnalyzeFileForITKImageSeriesReaderSpacingAssumption.
      */
      void UndoPrematureGrouping();

    protected:
      
      StringContainer m_GroupedFiles;
      StringContainer m_UnsortedFiles;

      bool m_GantryTilt;
  };

  /**
    \brief Gantry tilt analysis result.

    Takes geometry information for two slices of a DICOM series and
    calculates whether these fit into an orthogonal block or not.
    If NOT, they can either be the result of an acquisition with
    gantry tilt OR completly broken by some shearing transformation.

    All calculations are done in the constructor, results can then
    be read via the remaining methods.
  */
  class GantryTiltInformation
  {
    public:

      /**
        \brief Just so we can create empty instances for assigning results later.
      */
      GantryTiltInformation();

      /**
        \brief THE constructor, which does all the calculations.

        See code comments for explanation.
      */
      GantryTiltInformation( const Point3D& origin1, 
                             const Point3D& origin2,
                             const Vector3D& right, 
                             const Vector3D& up,
                             unsigned int numberOfSlicesApart);

      /**
        \brief Whether the slices were sheared.
      */
      bool IsSheared() const;

      /**
        \brief Whether the shearing is a gantry tilt or more complicated.
      */
      bool IsRegularGantryTilt() const;

      /**
        \brief The offset distance in Y direction for each slice (describes the tilt result).
      */
      ScalarType GetMatrixCoefficientForCorrectionInWorldCoordinates() const;


      /**
        \brief The z / inter-slice spacing. Needed to correct ImageSeriesReader's result.
      */
      ScalarType GetRealZSpacing() const;

      /**
        \brief The shift between first and last slice in mm.

        Needed to resize an orthogonal image volume.
      */
      ScalarType GetTiltCorrectedAdditionalSize() const;

    protected:

      /**
        \brief Projection of point p onto line through lineOrigin in direction of lineDirection.
      */
      Point3D projectPointOnLine( Point3D p, Point3D lineOrigin, Vector3D lineDirection ); 

      ScalarType m_ShiftUp;
      ScalarType m_ShiftRight;
      ScalarType m_ShiftNormal;
      unsigned int m_NumberOfSlicesApart;
  };

  /**
    \brief for internal sorting.
  */
  typedef std::pair<StringContainer, StringContainer> TwoStringContainers;
 
  /**
    \brief Maps DICOM tags to MITK properties.
  */
  typedef std::map<std::string, std::string> TagToPropertyMapType;

  /**
    \brief Ensure an equal z-spacing for a group of files.
    
    Takes as input a number of images, which are all equally oriented and spatially sorted along their normal direction.

    Internally used by GetSeries. Returns two lists: the first one contins slices of equal inter-slice spacing.
    The second list contains remaining files, which need to be run through AnalyzeFileForITKImageSeriesReaderSpacingAssumption again.
    
    Relevant code that is matched here is in
    itkImageSeriesReader.txx (ImageSeriesReader<TOutputImage>::GenerateOutputInformation(void)), lines 176 to 245 (as of ITK 3.20)
  */
  static
  SliceGroupingAnalysisResult
  AnalyzeFileForITKImageSeriesReaderSpacingAssumption(const StringContainer& files, bool groupsOfSimilarImages, const gdcm::Scanner::MappingType& tagValueMappings_);
      
  static
  std::string
  ConstCharStarToString(const char* s);

  static
  Point3D
  DICOMStringToPoint3D(const std::string& s, bool& successful);

  static
  void
  DICOMStringToOrientationVectors(const std::string& s, Vector3D& right, Vector3D& up, bool& successful);

  template <typename ImageType>
  static
  typename ImageType::Pointer
  InPlaceFixUpTiltedGeometry( ImageType* input, const GantryTiltInformation& tiltInfo );


  /**
    \brief Sort a set of file names in an order that is meaningful for loading them into an mitk::Image.
   
    \warning This method assumes that input files are similar in basic properties such as 
             slice thicknes, image orientation, pixel spacing, rows, columns. 
             It should always be ok to put the result of a call to GetSeries(..) into this method.
   
    Sorting order is determined by
   
     1. image position along its normal (distance from world origin)
     2. acquisition time
   
    If P<n> denotes a position and T<n> denotes a time step, this method will order slices from three timesteps like this:
\verbatim
  P1T1 P1T2 P1T3 P2T1 P2T2 P2T3 P3T1 P3T2 P3T3
\endverbatim
   
   */
  static StringContainer SortSeriesSlices(const StringContainer &unsortedFilenames);

public:
  /**
   \brief Checks if a specific file is a Philips3D ultrasound DICOM file.
  */
  static bool IsPhilips3DDicom(const std::string &filename);
protected:

  /**
   \brief Read a Philips3D ultrasound DICOM file and put into an mitk::Image.
  */
  static bool ReadPhilips3DDicom(const std::string &filename, mitk::Image::Pointer output_image);
     
  /**
    \brief Construct a UID that takes into account sorting criteria from GetSeries().
  */
  static std::string CreateMoreUniqueSeriesIdentifier( gdcm::Scanner::TagToValue& tagValueMap );
  
  /**
    \brief Helper for CreateMoreUniqueSeriesIdentifier
  */
  static std::string CreateSeriesIdentifierPart( gdcm::Scanner::TagToValue& tagValueMap, const gdcm::Tag& tag );
  
  /**
    \brief Helper for CreateMoreUniqueSeriesIdentifier
  */
  static std::string IDifyTagValue(const std::string& value);

  typedef itk::GDCMImageIO DcmIoType;

  /**
    \brief Progress callback for DicomSeriesReader.
  */
  class CallbackCommand : public itk::Command
  {
  public:
    CallbackCommand(UpdateCallBackMethod callback)
      : m_Callback(callback)
    {
    }

    void Execute(const itk::Object *caller, const itk::EventObject&)
    {
      (*this->m_Callback)(static_cast<const itk::ProcessObject*>(caller)->GetProgress());
    }

    void Execute(itk::Object *caller, const itk::EventObject&)
    {
      (*this->m_Callback)(static_cast<itk::ProcessObject*>(caller)->GetProgress());
    }

  protected:

    UpdateCallBackMethod m_Callback;
  };

  /**
   \brief Scan for slice image information
  */
  static void ScanForSliceInformation( const StringContainer &filenames, gdcm::Scanner& scanner );

  /**
   \brief Performs actual loading of a series and creates an image having the specified pixel type.
  */
  template <typename PixelType>
  static 
  void 
  LoadDicom(const StringContainer &filenames, DataNode &node, bool sort, bool check_4d, bool correctTilt, UpdateCallBackMethod callback);

  /**
    \brief Feed files into itk::ImageSeriesReader and retrieve a 3D MITK image.

    \param command can be used for progress reporting
  */
  template <typename PixelType>
  static 
  Image::Pointer 
  LoadDICOMByITK( const StringContainer&, bool correctTilt, const GantryTiltInformation& tiltInfo, CallbackCommand* command = NULL);

  /**
    \brief Sort files into time step blocks of a 3D+t image.

    Called by LoadDicom. Expects to be fed a single list of filenames that have been sorted by
    GetSeries previously (one map entry). This method will check how many timestep can be filled
    with given files.

    Assumption is that the number of time steps is determined by how often the first position in
    space repeats. I.e. if the first three files in the input parameter all describe the same
    location in space, we'll construct three lists of files. and sort the remaining files into them.

    \todo We can probably remove this method if we somehow transfer 3D+t information from GetSeries to LoadDicomSeries.
  */
  static 
  std::list<StringContainer> 
  SortIntoBlocksFor3DplusT( const StringContainer& presortedFilenames, const gdcm::Scanner::MappingType& tagValueMappings_, bool sort, bool& canLoadAs4D);

  /**
   \brief Defines spatial sorting for sorting by GDCM 2.

   Sorts by image position along image normal (distance from world origin).
   In cases of conflict, acquisition time is used as a secondary sort criterium.
  */
  static 
  bool 
  GdcmSortFunction(const gdcm::DataSet &ds1, const gdcm::DataSet &ds2);


  /**
    \brief Copy information about files and DICOM tags from ITK's MetaDataDictionary
           and from the list of input files to the PropertyList of mitk::Image.
    \todo Tag copy must follow; image level will cause some additional files parsing, probably.
  */
  static void CopyMetaDataToImageProperties( StringContainer filenames, const gdcm::Scanner::MappingType& tagValueMappings_, DcmIoType* io, Image* image);
  static void CopyMetaDataToImageProperties( std::list<StringContainer> imageBlock, const gdcm::Scanner::MappingType& tagValueMappings_, DcmIoType* io, Image* image);
  
  /**
    \brief Map between DICOM tags and MITK properties.

    Uses as a positive list for copying specified DICOM tags (from ITK's ImageIO)
    to MITK properties. ITK provides MetaDataDictionary entries of form
    "gggg|eeee" (g = group, e = element), e.g. "0028,0109" (Largest Pixel in Series),
    which we want to sort as dicom.series.largest_pixel_in_series".
  */
  static const TagToPropertyMapType& GetDICOMTagsToMITKPropertyMap();
};

}


#endif /* MITKDICOMSERIESREADER_H_ */
