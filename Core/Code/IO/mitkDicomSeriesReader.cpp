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

// uncomment for learning more about the internal sorting mechanisms
//#define MBILOG_ENABLE_DEBUG 

#include <mitkDicomSeriesReader.h>

#include <itkGDCMSeriesFileNames.h>

#include <gdcmAttribute.h>
#include <gdcmPixmapReader.h>
#include <gdcmStringFilter.h>
#include <gdcmDirectory.h>
#include <gdcmScanner.h>

#include "mitkProperties.h"

namespace mitk
{

typedef itk::GDCMSeriesFileNames DcmFileNamesGeneratorType;


DicomSeriesReader::SliceGroupingAnalysisResult::SliceGroupingAnalysisResult()
:m_GantryTilt(false)
{
}

DicomSeriesReader::StringContainer DicomSeriesReader::SliceGroupingAnalysisResult::GetBlockFilenames()
{
  return m_GroupedFiles;
}

DicomSeriesReader::StringContainer DicomSeriesReader::SliceGroupingAnalysisResult::GetUnsortedFilenames()
{
  return m_UnsortedFiles;
}

bool DicomSeriesReader::SliceGroupingAnalysisResult::ContainsGantryTilt()
{
  return m_GantryTilt;
}

void DicomSeriesReader::SliceGroupingAnalysisResult::AddFileToSortedBlock(const std::string& filename)
{
  m_GroupedFiles.push_back( filename );
}

void DicomSeriesReader::SliceGroupingAnalysisResult::AddFileToUnsortedBlock(const std::string& filename)
{
  m_UnsortedFiles.push_back( filename );
}

void DicomSeriesReader::SliceGroupingAnalysisResult::FlagGantryTilt()
{
  m_GantryTilt = true;
}

void DicomSeriesReader::SliceGroupingAnalysisResult::UndoPrematureGrouping()
{
  assert( !m_GroupedFiles.empty() );
  m_UnsortedFiles.insert( m_UnsortedFiles.begin(), m_GroupedFiles.back() );
  m_GroupedFiles.pop_back();
}







const DicomSeriesReader::TagToPropertyMapType& DicomSeriesReader::GetDICOMTagsToMITKPropertyMap()
{
  static bool initialized = false;
  static TagToPropertyMapType dictionary;
  if (!initialized)
  {
    /*
       Selection criteria:
       - no sequences because we cannot represent that
       - nothing animal related (specied, breed registration number), MITK focusses on human medical image processing.
       - only general attributes so far

       When extending this, we should make use of a real dictionary (GDCM/DCMTK and lookup the tag names there)
    */

    // Patient module
    dictionary["0010|0010"] = "dicom.patient.PatientsName";
    dictionary["0010|0020"] = "dicom.patient.PatientID";
    dictionary["0010|0030"] = "dicom.patient.PatientsBirthDate";
    dictionary["0010|0040"] = "dicom.patient.PatientsSex";
    dictionary["0010|0032"] = "dicom.patient.PatientsBirthTime";
    dictionary["0010|1000"] = "dicom.patient.OtherPatientIDs";
    dictionary["0010|1001"] = "dicom.patient.OtherPatientNames";
    dictionary["0010|2160"] = "dicom.patient.EthnicGroup";
    dictionary["0010|4000"] = "dicom.patient.PatientComments";
    dictionary["0012|0062"] = "dicom.patient.PatientIdentityRemoved";
    dictionary["0012|0063"] = "dicom.patient.DeIdentificationMethod";

    // General Study module
    dictionary["0020|000d"] = "dicom.study.StudyInstanceUID";
    dictionary["0008|0020"] = "dicom.study.StudyDate";
    dictionary["0008|0030"] = "dicom.study.StudyTime";
    dictionary["0008|0090"] = "dicom.study.ReferringPhysiciansName";
    dictionary["0020|0010"] = "dicom.study.StudyID";
    dictionary["0008|0050"] = "dicom.study.AccessionNumber";
    dictionary["0008|1030"] = "dicom.study.StudyDescription";
    dictionary["0008|1048"] = "dicom.study.PhysiciansOfRecord";
    dictionary["0008|1060"] = "dicom.study.NameOfPhysicianReadingStudy";

    // General Series module
    dictionary["0008|0060"] = "dicom.series.Modality";
    dictionary["0020|000e"] = "dicom.series.SeriesInstanceUID";
    dictionary["0020|0011"] = "dicom.series.SeriesNumber";
    dictionary["0020|0060"] = "dicom.series.Laterality";
    dictionary["0008|0021"] = "dicom.series.SeriesDate";
    dictionary["0008|0031"] = "dicom.series.SeriesTime";
    dictionary["0008|1050"] = "dicom.series.PerformingPhysiciansName";
    dictionary["0018|1030"] = "dicom.series.ProtocolName";
    dictionary["0008|103e"] = "dicom.series.SeriesDescription";
    dictionary["0008|1070"] = "dicom.series.OperatorsName";
    dictionary["0018|0015"] = "dicom.series.BodyPartExamined";
    dictionary["0018|5100"] = "dicom.series.PatientPosition";
    dictionary["0028|0108"] = "dicom.series.SmallestPixelValueInSeries";
    dictionary["0028|0109"] = "dicom.series.LargestPixelValueInSeries";
   
    // VOI LUT module
    dictionary["0028|1050"] = "dicom.voilut.WindowCenter";
    dictionary["0028|1051"] = "dicom.voilut.WindowWidth";
    dictionary["0028|1055"] = "dicom.voilut.WindowCenterAndWidthExplanation";

    initialized = true;
  }

  return dictionary;
}


DataNode::Pointer 
DicomSeriesReader::LoadDicomSeries(const StringContainer &filenames, bool sort, bool check_4d, bool correctTilt, UpdateCallBackMethod callback)
{
  DataNode::Pointer node = DataNode::New();

  if (DicomSeriesReader::LoadDicomSeries(filenames, *node, sort, check_4d, correctTilt, callback))
  {
    if( filenames.empty() )
    {
      return NULL;
    }

    return node;
  }
  else
  {
    return NULL;
  }
}

bool 
DicomSeriesReader::LoadDicomSeries(const StringContainer &filenames, DataNode &node, bool sort, bool check_4d, bool correctTilt, UpdateCallBackMethod callback)
{
  if( filenames.empty() )
  {
    MITK_WARN << "Calling LoadDicomSeries with empty filename string container. Probably invalid application logic.";
    node.SetData(NULL);
    return true; // this is not actually an error but the result is very simple
  }

  DcmIoType::Pointer io = DcmIoType::New();

  try
  {
    if (io->CanReadFile(filenames.front().c_str()))
    {
      io->SetFileName(filenames.front().c_str());
      io->ReadImageInformation();

      switch (io->GetComponentType())
      {
      case DcmIoType::UCHAR:
        DicomSeriesReader::LoadDicom<unsigned char>(filenames, node, sort, check_4d, correctTilt, callback);
        break;
      case DcmIoType::CHAR:
        DicomSeriesReader::LoadDicom<char>(filenames, node, sort, check_4d, correctTilt, callback);
        break;
      case DcmIoType::USHORT:
        DicomSeriesReader::LoadDicom<unsigned short>(filenames, node, sort, check_4d, correctTilt, callback);
        break;
      case DcmIoType::SHORT:
        DicomSeriesReader::LoadDicom<short>(filenames, node, sort, check_4d, correctTilt, callback);
        break;
      case DcmIoType::UINT:
        DicomSeriesReader::LoadDicom<unsigned int>(filenames, node, sort, check_4d, correctTilt, callback);
        break;
      case DcmIoType::INT:
        DicomSeriesReader::LoadDicom<int>(filenames, node, sort, check_4d, correctTilt, callback);
        break;
      case DcmIoType::ULONG:
        DicomSeriesReader::LoadDicom<long unsigned int>(filenames, node, sort, check_4d, correctTilt, callback);
        break;
      case DcmIoType::LONG:
        DicomSeriesReader::LoadDicom<long int>(filenames, node, sort, check_4d, correctTilt, callback);
        break;
      case DcmIoType::FLOAT:
        DicomSeriesReader::LoadDicom<float>(filenames, node, sort, check_4d, correctTilt, callback);
        break;
      case DcmIoType::DOUBLE:
        DicomSeriesReader::LoadDicom<double>(filenames, node, sort, check_4d, correctTilt, callback);
        break;
      default:
        MITK_ERROR << "Found unsupported DICOM pixel type: (enum value) " << io->GetComponentType();
      }
      
      if (node.GetData())
      {
        return true;
      }
    }
  }
  catch(itk::MemoryAllocationError& e)
  {
    MITK_ERROR << "Out of memory. Cannot load DICOM series: " << e.what();
  }
  catch(std::exception& e)
  {
    MITK_ERROR << "Error encountered when loading DICOM series:" << e.what();
  }
  catch(...)
  {
    MITK_ERROR << "Unspecified error encountered when loading DICOM series.";
  }

  return false;
}


bool 
DicomSeriesReader::IsDicom(const std::string &filename)
{
  DcmIoType::Pointer io = DcmIoType::New();

  return io->CanReadFile(filename.c_str());
}

bool 
DicomSeriesReader::IsPhilips3DDicom(const std::string &filename)
{
  DcmIoType::Pointer io = DcmIoType::New();

  if (io->CanReadFile(filename.c_str()))
  {
    //Look at header Tag 3001,0010 if it is "Philips3D"
    gdcm::Reader reader;
    reader.SetFileName(filename.c_str());
    reader.Read();
    gdcm::DataSet &data_set = reader.GetFile().GetDataSet();
    gdcm::StringFilter sf;
    sf.SetFile(reader.GetFile());

    if (data_set.FindDataElement(gdcm::Tag(0x3001, 0x0010)) &&
      (sf.ToString(gdcm::Tag(0x3001, 0x0010)) == "Philips3D "))
    {
      return true;
    }
  }
  return false;
}

bool 
DicomSeriesReader::ReadPhilips3DDicom(const std::string &filename, mitk::Image::Pointer output_image)
{
  // Now get PhilipsSpecific Tags

  gdcm::PixmapReader reader;
  reader.SetFileName(filename.c_str());
  reader.Read();
  gdcm::DataSet &data_set = reader.GetFile().GetDataSet();
  gdcm::StringFilter sf;
  sf.SetFile(reader.GetFile());

  gdcm::Attribute<0x0028,0x0011> dimTagX; // coloumns || sagittal
  gdcm::Attribute<0x3001,0x1001, gdcm::VR::UL, gdcm::VM::VM1> dimTagZ; //I have no idea what is VM1. // (Philips specific) // transversal
  gdcm::Attribute<0x0028,0x0010> dimTagY; // rows || coronal
  gdcm::Attribute<0x0028,0x0008> dimTagT; // how many frames
  gdcm::Attribute<0x0018,0x602c> spaceTagX; // Spacing in X , unit is "physicalTagx" (usually centimeter)
  gdcm::Attribute<0x0018,0x602e> spaceTagY;
  gdcm::Attribute<0x3001,0x1003, gdcm::VR::FD, gdcm::VM::VM1> spaceTagZ; // (Philips specific)
  gdcm::Attribute<0x0018,0x6024> physicalTagX; // if 3, then spacing params are centimeter
  gdcm::Attribute<0x0018,0x6026> physicalTagY;
  gdcm::Attribute<0x3001,0x1002, gdcm::VR::US, gdcm::VM::VM1> physicalTagZ; // (Philips specific)

  dimTagX.Set(data_set);
  dimTagY.Set(data_set);
  dimTagZ.Set(data_set);
  dimTagT.Set(data_set);
  spaceTagX.Set(data_set);
  spaceTagY.Set(data_set);
  spaceTagZ.Set(data_set);
  physicalTagX.Set(data_set);
  physicalTagY.Set(data_set);
  physicalTagZ.Set(data_set);

  unsigned int
    dimX = dimTagX.GetValue(),
    dimY = dimTagY.GetValue(),
    dimZ = dimTagZ.GetValue(),
    dimT = dimTagT.GetValue(),
    physicalX = physicalTagX.GetValue(),
    physicalY = physicalTagY.GetValue(),
    physicalZ = physicalTagZ.GetValue();

  float
    spaceX = spaceTagX.GetValue(),
    spaceY = spaceTagY.GetValue(),
    spaceZ = spaceTagZ.GetValue();

  if (physicalX == 3) // spacing parameter in cm, have to convert it to mm.
    spaceX = spaceX * 10;

  if (physicalY == 3) // spacing parameter in cm, have to convert it to mm.
    spaceY = spaceY * 10;
  
  if (physicalZ == 3) // spacing parameter in cm, have to convert it to mm.
    spaceZ = spaceZ * 10;

  // Ok, got all necessary Tags!
  // Now read Pixeldata (7fe0,0010) X x Y x Z x T Elements

  const gdcm::Pixmap &pixels = reader.GetPixmap();
  gdcm::RAWCodec codec;

  codec.SetPhotometricInterpretation(gdcm::PhotometricInterpretation::MONOCHROME2);
  codec.SetPixelFormat(pixels.GetPixelFormat());
  codec.SetPlanarConfiguration(0);

  gdcm::DataElement out;
  codec.Decode(data_set.GetDataElement(gdcm::Tag(0x7fe0, 0x0010)), out);

  const gdcm::ByteValue *bv = out.GetByteValue();
  const char *new_pixels = bv->GetPointer();

  // Create MITK Image + Geometry
  typedef itk::Image<unsigned char, 4> ImageType;   //Pixeltype might be different sometimes? Maybe read it out from header
  ImageType::RegionType myRegion;
  ImageType::SizeType mySize;
  ImageType::IndexType myIndex;
  ImageType::SpacingType mySpacing;
  ImageType::Pointer imageItk = ImageType::New();

  mySpacing[0] = spaceX;
  mySpacing[1] = spaceY;
  mySpacing[2] = spaceZ;
  mySpacing[3] = 1;
  myIndex[0] = 0;
  myIndex[1] = 0;
  myIndex[2] = 0;
  myIndex[3] = 0;
  mySize[0] = dimX;
  mySize[1] = dimY;
  mySize[2] = dimZ;
  mySize[3] = dimT;
  myRegion.SetSize( mySize);
  myRegion.SetIndex( myIndex );
  imageItk->SetSpacing(mySpacing);
  imageItk->SetRegions( myRegion);
  imageItk->Allocate();
  imageItk->FillBuffer(0);

  itk::ImageRegionIterator<ImageType>  iterator(imageItk, imageItk->GetLargestPossibleRegion());
  iterator.GoToBegin();
  unsigned long pixCount = 0;
  unsigned long planeSize = dimX*dimY;
  unsigned long planeCount = 0;
  unsigned long timeCount = 0;
  unsigned long numberOfSlices = dimZ;

  while (!iterator.IsAtEnd())
  {
    unsigned long adressedPixel =
      pixCount
      + (numberOfSlices-1-planeCount)*planeSize // add offset to adress the first pixel of current plane
      + timeCount*numberOfSlices*planeSize; // add time offset

    iterator.Set( new_pixels[ adressedPixel ] );
    pixCount++;
    ++iterator;

    if (pixCount == planeSize)
    {
      pixCount = 0;
      planeCount++;
    }
    if (planeCount == numberOfSlices)
    {
      planeCount = 0;
      timeCount++;
    }
    if (timeCount == dimT)
    {
      break;
    }
  }
  mitk::CastToMitkImage(imageItk, output_image);
  return true; // actually never returns false yet.. but exception possible
}
      
DicomSeriesReader::GantryTiltInformation::GantryTiltInformation()
: m_ShiftUp(0.0)
, m_ShiftRight(0.0)
, m_ShiftNormal(0.0)
, m_NumberOfSlicesApart(1)
{
}

DicomSeriesReader::GantryTiltInformation::GantryTiltInformation( 
    const Point3D& origin1, const Point3D& origin2,
    const Vector3D& right, const Vector3D& up,
    unsigned int numberOfSlicesApart)
: m_ShiftUp(0.0)
, m_ShiftRight(0.0)
, m_ShiftNormal(0.0)
, m_NumberOfSlicesApart(numberOfSlicesApart)
{
  assert(numberOfSlicesApart);
  // determine if slice 1 (imagePosition1 and imageOrientation1) and slice 2 can be in one orthogonal slice stack:
  // calculate a line from origin 1, directed along the normal of slice (calculated as the cross product of orientation 1)
  // check if this line passes through origin 2
        
  /* 
     Determine if line (imagePosition2 + l * normal) contains imagePosition1.
     Done by calculating the distance of imagePosition1 from line (imagePosition2 + l *normal)

     E.g. http://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html

     squared distance = | (pointAlongNormal - origin2) x (origin2 - origin1) | ^ 2
     /
     |pointAlongNormal - origin2| ^ 2

     ( x meaning the cross product )
  */

  Vector3D normal = itk::CrossProduct(right, up);
  Point3D pointAlongNormal = origin2 + normal;

  double numerator = itk::CrossProduct( pointAlongNormal - origin2 , origin2 - origin1 ).GetSquaredNorm();
  double denominator = (pointAlongNormal - origin2).GetSquaredNorm();

  double distance = sqrt(numerator / denominator);

  if ( distance > 0.001 ) // mitk::eps is too small; 1/1000 of a mm should be enough to detect tilt
  {
    MITK_DEBUG << "  Series seems to contain a tilted (or sheared) geometry";
    MITK_DEBUG << "  Distance of expected slice origin from actual slice origin: " << distance;
    MITK_DEBUG << "    ==> storing this shift for later analysis:";
    MITK_DEBUG << "    v right: " << right;
    MITK_DEBUG << "    v up: " << up;
    MITK_DEBUG << "    v normal: " << normal;

    Point3D projectionRight = projectPointOnLine( origin1, origin2, right );
    Point3D projectionUp = projectPointOnLine( origin1, origin2, up );
    Point3D projectionNormal = projectPointOnLine( origin1, origin2, normal );

    m_ShiftRight = (projectionRight - origin2).GetNorm();
    m_ShiftUp = (projectionUp - origin2).GetNorm();
    m_ShiftNormal = (projectionNormal - origin2).GetNorm();

    MITK_DEBUG << "    shift normal: " << m_ShiftNormal;
    MITK_DEBUG << "    shift up: " << m_ShiftUp;
    MITK_DEBUG << "    shift right: " << m_ShiftRight;
    
    MITK_DEBUG << "    tilt angle (rad): " << tanh( m_ShiftUp / m_ShiftNormal );
    MITK_DEBUG << "    tilt angle (deg): " << tanh( m_ShiftUp / m_ShiftNormal ) * 180.0 / 3.1415926535;
  }
}
      
Point3D 
DicomSeriesReader::GantryTiltInformation::projectPointOnLine( Point3D p, Point3D lineOrigin, Vector3D lineDirection )
{
  /**
    See illustration at http://mo.mathematik.uni-stuttgart.de/inhalt/aussage/aussage472/

    vector(lineOrigin,p) = normal * ( innerproduct((p - lineOrigin),normal) / squared-length(normal) )
  */

  Vector3D lineOriginToP = p - lineOrigin;
  ScalarType innerProduct = lineOriginToP * lineDirection;
  
  ScalarType factor = innerProduct / lineDirection.GetSquaredNorm();
  Point3D projection = lineOrigin + factor * lineDirection;

  return projection;
}

ScalarType 
DicomSeriesReader::GantryTiltInformation::GetTiltCorrectedAdditionalSize() const
{
  // this seems to be a bit too much sometimes, but better too much than cutting off parts of the image
  return int(m_ShiftUp + 1.0); // to next bigger int: plus 1, then cut off after point
}

ScalarType 
DicomSeriesReader::GantryTiltInformation::GetMatrixCoefficientForCorrectionInWorldCoordinates() const
{
  // so many mm shifted per slice!
  return m_ShiftUp / static_cast<ScalarType>(m_NumberOfSlicesApart);
}

ScalarType 
DicomSeriesReader::GantryTiltInformation::GetRealZSpacing() const
{
  return m_ShiftNormal / static_cast<ScalarType>(m_NumberOfSlicesApart);
}


bool 
DicomSeriesReader::GantryTiltInformation::IsSheared() const
{
  return (   m_ShiftRight > 0.001
          ||    m_ShiftUp > 0.001);
}

      
bool 
DicomSeriesReader::GantryTiltInformation::IsRegularGantryTilt() const
{
  return (   m_ShiftRight < 0.001 
          &&    m_ShiftUp > 0.001);
}


std::string
DicomSeriesReader::ConstCharStarToString(const char* s)
{
  return s ?  std::string(s) : std::string();
}
  
Point3D
DicomSeriesReader::DICOMStringToPoint3D(const std::string& s, bool& successful)
{
  Point3D p;
  successful = true;

  std::istringstream originReader(s);
  std::string coordinate;
  unsigned int dim(0);
  while( std::getline( originReader, coordinate, '\\' ) ) 
  {
    if (dim > 4) break; // otherwise we access invalid array index

    p[dim++] = atof(coordinate.c_str());
  }

  if (dim != 3)
  {
    successful = false;
    MITK_ERROR << "Reader implementation made wrong assumption on tag (0020,0032). Found " << dim << " instead of 3 values.";
  }

  return p;
}
  
void
DicomSeriesReader::DICOMStringToOrientationVectors(const std::string& s, Vector3D& right, Vector3D& up, bool& successful)
{
  successful = true;

  std::istringstream orientationReader(s);
  std::string coordinate;
  unsigned int dim(0);
  while( std::getline( orientationReader, coordinate, '\\' ) )
  {
    if (dim > 6) break; // otherwise we access invalid array index

    if (dim<3) 
    {
      right[dim++] = atof(coordinate.c_str());
    }
    else
    {
      up[dim++ - 3] = atof(coordinate.c_str());
    }
  }

  if (dim != 6)
  {
    successful = false;
    MITK_ERROR << "Reader implementation made wrong assumption on tag (0020,0037). Found " << dim << " instead of 6 values.";
  }
}


DicomSeriesReader::SliceGroupingAnalysisResult
DicomSeriesReader::AnalyzeFileForITKImageSeriesReaderSpacingAssumption(
    const StringContainer& files,
    bool groupImagesWithGantryTilt,
    const gdcm::Scanner::MappingType& tagValueMappings_)
{
  // result.first = files that fit ITK's assumption
  // result.second = files that do not fit, should be run through AnalyzeFileForITKImageSeriesReaderSpacingAssumption() again
  SliceGroupingAnalysisResult result;
 
  // we const_cast here, because I could not use a map.at(), which would make the code much more readable
  gdcm::Scanner::MappingType& tagValueMappings = const_cast<gdcm::Scanner::MappingType&>(tagValueMappings_);
  const gdcm::Tag tagImagePositionPatient(0x0020,0x0032); // Image Position (Patient)
  const gdcm::Tag    tagImageOrientation(0x0020, 0x0037); // Image Orientation

  Vector3D fromFirstToSecondOrigin; fromFirstToSecondOrigin.Fill(0.0);
  bool fromFirstToSecondOriginInitialized(false);
  Point3D thisOrigin;
  Point3D lastOrigin;
  lastOrigin.Fill(0.0f);
  Point3D lastDifferentOrigin;
  lastDifferentOrigin.Fill(0.0f);

  bool lastOriginInitialized(false);

  MITK_DEBUG << "--------------------------------------------------------------------------------";
  MITK_DEBUG << "Analyzing files for z-spacing assumption of ITK's ImageSeriesReader (group tilted: " << groupImagesWithGantryTilt << ")";
  unsigned int fileIndex(0);
  for (StringContainer::const_iterator fileIter = files.begin();
       fileIter != files.end();
       ++fileIter, ++fileIndex)
  {
    bool fileFitsIntoPattern(false);
    std::string thisOriginString;
    // Read tag value into point3D. PLEASE replace this by appropriate GDCM code if you figure out how to do that
    thisOriginString = ConstCharStarToString( tagValueMappings[fileIter->c_str()][tagImagePositionPatient] );

    bool ignoredConversionError(-42); // hard to get here, no graceful way to react
    thisOrigin = DICOMStringToPoint3D( thisOriginString, ignoredConversionError );

    MITK_DEBUG << "  " << fileIndex << " " << *fileIter
                       << " at "
                       << thisOriginString << "(" << thisOrigin[0] << "," << thisOrigin[1] << "," << thisOrigin[2] << ")";

    if ( lastOriginInitialized && (thisOrigin == lastOrigin) )
    {
      MITK_DEBUG << "    ==> Sort away " << *fileIter << " for separate time step"; // we already have one occupying this position
      result.AddFileToUnsortedBlock( *fileIter );
      fileFitsIntoPattern = false;
    }
    else
    {
      if (!fromFirstToSecondOriginInitialized && lastOriginInitialized) // calculate vector as soon as possible when we get a new position
      {
        fromFirstToSecondOrigin = thisOrigin - lastDifferentOrigin;
        fromFirstToSecondOriginInitialized = true;

        // Here we calculate if this slice and the previous one are well aligned,
        // i.e. we test if the previous origin is on a line through the current
        // origin, directed into the normal direction of the current slice.

        // If this is NOT the case, then we have a data set with a TILTED GANTRY geometry,
        // which cannot be simply loaded into a single mitk::Image at the moment.
        // For this case, we flag this finding in the result and DicomSeriesReader
        // can correct for that later.

        Vector3D right; right.Fill(0.0);
        Vector3D up; right.Fill(0.0); // might be down as well, but it is just a name at this point
        DICOMStringToOrientationVectors( tagValueMappings[fileIter->c_str()][tagImageOrientation], right, up, ignoredConversionError );
       
        GantryTiltInformation tiltInfo( lastDifferentOrigin, thisOrigin, right, up, 1 );

        if ( tiltInfo.IsSheared() ) // mitk::eps is too small; 1/1000 of a mm should be enough to detect tilt
        {
          /* optimistic approach, accepting gantry tilt: save file for later, check all further files */
          
          // at this point we have TWO slices analyzed! if they are the only two files, we still split, because there is no third to verify our tilting assumption.
          // later with a third being available, we must check if the initial tilting vector is still valid. if yes, continue. 
          // if NO, we need to split the already sorted part (result.first) and the currently analyzed file (*fileIter)

          // tell apart gantry tilt from overall skewedness
          // sort out irregularly sheared slices, that IS NOT tilting

          if ( groupImagesWithGantryTilt && tiltInfo.IsRegularGantryTilt() )
          {
            result.FlagGantryTilt();
            result.AddFileToSortedBlock(*fileIter); // this file is good for current block
            fileFitsIntoPattern = true;
          }
          else
          {
            result.AddFileToUnsortedBlock( *fileIter ); // sort away for further analysis
            fileFitsIntoPattern = false;
          }
        }
        else
        {
          result.AddFileToSortedBlock(*fileIter); // this file is good for current block
          fileFitsIntoPattern = true;
        }
      }
      else if (fromFirstToSecondOriginInitialized) // we already know the offset between slices
      {
        Point3D assumedOrigin = lastDifferentOrigin + fromFirstToSecondOrigin;

        Vector3D originError = assumedOrigin - thisOrigin;
        double norm = originError.GetNorm();
        double toleratedError(0.005); // max. 1/10mm error when measurement crosses 20 slices in z direction

        if (norm > toleratedError)
        {
          MITK_DEBUG << "  File does not fit into the inter-slice distance pattern (diff = " 
                               << norm << ", allowed " 
                               << toleratedError << ").";
          MITK_DEBUG << "  Expected position (" << assumedOrigin[0] << ","
                                            << assumedOrigin[1] << ","
                                            << assumedOrigin[2] << "), got position ("
                                            << thisOrigin[0] << ","
                                            << thisOrigin[1] << ","
                                            << thisOrigin[2] << ")";
          MITK_DEBUG  << "    ==> Sort away " << *fileIter << " for later analysis";

          // At this point we know we deviated from the expectation of ITK's ImageSeriesReader
          // We split the input file list at this point, i.e. all files up to this one (excluding it)
          // are returned as group 1, the remaining files (including the faulty one) are group 2
          
          /* Optimistic approach: check if any of the remaining slices fits in */
          result.AddFileToUnsortedBlock( *fileIter ); // sort away for further analysis
          fileFitsIntoPattern = false;
        }
        else
        {
          result.AddFileToSortedBlock(*fileIter); // this file is good for current block
          fileFitsIntoPattern = true;
        }
      }
      else // this should be the very first slice
      {
        result.AddFileToSortedBlock(*fileIter); // this file is good for current block
        fileFitsIntoPattern = true;
      }
    }

    // record current origin for reference in later iterations
    if ( !lastOriginInitialized || ( fileFitsIntoPattern && (thisOrigin != lastOrigin) ) )
    {
      lastDifferentOrigin = thisOrigin;
    }

    lastOrigin = thisOrigin; 
    lastOriginInitialized = true;
  }

  if ( result.ContainsGantryTilt() )
  {
    // check here how many files were grouped.
    // IF it was only two files AND we assume tiltedness (e.g. save "distance")
    // THEN we would want to also split the two previous files (simple) because
    // we don't have any reason to assume they belong together

    if ( result.GetBlockFilenames().size() == 2 )
    {
      result.UndoPrematureGrouping();
    }
  }
  
  return result;
}

DicomSeriesReader::UidFileNamesMap 
DicomSeriesReader::GetSeries(const StringContainer& files, bool groupImagesWithGantryTilt, const StringContainer &restrictions)
{
  return GetSeries(files, true, groupImagesWithGantryTilt, restrictions);
}
  
DicomSeriesReader::UidFileNamesMap 
DicomSeriesReader::GetSeries(const StringContainer& files, bool sortTo3DPlust, bool groupImagesWithGantryTilt, const StringContainer& /*restrictions*/)
{
  /**
    assumption about this method:
      returns a map of uid-like-key --> list(filename)
      each entry should contain filenames that have images of same
        - series instance uid (automatically done by GDCMSeriesFileNames
        - 0020,0037 image orientation (patient)
        - 0028,0030 pixel spacing (x,y)
        - 0018,0050 slice thickness
  */

  UidFileNamesMap groupsOfSimilarImages; // preliminary result, refined into the final result mapOf3DPlusTBlocks

  // use GDCM directly, itk::GDCMSeriesFileNames does not work with GDCM 2

  // PART I: scan files for sorting relevant DICOM tags, 
  //         separate images that differ in any of those 
  //         attributes (they cannot possibly form a 3D block)

  // scan for relevant tags in dicom files
  gdcm::Scanner scanner;
  const gdcm::Tag tagSeriesInstanceUID(0x0020,0x000e); // Series Instance UID
    scanner.AddTag( tagSeriesInstanceUID );

  const gdcm::Tag tagImageOrientation(0x0020, 0x0037); // image orientation
    scanner.AddTag( tagImageOrientation );

  const gdcm::Tag tagPixelSpacing(0x0028, 0x0030); // pixel spacing
    scanner.AddTag( tagPixelSpacing );
    
  const gdcm::Tag tagSliceThickness(0x0018, 0x0050); // slice thickness
    scanner.AddTag( tagSliceThickness );

  const gdcm::Tag tagNumberOfRows(0x0028, 0x0010); // number rows
    scanner.AddTag( tagNumberOfRows );

  const gdcm::Tag tagNumberOfColumns(0x0028, 0x0011); // number cols
    scanner.AddTag( tagNumberOfColumns );
  
  // additional tags read in this scan to allow later analysis
  // THESE tag are not used for initial separating of files
  const gdcm::Tag tagImagePositionPatient(0x0020,0x0032); // Image Position (Patient)
    scanner.AddTag( tagImagePositionPatient );

  // TODO add further restrictions from arguments

  // let GDCM scan files
  if ( !scanner.Scan( files ) )
  {
    MITK_ERROR << "gdcm::Scanner failed when scanning " << files.size() << " input files.";
    return groupsOfSimilarImages;
  }

  // assign files IDs that will separate them for loading into image blocks
  for (gdcm::Scanner::ConstIterator fileIter = scanner.Begin();
       fileIter != scanner.End();
       ++fileIter)
  {
    //MITK_DEBUG << "Scan file " << fileIter->first << std::endl;
    if ( std::string(fileIter->first).empty() ) continue; // TODO understand why Scanner has empty string entries
    if ( std::string(fileIter->first) == std::string("DICOMDIR") ) continue;

    // we const_cast here, because I could not use a map.at() function in CreateMoreUniqueSeriesIdentifier.
    // doing the same thing with find would make the code less readable. Since we forget the Scanner results
    // anyway after this function, we can simply tolerate empty map entries introduced by bad operator[] access
    std::string moreUniqueSeriesId = CreateMoreUniqueSeriesIdentifier( const_cast<gdcm::Scanner::TagToValue&>(fileIter->second) );
    groupsOfSimilarImages [ moreUniqueSeriesId ].push_back( fileIter->first );
  }
  
  // PART II: sort slices spatially

  for ( UidFileNamesMap::const_iterator groupIter = groupsOfSimilarImages.begin(); groupIter != groupsOfSimilarImages.end(); ++groupIter )
  {
    try
    {
    groupsOfSimilarImages[ groupIter->first ] = SortSeriesSlices( groupIter->second  ); // sort each slice group spatially
    } catch(...)
    {
       MITK_ERROR << "Catched something.";
    }
  }

  // PART III: analyze pre-sorted images for valid blocks (i.e. blocks of equal z-spacing), 
  //          separate into multiple blocks if necessary.
  //          
  //          Analysis performs the following steps:
  //            * imitate itk::ImageSeriesReader: use the distance between the first two images as z-spacing
  //            * check what images actually fulfill ITK's z-spacing assumption
  //            * separate all images that fail the test into new blocks, re-iterate analysis for these blocks

  UidFileNamesMap mapOf3DPlusTBlocks; // final result of this function
  for ( UidFileNamesMap::const_iterator groupIter = groupsOfSimilarImages.begin(); groupIter != groupsOfSimilarImages.end(); ++groupIter )
  {
    UidFileNamesMap mapOf3DBlocks;      // intermediate result for only this group(!)
    std::map<std::string, SliceGroupingAnalysisResult> mapOf3DBlockAnalysisResults;
    StringContainer filesStillToAnalyze = groupIter->second;
    std::string groupUID = groupIter->first;
    unsigned int subgroup(0);
    MITK_DEBUG << "Analyze group " << groupUID;

    while (!filesStillToAnalyze.empty()) // repeat until all files are grouped somehow 
    {
      SliceGroupingAnalysisResult analysisResult = 
        AnalyzeFileForITKImageSeriesReaderSpacingAssumption( filesStillToAnalyze, 
                                                             groupImagesWithGantryTilt, 
                                                             scanner.GetMappings() );

      // enhance the UID for additional groups
      std::stringstream newGroupUID;
      newGroupUID << groupUID << '.' << subgroup;
      mapOf3DBlocks[ newGroupUID.str() ] = analysisResult.GetBlockFilenames();
      MITK_DEBUG << "Result: sorted 3D group " << newGroupUID.str() << " with " << mapOf3DBlocks[ newGroupUID.str() ].size() << " files";
      
      ++subgroup;
        
      filesStillToAnalyze = analysisResult.GetUnsortedFilenames(); // remember what needs further analysis
    }

    // end of grouping, now post-process groups
  
    // PART IV: attempt to group blocks to 3D+t blocks if requested
    //           inspect entries of mapOf3DBlocks
    //            - if number of files is identical to previous entry, collect for 3D+t block
    //            - as soon as number of files changes from previous entry, record collected blocks as 3D+t block, start a new one, continue

    // decide whether or not to group 3D blocks into 3D+t blocks where possible
    if ( !sortTo3DPlust )
    {
      // copy 3D blocks to output
      // TODO avoid collisions (or prove impossibility)
      mapOf3DPlusTBlocks.insert( mapOf3DBlocks.begin(), mapOf3DBlocks.end() );
    }
    else
    {
      // sort 3D+t (as described in "PART IV")
      
      MITK_DEBUG << "================================================================================";
      MITK_DEBUG << "3D+t analysis:";

      unsigned int numberOfFilesInPreviousBlock(0);
      std::string previousBlockKey;

      for ( UidFileNamesMap::const_iterator block3DIter = mapOf3DBlocks.begin();
            block3DIter != mapOf3DBlocks.end();
            ++block3DIter )
      {
        unsigned int numberOfFilesInThisBlock = block3DIter->second.size();
        std::string thisBlockKey = block3DIter->first;
        
        if (numberOfFilesInPreviousBlock == 0)
        {
          numberOfFilesInPreviousBlock = numberOfFilesInThisBlock;
          mapOf3DPlusTBlocks[thisBlockKey].insert( mapOf3DPlusTBlocks[thisBlockKey].end(), 
                                                   block3DIter->second.begin(), 
                                                   block3DIter->second.end() );
          MITK_DEBUG << "  3D+t group " << thisBlockKey << " started";
          previousBlockKey = thisBlockKey;
        }
        else
        {
          bool identicalOrigins;
          try {
            // check whether this and the previous block share a comon origin
            // TODO should be safe, but a little try/catch or other error handling wouldn't hurt
            
            const char 
                *origin_value = scanner.GetValue( mapOf3DBlocks[thisBlockKey].front().c_str(), tagImagePositionPatient ),
                *previous_origin_value = scanner.GetValue( mapOf3DBlocks[previousBlockKey].front().c_str(), tagImagePositionPatient ),
                *destination_value = scanner.GetValue( mapOf3DBlocks[thisBlockKey].back().c_str(), tagImagePositionPatient ),
                *previous_destination_value = scanner.GetValue( mapOf3DBlocks[previousBlockKey].back().c_str(), tagImagePositionPatient );
         
            if (!origin_value || !previous_origin_value || !destination_value || !previous_destination_value)
            {
              identicalOrigins = false;
            } 
            else
            { 
              std::string thisOriginString = ConstCharStarToString( origin_value );
              std::string previousOriginString = ConstCharStarToString( previous_origin_value );

              // also compare last origin, because this might differ if z-spacing is different
              std::string thisDestinationString = ConstCharStarToString( destination_value );
              std::string previousDestinationString = ConstCharStarToString( previous_destination_value );

              identicalOrigins =  ( (thisOriginString == previousOriginString) && (thisDestinationString == previousDestinationString) );
            }

          } catch(...)
          {
            identicalOrigins = false;
          }

          if (identicalOrigins && (numberOfFilesInPreviousBlock == numberOfFilesInThisBlock))
          {
            // group with previous block
            mapOf3DPlusTBlocks[previousBlockKey].insert( mapOf3DPlusTBlocks[previousBlockKey].end(), 
                                                         block3DIter->second.begin(), 
                                                         block3DIter->second.end() );
            MITK_DEBUG << "  --> group enhanced with another timestep";
          }
          else
          {
            // start a new block
            mapOf3DPlusTBlocks[thisBlockKey].insert( mapOf3DPlusTBlocks[thisBlockKey].end(), 
                                                     block3DIter->second.begin(), 
                                                     block3DIter->second.end() );
            MITK_DEBUG << "  ==> group closed with " << mapOf3DPlusTBlocks[previousBlockKey].size() / numberOfFilesInPreviousBlock << " time steps";
            previousBlockKey = thisBlockKey;
            MITK_DEBUG << "  3D+t group " << thisBlockKey << " started";
          }
        }
        
        numberOfFilesInPreviousBlock = numberOfFilesInThisBlock;
      }
    }
  }

  MITK_DEBUG << "================================================================================";
  MITK_DEBUG << "Summary: ";
  for ( UidFileNamesMap::const_iterator groupIter = mapOf3DPlusTBlocks.begin(); groupIter != mapOf3DPlusTBlocks.end(); ++groupIter )
  {
    MITK_DEBUG << "  Image volume " << groupIter->first << " with " << groupIter->second.size() << " files";
  }
  MITK_DEBUG << "Done. ";
  MITK_DEBUG << "================================================================================";

  return mapOf3DPlusTBlocks;
}

DicomSeriesReader::UidFileNamesMap 
DicomSeriesReader::GetSeries(const std::string &dir, bool groupImagesWithGantryTilt, const StringContainer &restrictions)
{
  gdcm::Directory directoryLister;
  directoryLister.Load( dir.c_str(), false ); // non-recursive
  return GetSeries(directoryLister.GetFilenames(), groupImagesWithGantryTilt, restrictions);
}

std::string
DicomSeriesReader::CreateSeriesIdentifierPart( gdcm::Scanner::TagToValue& tagValueMap, const gdcm::Tag& tag )
{
  std::string result;
  try
  {
    result = IDifyTagValue( tagValueMap[ tag ] ? tagValueMap[ tag ] : std::string("") );
  }
  catch (std::exception& e)
  {
    MITK_WARN << "Could not access tag " << tag << ": " << e.what();
  }
   
  return result;
}

std::string 
DicomSeriesReader::CreateMoreUniqueSeriesIdentifier( gdcm::Scanner::TagToValue& tagValueMap )
{
  const gdcm::Tag tagSeriesInstanceUID(0x0020,0x000e); // Series Instance UID
  const gdcm::Tag tagImageOrientation(0x0020, 0x0037); // image orientation
  const gdcm::Tag tagPixelSpacing(0x0028, 0x0030); // pixel spacing
  const gdcm::Tag tagSliceThickness(0x0018, 0x0050); // slice thickness
  const gdcm::Tag tagNumberOfRows(0x0028, 0x0010); // number rows
  const gdcm::Tag tagNumberOfColumns(0x0028, 0x0011); // number cols

  const char* tagSeriesInstanceUid = tagValueMap[tagSeriesInstanceUID];
  if (!tagSeriesInstanceUid)
  {
    mitkThrow() << "CreateMoreUniqueSeriesIdentifier() could not access series instance UID. Something is seriously wrong with this image, so stopping here.";
  }
  std::string constructedID = tagSeriesInstanceUid;

  constructedID += CreateSeriesIdentifierPart( tagValueMap, tagNumberOfRows );
  constructedID += CreateSeriesIdentifierPart( tagValueMap, tagNumberOfColumns );
  constructedID += CreateSeriesIdentifierPart( tagValueMap, tagPixelSpacing );
  constructedID += CreateSeriesIdentifierPart( tagValueMap, tagSliceThickness );
  
  // be a bit tolerant for orienatation, let only the first few digits matter (http://bugs.mitk.org/show_bug.cgi?id=12263)
  // NOT constructedID += CreateSeriesIdentifierPart( tagValueMap, tagImageOrientation );
  try
  {
    bool conversionError(false);
    Vector3D right; right.Fill(0.0);
    Vector3D up; right.Fill(0.0);
    DICOMStringToOrientationVectors( tagValueMap[tagImageOrientation], right, up, conversionError );
    //string  newstring sprintf(simplifiedOrientationString, "%.3f\\%.3f\\%.3f\\%.3f\\%.3f\\%.3f", right[0], right[1], right[2], up[0], up[1], up[2]);
    std::ostringstream ss;
    ss.setf(std::ios::fixed, std::ios::floatfield);
    ss.precision(5);
    ss << right[0] << "\\"
       << right[1] << "\\"
       << right[2] << "\\"
       << up[0] << "\\"
       << up[1] << "\\"
       << up[2];
    std::string simplifiedOrientationString(ss.str());

    constructedID += IDifyTagValue( simplifiedOrientationString );
  }
  catch (std::exception& e)
  {
    MITK_WARN << "Could not access tag " << tagImageOrientation << ": " << e.what();
  }
 

  constructedID.resize( constructedID.length() - 1 ); // cut of trailing '.'

  return constructedID; 
}

std::string 
DicomSeriesReader::IDifyTagValue(const std::string& value)
{
  std::string IDifiedValue( value );
  if (value.empty()) throw std::logic_error("IDifyTagValue() illegaly called with empty tag value");

  // Eliminate non-alnum characters, including whitespace...
  //   that may have been introduced by concats.
  for(std::size_t i=0; i<IDifiedValue.size(); i++)
  {
    while(i<IDifiedValue.size() 
      && !( IDifiedValue[i] == '.'
        || (IDifiedValue[i] >= 'a' && IDifiedValue[i] <= 'z')
        || (IDifiedValue[i] >= '0' && IDifiedValue[i] <= '9')
        || (IDifiedValue[i] >= 'A' && IDifiedValue[i] <= 'Z')))
    {
      IDifiedValue.erase(i, 1);
    }
  }


  IDifiedValue += ".";
  return IDifiedValue;
}

DicomSeriesReader::StringContainer 
DicomSeriesReader::GetSeries(const std::string &dir, const std::string &series_uid, bool groupImagesWithGantryTilt, const StringContainer &restrictions)
{
  UidFileNamesMap allSeries = GetSeries(dir, groupImagesWithGantryTilt, restrictions);
  StringContainer resultingFileList;

  for ( UidFileNamesMap::const_iterator idIter = allSeries.begin(); 
        idIter != allSeries.end(); 
        ++idIter )
  {
    if ( idIter->first.find( series_uid ) == 0 ) // this ID starts with given series_uid
    {
      resultingFileList.insert( resultingFileList.end(), idIter->second.begin(), idIter->second.end() ); // append
    }
  }

  return resultingFileList;
}

DicomSeriesReader::StringContainer 
DicomSeriesReader::SortSeriesSlices(const StringContainer &unsortedFilenames)
{
  gdcm::Sorter sorter;

  sorter.SetSortFunction(DicomSeriesReader::GdcmSortFunction);
  try
  {
    sorter.Sort(unsortedFilenames);
    return sorter.GetFilenames();
  }
  catch(std::logic_error&)
  {
    MITK_WARN << "Sorting error. Leaving series unsorted."; 
    return unsortedFilenames;
  }
}

bool 
DicomSeriesReader::GdcmSortFunction(const gdcm::DataSet &ds1, const gdcm::DataSet &ds2)
{
  // make sure we have Image Position and Orientation
  if ( ! ( 
      ds1.FindDataElement(gdcm::Tag(0x0020,0x0032)) &&
      ds1.FindDataElement(gdcm::Tag(0x0020,0x0037)) &&
      ds2.FindDataElement(gdcm::Tag(0x0020,0x0032)) &&
      ds2.FindDataElement(gdcm::Tag(0x0020,0x0037)) 
        )
      )
      {
        MITK_WARN << "Dicom images are missing attributes for a meaningful sorting.";
        throw std::logic_error("Dicom images are missing attributes for a meaningful sorting.");
      }
  
  gdcm::Attribute<0x0020,0x0032> image_pos1; // Image Position (Patient)
  gdcm::Attribute<0x0020,0x0037> image_orientation1; // Image Orientation (Patient)

  image_pos1.Set(ds1);
  image_orientation1.Set(ds1);

  gdcm::Attribute<0x0020,0x0032> image_pos2;
  gdcm::Attribute<0x0020,0x0037> image_orientation2;

  image_pos2.Set(ds2);
  image_orientation2.Set(ds2);

  /*
     we tolerate very small differences in image orientation, since we got to know about
     acquisitions where these values change across a single series (7th decimal digit)
     (http://bugs.mitk.org/show_bug.cgi?id=12263)

     still, we want to check if our assumption of 'almost equal' orientations is valid
  */
  for (unsigned int dim = 0; dim < 6; ++dim)
  {
    if ( fabs(image_orientation2[dim] - image_orientation1[dim]) > 0.0001 )
    {
      MITK_ERROR << "Dicom images have different orientations.";
      throw std::logic_error("Dicom images have different orientations. Call GetSeries() first to separate images.");
    }
  }

  double normal[3];

  normal[0] = image_orientation1[1] * image_orientation1[5] - image_orientation1[2] * image_orientation1[4];
  normal[1] = image_orientation1[2] * image_orientation1[3] - image_orientation1[0] * image_orientation1[5];
  normal[2] = image_orientation1[0] * image_orientation1[4] - image_orientation1[1] * image_orientation1[3];

  double
      dist1 = 0.0,
      dist2 = 0.0;

  for (unsigned char i = 0u; i < 3u; ++i)
  {
    dist1 += normal[i] * image_pos1[i];
    dist2 += normal[i] * image_pos2[i];
  }

  if ( fabs(dist1 - dist2) < mitk::eps)
  {
    gdcm::Attribute<0x0008,0x0032> acq_time1;   // Acquisition time (may be missing, so we check existence first)
    gdcm::Attribute<0x0008,0x0032> acq_time2;

    gdcm::Attribute<0x0020,0x0012> acq_number1; // Acquisition number (may also be missing, so we check existence first)
    gdcm::Attribute<0x0020,0x0012> acq_number2;

    if (ds1.FindDataElement(gdcm::Tag(0x0008,0x0032)) && ds2.FindDataElement(gdcm::Tag(0x0008,0x0032)))
    {
      acq_time1.Set(ds1);
      acq_time2.Set(ds2);

      return acq_time1 < acq_time2;
    }
    else if (ds1.FindDataElement(gdcm::Tag(0x0020,0x0012)) && ds2.FindDataElement(gdcm::Tag(0x0020,0x0012)))
    {
      acq_number1.Set(ds1);
      acq_number2.Set(ds2);

      return acq_number1 < acq_number2;
    }
    else
    {
      return true;
    }
  }
  else
  {
    // default: compare position
    return dist1 < dist2;
  }
}
  
std::string DicomSeriesReader::GetConfigurationString()
{
  std::stringstream configuration;
  configuration << "MITK_USE_GDCMIO: ";
  configuration << "true";
  configuration << "\n";

  configuration << "GDCM_VERSION: ";
#ifdef GDCM_MAJOR_VERSION
  configuration << GDCM_VERSION;
#endif
  //configuration << "\n";

  return configuration.str();
}

void DicomSeriesReader::CopyMetaDataToImageProperties(StringContainer filenames, const gdcm::Scanner::MappingType &tagValueMappings_, DcmIoType *io, Image *image)
{
  std::list<StringContainer> imageBlock;
  imageBlock.push_back(filenames);
  CopyMetaDataToImageProperties(imageBlock, tagValueMappings_, io, image);
}
  
void DicomSeriesReader::CopyMetaDataToImageProperties( std::list<StringContainer> imageBlock, const gdcm::Scanner::MappingType& tagValueMappings_,  DcmIoType* io, Image* image)
{
  if (!io || !image) return;

  StringLookupTable filesForSlices;
  StringLookupTable sliceLocationForSlices;
  StringLookupTable instanceNumberForSlices;
  StringLookupTable SOPInstanceNumberForSlices;

  gdcm::Scanner::MappingType& tagValueMappings = const_cast<gdcm::Scanner::MappingType&>(tagValueMappings_);

  //DICOM tags which should be added to the image properties
  const gdcm::Tag tagSliceLocation(0x0020, 0x1041); // slice location

  const gdcm::Tag tagInstanceNumber(0x0020, 0x0013); // (image) instance number

  const gdcm::Tag tagSOPInstanceNumber(0x0008, 0x0018); // SOP instance number
  unsigned int timeStep(0);

  std::string propertyKeySliceLocation = "dicom.image.0020.1041";
  std::string propertyKeyInstanceNumber = "dicom.image.0020.0013";
  std::string propertyKeySOPInstanceNumber = "dicom.image.0008.0018";

  // tags for each image
  for ( std::list<StringContainer>::iterator i = imageBlock.begin(); i != imageBlock.end(); i++, timeStep++ )
  {

    const StringContainer& files = (*i);
    unsigned int slice(0);
    for ( StringContainer::const_iterator fIter = files.begin();
          fIter != files.end();
          ++fIter, ++slice )
    {
      filesForSlices.SetTableValue( slice, *fIter );
      gdcm::Scanner::TagToValue tagValueMapForFile = tagValueMappings[fIter->c_str()];
      if(tagValueMapForFile.find(tagSliceLocation) != tagValueMapForFile.end())
        sliceLocationForSlices.SetTableValue(slice, tagValueMapForFile[tagSliceLocation]);
      if(tagValueMapForFile.find(tagInstanceNumber) != tagValueMapForFile.end())
        instanceNumberForSlices.SetTableValue(slice, tagValueMapForFile[tagInstanceNumber]);
      if(tagValueMapForFile.find(tagSOPInstanceNumber) != tagValueMapForFile.end())
        SOPInstanceNumberForSlices.SetTableValue(slice, tagValueMapForFile[tagSOPInstanceNumber]);
    }

    image->SetProperty( "files", StringLookupTableProperty::New( filesForSlices ) );

    //If more than one time step add postfix ".t" + timestep
    if(timeStep != 0)
    {
      propertyKeySliceLocation.append(".t" + timeStep);
      propertyKeyInstanceNumber.append(".t" + timeStep);
      propertyKeySOPInstanceNumber.append(".t" + timeStep);
    }
    image->SetProperty( propertyKeySliceLocation.c_str(), StringLookupTableProperty::New( sliceLocationForSlices ) );
    image->SetProperty( propertyKeyInstanceNumber.c_str(), StringLookupTableProperty::New( instanceNumberForSlices ) );
    image->SetProperty( propertyKeySOPInstanceNumber.c_str(), StringLookupTableProperty::New( SOPInstanceNumberForSlices ) );
  }

  // Copy tags for series, study, patient level (leave interpretation to application).
  // These properties will be copied to the DataNode by DicomSeriesReader.

  // tags for the series (we just use the one that ITK copied to its dictionary (proably that of the last slice)
  const itk::MetaDataDictionary& dict = io->GetMetaDataDictionary();
  const TagToPropertyMapType& propertyLookup = DicomSeriesReader::GetDICOMTagsToMITKPropertyMap();

  itk::MetaDataDictionary::ConstIterator dictIter = dict.Begin();
  while ( dictIter != dict.End() )
  {
    //MITK_DEBUG << "Key " << dictIter->first;
    std::string value;
    if ( itk::ExposeMetaData<std::string>( dict, dictIter->first, value ) )
    {
      //MITK_DEBUG << "Value " << value;

      TagToPropertyMapType::const_iterator valuePosition = propertyLookup.find( dictIter->first );
      if ( valuePosition != propertyLookup.end() )
      {
        std::string propertyKey = valuePosition->second;
        //MITK_DEBUG << "--> " << propertyKey;
        
        image->SetProperty( propertyKey.c_str(), StringProperty::New(value) );
      }
    }
    else
    {
      MITK_WARN << "Tag " << dictIter->first << " not read as string as expected. Ignoring..." ;
    }
    ++dictIter;
  }
}
  

} // end namespace mitk

#include <mitkDicomSeriesReader.txx>

