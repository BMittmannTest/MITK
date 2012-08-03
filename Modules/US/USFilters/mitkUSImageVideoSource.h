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


#ifndef MITKUSImageVideoSource_H_HEADER_INCLUDED_
#define MITKUSImageVideoSource_H_HEADER_INCLUDED_

#include <itkProcessObject.h>
#include "mitkUSImage.h"
#include "mitkOpenCVVideoSource.h"
#include "mitkOpenCVToMitkImageFilter.h"

namespace mitk {

  /**Documentation
  * \brief This class can be pointed to a video file or a videodevice and delivers USImages with default metadata Sets
  *
  * \ingroup US
  */
  class MitkUS_EXPORT USImageVideoSource : public itk::Object
  {
  public:
    mitkClassMacro(USImageVideoSource, itk::ProcessObject);
    itkNewMacro(Self);

    /**
    *\brief Opens a video file for streaming. If nothing goes wrong, the 
    * VideoSource is ready to deliver images after calling this function.
    */
    void SetVideoFileInput(std::string path);

    /**
    *\brief Opens a video device for streaming. Takes the Device id. Try -1 for "grab the first you can get"
    * which works quite well if only one device is available. If nothing goes wrong, the 
    * VideoSource is ready to deliver images after calling this function.
    */
    void SetCameraInput(int deviceID);

    /**
    *\brief Sets the output image to rgb or grayscale. Output is grayscale by default
    * and can be set to color by passing true, or to grayscale again by passing false.
    */
    void SetColorOutput(bool isColor);

    /**
    * /brief Defines the cropping area. The rectangle will be justified to the image borders
    * if the given rectangle is larger than the video source. If a correct rectangle is given,
    * The dimensions of the output image will be equal to those of the rectangle.
    */
    void SetRegionOfInterest(int topLeftX, int topLeftY, int bottomRightX, int bottomRightY);

    /**
    * /brief Removes the region of interest. Produced images will be uncropped after call.
    */
    void RemoveRegionOfInterest();

    /**
    *\brief Retrieves the next frame. This will typically be the next frame in a file
    * or the last cahced file in a devcie.
    */
    mitk::USImage::Pointer GetNextImage();
    
    // Getter & Setter
    itkGetMacro(OpenCVVideoSource, mitk::OpenCVVideoSource::Pointer);
    itkSetMacro(OpenCVVideoSource, mitk::OpenCVVideoSource::Pointer);
    itkGetMacro(IsVideoReady, bool);


  protected:
    USImageVideoSource();
    virtual ~USImageVideoSource();

    /**
    * \brief The source of the video
    */
    mitk::OpenCVVideoSource::Pointer m_OpenCVVideoSource;

   bool m_IsVideoReady;
   bool m_IsGreyscale;
   cv::Rect m_CropRegion;
   mitk::OpenCVToMitkImageFilter::Pointer m_OpenCVToMitkFilter;
    

  };
} // namespace mitk
#endif /* MITKUSImageVideoSource_H_HEADER_INCLUDED_ */
