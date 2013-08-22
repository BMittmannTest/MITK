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

#ifndef mitkContourUtilshIncludett
#define mitkContourUtilshIncludett

#include "mitkImage.h"
#include "SegmentationExports.h"
#include "mitkContourModel.h"
#include "mitkLabelSetImage.h"
#include "mitkLegacyAdaptors.h"

#include <itkImage.h>

namespace mitk
{

/**
 * \brief Helpful methods for working with contours and images
 */
class Segmentation_EXPORT ContourUtils : public itk::Object
{
  public:

    mitkClassMacro(ContourUtils, itk::Object);
    itkNewMacro(ContourUtils);

    /**
      \brief Projects a contour onto an image point by point. Converts from world to index coordinates.

      \param correctionForIpSegmentation adds 0.5 to x and y index coordinates (difference between ipSegmentation and MITK contours)
    */
    static ContourModel::Pointer ProjectContourTo2DSlice(Image* slice, ContourModel* contourIn3D, bool constrainToInside);

    /**
      \brief Projects a slice index coordinates of a contour back into world coordinates.

      \param correctionForIpSegmentation subtracts 0.5 to x and y index coordinates (difference between ipSegmentation and MITK contours)
    */
    static ContourModel::Pointer BackProjectContourFrom2DSlice(const Geometry3D* sliceGeometry, ContourModel* contourIn2D);

    /**
      \brief Fill a contour in a 2D slice with a specified pixel value.
    */
    static void FillContourInSlice( ContourModel* projectedContour, Image* slice, const LabelSet* labelSet, int paintingPixelValue = 1 );

  protected:

    ContourUtils();
    virtual ~ContourUtils();

    /**
      \brief Paint a filled contour (e.g. of an ipSegmentation pixel type) into a mitk::Image (or arbitraty pixel type).
      Will not copy the whole filledContourSlice, but only set those pixels in originalSlice to overwritevalue, where the corresponding pixel
      in filledContourSlice is non-zero.
    */
    template<typename TPixel, unsigned int VImageDimension>
    void ItkCopyFilledContourToSlice( itk::Image<TPixel,VImageDimension>* originalSlice, const Image* filledContourSlice, int overwritevalue = 1 );

    /**
      \brief Paint a filled contour (e.g. of an ipSegmentation pixel type) into a mitk::Image (or arbitraty pixel type).
      Will not copy the whole filledContourSlice, but only set those pixels in originalSlice that can be overwritten acording to
      the label "locked" property to overwritevalue, where the corresponding pixel in filledContourSlice is non-zero.
    */
    template<typename TPixel, unsigned int VImageDimension>
    static void ItkCopyFilledContourToSlice2( itk::Image<TPixel,VImageDimension>* originalSlice, const LabelSetImage* filledContourSlice, int overwritevalue = 1 );
};

}

#endif

