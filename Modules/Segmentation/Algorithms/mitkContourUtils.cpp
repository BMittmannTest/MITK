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

#include "mitkContourUtils.h"
#include "mitkImageCast.h"
#include "mitkImageAccessByItk.h"
#include "mitkInstantiateAccessFunctions.h"
#include "ipSegmentation.h"

#define InstantiateAccessFunction_ItkCopyFilledContourToSlice2(pixelType, dim) \
  template void mitk::ContourUtils::ItkCopyFilledContourToSlice2(itk::Image<pixelType,dim>*, const mitk::LabelSetImage*, int pixelvalue);

// explicitly instantiate the 2D version of this method
InstantiateAccessFunctionForFixedDimension(ItkCopyFilledContourToSlice2, 2);

mitk::ContourUtils::ContourUtils()
{
}

mitk::ContourUtils::~ContourUtils()
{
}

mitk::ContourModel::Pointer mitk::ContourUtils::ProjectContourTo2DSlice(Image* slice, ContourModel* contourIn3D, int timestep)
{
  if ( !slice || !contourIn3D ) return NULL;

  ContourModel::Pointer projectedContour = ContourModel::New();
  projectedContour->Initialize();
  projectedContour->Expand(timestep+1);

  // fixme: we may consider the timestamp in the next line
  const Geometry3D* sliceGeometry = slice->GetGeometry(0);
  if (!sliceGeometry) return NULL;

  for ( ContourModel::ConstVertexIterator iter = contourIn3D->Begin(timestep);
        iter != contourIn3D->End(timestep);
        ++iter )
  {
    Point3D currentPointIn3D;
    for (int i = 0; i < 3; ++i)
      currentPointIn3D[i] = (*iter)->Coordinates[i];

    Point3D projectedPointIn2D;
    projectedPointIn2D.Fill(0.0);
    sliceGeometry->WorldToIndex( currentPointIn3D, projectedPointIn2D );
    projectedContour->AddVertex( projectedPointIn2D, timestep );
  }

  return projectedContour;
}

void mitk::ContourUtils::BackProjectContourFrom2DSlice(const Geometry3D* sliceGeometry, ContourModel* contourIn2D, ContourModel* contourIn3D, int timestep )
{
  if ( !sliceGeometry || !contourIn2D ) return;

  contourIn3D->Clear(timestep);
  contourIn3D->Initialize();
  contourIn3D->Expand(timestep+1);
/*
  ContourModel::Pointer worldContour = ContourModel::New();
  worldContour->Initialize();
  worldContour->Expand(timestep+1);
  worldContour->Close(timestep);
*/
  ContourModel::ConstVertexIterator iter = contourIn2D->Begin(timestep);
  for ( ; iter != contourIn2D->End(timestep); ++iter)
  {
    Point3D worldPointIn3D;
    sliceGeometry->IndexToWorld( (*iter)->Coordinates, worldPointIn3D );
    contourIn3D->AddVertex(worldPointIn3D,timestep);
    //worldContour->AddVertex( worldPointIn3D, timestep );
  }

 // return worldContour;
}

void mitk::ContourUtils::FillContourInSlice( ContourModel* projectedContour, LabelSetImage* slice, int paintingPixelValue, int timestep )
{
  // 1. Use ipSegmentation to draw a filled(!) contour into a new 8 bit 2D image, which will later be copied back to the slice.
  //    We don't work on the "real" working data, because ipSegmentation would restrict us to 8 bit images

  // convert the projected contour into a ipSegmentation format
  mitkIpInt4_t* picContour = new mitkIpInt4_t[2 * projectedContour->GetNumberOfVertices(timestep)];
//  const Contour::PathType::VertexListType* pointsIn2D = projectedContour->GetContourPath()->GetVertexList();
  unsigned int index(0);
  for ( ContourModel::ConstVertexIterator iter = projectedContour->Begin(timestep);
        iter != projectedContour->End(timestep);
        ++iter, ++index )
  {
    picContour[ 2 * index + 0 ] = static_cast<mitkIpInt4_t>( (*iter)->Coordinates[0] + 1.0 ); // +0.5 wahrscheinlich richtiger
    picContour[ 2 * index + 1 ] = static_cast<mitkIpInt4_t>( (*iter)->Coordinates[1] + 1.0 );
    //MITK_INFO << "mitk 2d [" << (*iter)[0] << ", " << (*iter)[1] << "]  pic [" << picContour[ 2*index+0] << ", " << picContour[ 2*index+1] << "]";
  }

  assert( slice->GetSliceData() );
  mitkIpPicDescriptor* originalPicSlice = mitkIpPicNew();
  CastToIpPicDescriptor( slice, originalPicSlice);
  mitkIpPicDescriptor* picSlice = ipMITKSegmentationNew( originalPicSlice );
  ipMITKSegmentationClear( picSlice );

  assert( originalPicSlice && picSlice );

  // here comes the actual contour filling algorithm (from ipSegmentation/Graphics Gems)
  ipMITKSegmentationCombineRegion ( picSlice, picContour, projectedContour->GetNumberOfVertices(timestep), NULL, IPSEGMENTATION_OR,  1); // set to 1

  delete[] picContour;

  // 2. Copy the filled contour to the working data slice
  //    copy all pixels that are non-zero to the original image (not caring for the actual type of the working image). perhaps make the replace value a parameter ( -> general painting tool ).
  //    make the pic slice an mitk/itk image (as little ipPic code as possible), call a templated method with accessbyitk, iterate over the original and the modified slice

  LabelSetImage::Pointer ipsegmentationModifiedSlice = LabelSetImage::New();
  ipsegmentationModifiedSlice->Initialize( CastToImageDescriptor( picSlice ) );
  ipsegmentationModifiedSlice->SetSlice( picSlice->data );
  ipsegmentationModifiedSlice->SetLabelSet( slice->GetLabelSet() );

  //AccessFixedDimensionByItk_2( slice, ItkCopyFilledContourToSlice, 2, ipsegmentationModifiedSlice, paintingPixelValue );
  AccessFixedDimensionByItk_2( slice, ItkCopyFilledContourToSlice2, 2, ipsegmentationModifiedSlice, paintingPixelValue );

  ipsegmentationModifiedSlice = NULL; // free MITK header information
  ipMITKSegmentationFree( picSlice ); // free actual memory
}

template<typename TPixel, unsigned int VImageDimension>
void mitk::ContourUtils::ItkCopyFilledContourToSlice( itk::Image<TPixel,VImageDimension>* originalSlice, const Image* filledContourSlice, int overwritevalue )
{
  typedef itk::Image<TPixel,VImageDimension> SliceType;

  typename SliceType::Pointer filledContourSliceITK;
  CastToItkImage( filledContourSlice, filledContourSliceITK );

  // now the original slice and the ipSegmentation-painted slice are in the same format, and we can just copy all pixels that are non-zero
  typedef itk::ImageRegionIterator< SliceType >       OutputIteratorType;
  typedef itk::ImageRegionConstIterator< SliceType >   InputIteratorType;

  InputIteratorType inputIterator( filledContourSliceITK, filledContourSliceITK->GetLargestPossibleRegion() );
  OutputIteratorType outputIterator( originalSlice, originalSlice->GetLargestPossibleRegion() );

  outputIterator.GoToBegin();
  inputIterator.GoToBegin();

  while ( !outputIterator.IsAtEnd() )
  {
    if ( inputIterator.Get() != 0 )
    {
      outputIterator.Set( overwritevalue );
    }

    ++outputIterator;
    ++inputIterator;
  }
}

template<typename TPixel, unsigned int VImageDimension>
void mitk::ContourUtils::ItkCopyFilledContourToSlice2( itk::Image<TPixel,VImageDimension>* originalSlice, const LabelSetImage* filledContourSlice, int overwritevalue )
{
  typedef itk::Image<TPixel,VImageDimension> SliceType;

  typename SliceType::Pointer filledContourSliceITK;
  CastToItkImage( filledContourSlice, filledContourSliceITK );

  // now the original slice and the ipSegmentation-painted slice are in the same format, and we can just copy all pixels that are non-zero
  typedef itk::ImageRegionIterator< SliceType >       OutputIteratorType;
  typedef itk::ImageRegionConstIterator< SliceType >   InputIteratorType;

  InputIteratorType inputIterator( filledContourSliceITK, filledContourSliceITK->GetLargestPossibleRegion() );
  OutputIteratorType outputIterator( originalSlice, originalSlice->GetLargestPossibleRegion() );

  outputIterator.GoToBegin();
  inputIterator.GoToBegin();

  const int& activePixelValue = filledContourSlice->GetActiveLabelIndex();

  if (overwritevalue != 0)
  {
      while ( !outputIterator.IsAtEnd() )
      {
        const int targetValue = outputIterator.Get();
        if ( inputIterator.Get() != 0 )
        {
          if (!filledContourSlice->GetLabelLocked(targetValue))
            outputIterator.Set( overwritevalue );
        }

        ++outputIterator;
        ++inputIterator;
      }
  }
  else // we are erasing
  {
      while ( !outputIterator.IsAtEnd() )
      {
        const int targetValue = outputIterator.Get();
        if (inputIterator.Get() != 0)
        {
          if (targetValue == activePixelValue)
            outputIterator.Set( overwritevalue );
        }

        ++outputIterator;
        ++inputIterator;
      }
  }
}
