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

#include "mitkMedianTool3D.h"

#include "mitkBaseRenderer.h"
#include "mitkRenderingManager.h"
#include "mitkInteractionConst.h"
#include "mitkApplyDiffImageOperation.h"
#include "mitkOperationEvent.h"
#include "mitkDiffImageApplier.h"
#include "mitkUndoController.h"
#include "mitkImageAccessByItk.h"
#include "mitkToolManager.h"
#include "mitkImageCast.h"
#include "mitkImageTimeSelector.h"

#include <itkBinaryThresholdImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryMedianImageFilter.h>
#include <itkLabelObject.h>
#include <itkLabelMap.h>
#include <itkLabelImageToLabelMapFilter.h>
#include <itkAutoCropLabelMapFilter.h>
#include <itkLabelMapToLabelImageFilter.h>

// us
#include <usModule.h>
#include <usModuleResource.h>
#include <usGetModuleContext.h>
#include <usModuleContext.h>

namespace mitk {
  MITK_TOOL_MACRO(Segmentation_EXPORT, MedianTool3D, "MedianTool3D tool");
}

mitk::MedianTool3D::MedianTool3D()
{
}

mitk::MedianTool3D::~MedianTool3D()
{
}

const char** mitk::MedianTool3D::GetXPM() const
{
  return NULL;//mitkMedianTool3D_xpm;
}

us::ModuleResource mitk::MedianTool3D::GetIconResource() const
{
  us::Module* module = us::GetModuleContext()->GetModule();
  us::ModuleResource resource = module->GetResource("MedianTool3D_48x48.png");
  return resource;
}

const char* mitk::MedianTool3D::GetName() const
{
  return "MedianTool3D";
}

void mitk::MedianTool3D::Run()
{
  mitk::DataNode* workingNode = m_ToolManager->GetWorkingData(0);
  assert(workingNode);

  mitk::LabelSetImage* workingImage = dynamic_cast< mitk::LabelSetImage* >( workingNode->GetData() );
  assert(workingImage);

  // todo: use it later
  unsigned int timestep = mitk::RenderingManager::GetInstance()->GetTimeNavigationController()->GetTime()->GetPos();

  m_ProgressCommand = mitk::ToolCommand::New();

  CurrentlyBusy.Send(true);

  try
  {
    AccessByItk(workingImage, ITKProcessing);
  }
  catch( itk::ExceptionObject& e )
  {
   MITK_ERROR << "Exception caught: " << e.GetDescription();
   m_ProgressCommand->Reset();
   CurrentlyBusy.Send(false);
   return;
  }

  CurrentlyBusy.Send(false);

  workingImage->Modified();
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}


template < typename TPixel, unsigned int VDimension >
void mitk::MedianTool3D::ITKProcessing( itk::Image< TPixel, VDimension>* input )
{
  typedef itk::Image<TPixel, VDimension> ImageType;

  typedef itk::LabelObject< TPixel, VDimension > LabelObjectType;
  typedef itk::LabelMap< LabelObjectType > LabelMapType;
  typedef itk::LabelImageToLabelMapFilter< ImageType, LabelMapType > Image2LabelMapType;
  typedef itk::AutoCropLabelMapFilter< LabelMapType > AutoCropType;
  typedef itk::LabelMapToLabelImageFilter< LabelMapType, ImageType > LabelMap2ImageType;
  typedef itk::BinaryBallStructuringElement<TPixel, VDimension> BallType;
  typedef itk::BinaryThresholdImageFilter< ImageType, ImageType > ThresholdFilterType;
  typedef itk::BinaryMedianImageFilter< ImageType, ImageType > MedianFilterType;

  mitk::DataNode* workingNode = m_ToolManager->GetWorkingData(0);
  assert(workingNode);

  mitk::LabelSetImage* workingImage = dynamic_cast< mitk::LabelSetImage* >( workingNode->GetData() );
  assert(workingImage);

  int pixelValue = workingImage->GetActiveLabelIndex();

  typename ThresholdFilterType::Pointer thresholdFilter = ThresholdFilterType::New();
  thresholdFilter->SetInput(input);
  thresholdFilter->SetLowerThreshold(pixelValue);
  thresholdFilter->SetUpperThreshold(pixelValue);
  thresholdFilter->SetOutsideValue(0);
  thresholdFilter->SetInsideValue(pixelValue);

  typename Image2LabelMapType::Pointer image2label = Image2LabelMapType::New();
  image2label->SetInput(thresholdFilter->GetOutput());

  typename AutoCropType::Pointer autoCropFilter = AutoCropType::New();
  autoCropFilter->SetInput( image2label->GetOutput() );
  autoCropFilter->InPlaceOn();

  typename LabelMap2ImageType::Pointer label2image = LabelMap2ImageType::New();
  label2image->SetInput( autoCropFilter->GetOutput() );

  typename ImageType::SizeType radius;
  radius.Fill(1);

  typename MedianFilterType::Pointer medianFilter = MedianFilterType::New();
  medianFilter->SetInput( label2image->GetOutput() );
  medianFilter->SetForegroundValue( pixelValue );
  medianFilter->SetRadius( radius );

  if (m_ProgressCommand.IsNotNull())
  {
    thresholdFilter->AddObserver( itk::AnyEvent(), m_ProgressCommand );
    autoCropFilter->AddObserver( itk::AnyEvent(), m_ProgressCommand );
    medianFilter->AddObserver( itk::AnyEvent(), m_ProgressCommand );
    m_ProgressCommand->AddStepsToDo(100);
  }

  medianFilter->Update();

  if (m_ProgressCommand.IsNotNull())
  {
    m_ProgressCommand->Reset();
  }

  typename ImageType::Pointer result = medianFilter->GetOutput();
  result->DisconnectPipeline();

  typedef itk::ImageRegionConstIterator< ImageType > SourceIteratorType;
  typedef itk::ImageRegionIterator< ImageType > TargetIteratorType;

  typename ImageType::RegionType cropRegion;
  cropRegion = autoCropFilter->GetOutput()->GetLargestPossibleRegion();

  typename const ImageType::SizeType& cropSize = cropRegion.GetSize();
  typename const ImageType::IndexType& cropIndex = cropRegion.GetIndex();

  SourceIteratorType sourceIter( result, result->GetLargestPossibleRegion() );
  sourceIter.GoToBegin();

  TargetIteratorType targetIter( input, cropRegion );
  targetIter.GoToBegin();

  while ( !sourceIter.IsAtEnd() )
  {
    int targetValue = static_cast< int >( targetIter.Get() );
    int sourceValue = static_cast< int >( sourceIter.Get() );

    if ( (targetValue == pixelValue) || ( sourceValue && (!workingImage->GetLabelLocked(targetValue))) )
    {
      targetIter.Set( sourceValue );
    }

    ++sourceIter;
    ++targetIter;
  }
}
