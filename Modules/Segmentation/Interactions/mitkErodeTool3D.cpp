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

#include "mitkErodeTool3D.h"

#include "mitkBaseRenderer.h"
#include "mitkRenderingManager.h"
#include "mitkInteractionConst.h"
#include "mitkImageAccessByItk.h"
#include "mitkToolManager.h"
#include "mitkImageCast.h"

// itk
#include <itkBinaryThresholdImageFilter.h>
#include <itkBinaryBallStructuringElement.h>
#include <itkBinaryErodeImageFilter.h>
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
  MITK_TOOL_MACRO(Segmentation_EXPORT, ErodeTool3D, "ErodeTool3D tool");
}


mitk::ErodeTool3D::ErodeTool3D()
{
  m_ProgressCommand = mitk::ToolCommand::New();
  m_Radius = 1;
}

mitk::ErodeTool3D::~ErodeTool3D()
{
}

const char** mitk::ErodeTool3D::GetXPM() const
{
  return NULL;//mitkErodeTool3D_xpm;
}

us::ModuleResource mitk::ErodeTool3D::GetIconResource() const
{
  us::Module* module = us::GetModuleContext()->GetModule();
  assert(module);
  us::ModuleResource resource = module->GetResource("ErodeTool3D_48x48.png");
  return resource;
}

const char* mitk::ErodeTool3D::GetName() const
{
  return "ErodeTool3D";
}

void mitk::ErodeTool3D::SetRadius(int value)
{
  m_Radius = value;
}

void mitk::ErodeTool3D::Run()
{
//  this->InitializeUndoController();

  mitk::DataNode* workingNode = m_ToolManager->GetWorkingData(0);
  assert(workingNode);

  mitk::LabelSetImage* workingImage = dynamic_cast< mitk::LabelSetImage* >( workingNode->GetData() );
  assert(workingImage);

  m_OverwritePixelValue = workingImage->GetActiveLabelIndex();

  // todo: use it later
  unsigned int timestep = mitk::RenderingManager::GetInstance()->GetTimeNavigationController()->GetTime()->GetPos();

  CurrentlyBusy.Send(true);

  try
  {
    AccessByItk(workingImage, InternalProcessing);
  }
  catch( itk::ExceptionObject& e )
  {
   MITK_ERROR << "Exception caught: " << e.GetDescription();
   m_ProgressCommand->Reset();
   CurrentlyBusy.Send(false);
   return;
  }

  CurrentlyBusy.Send(false);

  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}


template < typename TPixel, unsigned int VDimension >
void mitk::ErodeTool3D::InternalProcessing( itk::Image< TPixel, VDimension>* input )
{
  typedef itk::Image<TPixel, VDimension> ImageType;
  typedef itk::LabelObject< TPixel, VDimension > LabelObjectType;
  typedef itk::LabelMap< LabelObjectType > LabelMapType;
  typedef itk::LabelImageToLabelMapFilter< ImageType, LabelMapType > Image2LabelMapType;
  typedef itk::AutoCropLabelMapFilter< LabelMapType > AutoCropType;
  typedef itk::LabelMapToLabelImageFilter< LabelMapType, ImageType > LabelMap2ImageType;
  typedef itk::BinaryBallStructuringElement<TPixel, VDimension> BallType;
  typedef itk::BinaryErodeImageFilter<ImageType, ImageType, BallType> ErodeFilterType;
  typedef itk::BinaryThresholdImageFilter< ImageType, ImageType > ThresholdFilterType;

  mitk::DataNode* workingNode = m_ToolManager->GetWorkingData(0);
  assert(workingNode);

  mitk::LabelSetImage* workingImage = dynamic_cast< mitk::LabelSetImage* >( workingNode->GetData() );
  assert(workingImage);

  typename ThresholdFilterType::Pointer thresholdFilter = ThresholdFilterType::New();
  thresholdFilter->SetInput(input);
  thresholdFilter->SetLowerThreshold(m_OverwritePixelValue);
  thresholdFilter->SetUpperThreshold(m_OverwritePixelValue);
  thresholdFilter->SetOutsideValue(0);
  thresholdFilter->SetInsideValue(1);

  typename Image2LabelMapType::Pointer image2label = Image2LabelMapType::New();
  image2label->SetInput(thresholdFilter->GetOutput());

  typename AutoCropType::SizeType border;
  border[0] = 3;
  border[1] = 3;
  border[2] = 3;

  typename AutoCropType::Pointer autoCropFilter = AutoCropType::New();
  autoCropFilter->SetInput( image2label->GetOutput() );
  autoCropFilter->SetCropBorder(border);
  autoCropFilter->InPlaceOn();

  typename LabelMap2ImageType::Pointer label2image = LabelMap2ImageType::New();
  label2image->SetInput( autoCropFilter->GetOutput() );

  BallType ball;
  ball.SetRadius(m_Radius);
  ball.CreateStructuringElement();

  typename ErodeFilterType::Pointer erodeFilter = ErodeFilterType::New();
  erodeFilter->SetKernel(ball);
  erodeFilter->SetInput(label2image->GetOutput());
  erodeFilter->SetErodeValue(1);

  if (m_ProgressCommand.IsNotNull())
  {
    thresholdFilter->AddObserver( itk::AnyEvent(), m_ProgressCommand );
    autoCropFilter->AddObserver( itk::AnyEvent(), m_ProgressCommand );
    erodeFilter->AddObserver( itk::AnyEvent(), m_ProgressCommand );
    m_ProgressCommand->AddStepsToDo(200);
  }

  erodeFilter->Update();

  if (m_ProgressCommand.IsNotNull())
  {
    m_ProgressCommand->Reset();
  }

  typename ImageType::Pointer result = erodeFilter->GetOutput();
  result->DisconnectPipeline();

  m_PreviewImage = mitk::Image::New();
  m_PreviewImage->InitializeByItk(result.GetPointer());
  m_PreviewImage->SetChannel(result->GetBufferPointer());

  typename ImageType::RegionType cropRegion;
  cropRegion = autoCropFilter->GetOutput()->GetLargestPossibleRegion();

  const typename ImageType::SizeType& cropSize = cropRegion.GetSize();
  const typename ImageType::IndexType& cropIndex = cropRegion.GetIndex();

  mitk::SlicedGeometry3D* slicedGeometry = m_PreviewImage->GetSlicedGeometry();
  mitk::Point3D origin;
  vtk2itk(cropIndex, origin);
  workingImage->GetSlicedGeometry()->IndexToWorld(origin, origin);
  slicedGeometry->SetOrigin(origin);

  m_PreviewNode->SetData(m_PreviewImage);

  m_RequestedRegion = workingImage->GetLargestPossibleRegion();

  m_RequestedRegion.SetIndex(0,cropIndex[0]);
  m_RequestedRegion.SetIndex(1,cropIndex[1]);
  m_RequestedRegion.SetIndex(2,cropIndex[2]);

  m_RequestedRegion.SetSize(0,cropSize[0]);
  m_RequestedRegion.SetSize(1,cropSize[1]);
  m_RequestedRegion.SetSize(2,cropSize[2]);
}
