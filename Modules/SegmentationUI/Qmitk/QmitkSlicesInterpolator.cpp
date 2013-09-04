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

#include "QmitkSlicesInterpolator.h"

#include "QmitkStdMultiWidget.h"
#include "QmitkSelectableGLWidget.h"

#include "mitkColorProperty.h"
#include "mitkProperties.h"
#include "mitkRenderingManager.h"
#include "mitkProgressBar.h"
#include "mitkGlobalInteraction.h"
#include "mitkOperationEvent.h"
#include "mitkInteractionConst.h"
#include "mitkApplyDiffImageOperation.h"
#include "mitkDiffImageApplier.h"
#include <mitkDiffSliceOperationApplier.h>
#include "mitkUndoController.h"
#include "mitkSegTool2D.h"
#include "mitkSurfaceToImageFilter.h"
#include "mitkSliceNavigationController.h"
#include <mitkVtkImageOverwrite.h>
#include <mitkExtractSliceFilter.h>
#include <mitkLabelSetImage.h>
#include <mitkImageReadAccessor.h>
#include <mitkImageTimeSelector.h>
#include <mitkImageToContourModelFilter.h>
#include <mitkContourUtils.h>
#include <itkCommand.h>

#include <QCheckBox>
#include <QPushButton>
#include <QMenu>
#include <QCursor>
#include <QVBoxLayout>
#include <QMessageBox>

QmitkSlicesInterpolator::QmitkSlicesInterpolator(QWidget* parent, const char*  /*name*/)
  :QWidget(parent),
//    ACTION_TO_SLICEDIMENSION( createActionToSliceDimension() ),
    m_SliceInterpolatorController( mitk::SegmentationInterpolationController::New() ),
    m_SurfaceInterpolator(mitk::SurfaceInterpolationController::GetInstance()),
    m_ToolManager(NULL),
    m_Initialized(false),
    m_LastSNC(0),
    m_LastSliceIndex(0),
    m_2DInterpolationEnabled(false),
    m_3DInterpolationEnabled(false)
{
  m_GroupBoxEnableExclusiveInterpolationMode = new QGroupBox("Interpolation", this);
  m_GroupBoxEnableExclusiveInterpolationMode->setCheckable(true);
  m_GroupBoxEnableExclusiveInterpolationMode->setChecked(false);

  QVBoxLayout* vboxLayout = new QVBoxLayout(m_GroupBoxEnableExclusiveInterpolationMode);

  m_CmbInterpolation = new QComboBox(m_GroupBoxEnableExclusiveInterpolationMode);
  m_CmbInterpolation->addItem("2-Dimensional");
  m_CmbInterpolation->addItem("3-Dimensional");
  vboxLayout->addWidget(m_CmbInterpolation);

  m_BtnApply2D = new QPushButton("Apply", m_GroupBoxEnableExclusiveInterpolationMode);
  vboxLayout->addWidget(m_BtnApply2D);

  m_BtnApplyForAllSlices2D = new QPushButton("Apply for all slices", m_GroupBoxEnableExclusiveInterpolationMode);
  vboxLayout->addWidget(m_BtnApplyForAllSlices2D);

  m_BtnApply3D = new QPushButton("Apply", m_GroupBoxEnableExclusiveInterpolationMode);
  vboxLayout->addWidget(m_BtnApply3D);

  m_ChkShowPositionNodes = new QCheckBox("Show Position Nodes", m_GroupBoxEnableExclusiveInterpolationMode);
  vboxLayout->addWidget(m_ChkShowPositionNodes);

  this->HideAllInterpolationControls();

  connect(m_GroupBoxEnableExclusiveInterpolationMode, SIGNAL(toggled(bool)), this, SLOT(ActivateInterpolation(bool)));
  connect(m_CmbInterpolation, SIGNAL(currentIndexChanged(int)), this, SLOT(OnInterpolationMethodChanged(int)));
  connect(m_BtnApply2D, SIGNAL(clicked()), this, SLOT(OnAcceptInterpolationClicked()));
  connect(m_BtnApplyForAllSlices2D, SIGNAL(clicked()), this, SLOT(OnAcceptAllInterpolationsClicked()));
  connect(m_BtnApply3D, SIGNAL(clicked()), this, SLOT(OnAccept3DInterpolationClicked()));
  connect(m_ChkShowPositionNodes, SIGNAL(toggled(bool)), this, SLOT(OnShowMarkers(bool)));
  connect(m_ChkShowPositionNodes, SIGNAL(toggled(bool)), this, SIGNAL(SignalShowMarkerNodes(bool)));

  QHBoxLayout* layout = new QHBoxLayout(this);
  layout->addWidget(m_GroupBoxEnableExclusiveInterpolationMode);
  this->setLayout(layout);

  itk::ReceptorMemberCommand<QmitkSlicesInterpolator>::Pointer command = itk::ReceptorMemberCommand<QmitkSlicesInterpolator>::New();
  command->SetCallbackFunction( this, &QmitkSlicesInterpolator::OnSliceInterpolationInfoChanged );
  m_InterpolationInfoChangedObserverTag = m_SliceInterpolatorController->AddObserver( itk::ModifiedEvent(), command );

  itk::ReceptorMemberCommand<QmitkSlicesInterpolator>::Pointer command2 = itk::ReceptorMemberCommand<QmitkSlicesInterpolator>::New();
  command2->SetCallbackFunction( this, &QmitkSlicesInterpolator::OnSurfaceInterpolationInfoChanged );
  m_SurfaceInterpolationInfoChangedObserverTag = m_SurfaceInterpolator->AddObserver( itk::ModifiedEvent(), command2 );

  // feedback node and its visualization properties
  m_FeedbackNode = mitk::DataNode::New();
  m_FeedbackContour = mitk::ContourModel::New();
  m_FeedbackNode->SetData( m_FeedbackContour );
  m_FeedbackNode->SetName( "Interpolation feedback" );
  m_FeedbackNode->SetProperty( "helper object", mitk::BoolProperty::New(false) );
  m_FeedbackNode->SetProperty( "contour.width", mitk::FloatProperty::New( 3.0 ) );

  m_InterpolatedSurfaceNode = mitk::DataNode::New();
  m_InterpolatedSurfaceNode->SetName( "Surface Interpolation feedback" );
  m_InterpolatedSurfaceNode->SetProperty( "color", mitk::ColorProperty::New(255.0,255.0,0.0) );
  m_InterpolatedSurfaceNode->SetProperty( "opacity", mitk::FloatProperty::New(0.5) );
  m_InterpolatedSurfaceNode->SetProperty( "includeInBoundingBox", mitk::BoolProperty::New(false));
  m_InterpolatedSurfaceNode->SetProperty( "helper object", mitk::BoolProperty::New(true) );
  m_InterpolatedSurfaceNode->SetVisibility(false);

  m_3DContourNode = mitk::DataNode::New();
  m_3DContourNode->SetName( "Drawn Contours" );
  m_3DContourNode->SetProperty( "color", mitk::ColorProperty::New(0.0, 0.0, 0.0) );
  m_3DContourNode->SetProperty( "helper object", mitk::BoolProperty::New(true));
  m_3DContourNode->SetProperty( "material.representation", mitk::VtkRepresentationProperty::New(VTK_WIREFRAME));
  m_3DContourNode->SetProperty( "material.wireframeLineWidth", mitk::FloatProperty::New(2.0f));
  m_3DContourNode->SetProperty( "3DContourContainer", mitk::BoolProperty::New(true));
  m_3DContourNode->SetProperty( "includeInBoundingBox", mitk::BoolProperty::New(false));
  m_3DContourNode->SetVisibility(false, mitk::BaseRenderer::GetInstance( mitk::BaseRenderer::GetRenderWindowByName("stdmulti.widget1")));
  m_3DContourNode->SetVisibility(false, mitk::BaseRenderer::GetInstance( mitk::BaseRenderer::GetRenderWindowByName("stdmulti.widget2")));
  m_3DContourNode->SetVisibility(false, mitk::BaseRenderer::GetInstance( mitk::BaseRenderer::GetRenderWindowByName("stdmulti.widget3")));
  m_3DContourNode->SetVisibility(false, mitk::BaseRenderer::GetInstance( mitk::BaseRenderer::GetRenderWindowByName("stdmulti.widget4")));

  QWidget::setContentsMargins(0, 0, 0, 0);
  if ( QWidget::layout() != NULL )
  {
    QWidget::layout()->setContentsMargins(0, 0, 0, 0);
  }

  //For running 3D Interpolation in background
  // create a QFuture and a QFutureWatcher

  connect(&m_Watcher, SIGNAL(started()), this, SLOT(StartUpdateInterpolationTimer()));
  connect(&m_Watcher, SIGNAL(finished()), this, SLOT(OnSurfaceInterpolationFinished()));
  connect(&m_Watcher, SIGNAL(finished()), this, SLOT(StopUpdateInterpolationTimer()));
  m_Timer = new QTimer(this);
  connect(m_Timer, SIGNAL(timeout()), this, SLOT(ChangeSurfaceColor()));
}

const QmitkSlicesInterpolator::ActionToSliceDimensionMapType QmitkSlicesInterpolator::CreateActionToSliceDimension()
{
  ActionToSliceDimensionMapType actionToSliceDimension;
  foreach(mitk::SliceNavigationController* slicer, m_ControllerToDeleteObserverTag.keys())
  {
    std::string name = slicer->GetRenderer()->GetName();
    if (name == "stdmulti.widget1")
      name = "Axial (red window)";
    else if (name == "stdmulti.widget2")
      name = "Sagittal (green window)";
    else if (name == "stdmulti.widget3")
      name = "Coronal (blue window)";
    actionToSliceDimension[new QAction(QString::fromStdString(name),0)] = slicer;
  }

  return actionToSliceDimension;
}

void QmitkSlicesInterpolator::SetDataStorage( mitk::DataStorage::Pointer storage )
{
  m_DataStorage = storage;
  m_SurfaceInterpolator->SetDataStorage(storage);
}

mitk::DataStorage* QmitkSlicesInterpolator::GetDataStorage()
{
  if ( m_DataStorage.IsNotNull() )
  {
    return m_DataStorage;
  }
  else
  {
    return NULL;
  }
}


void QmitkSlicesInterpolator::Initialize(mitk::ToolManager* toolManager, const QList<mitk::SliceNavigationController *> &controllers)
{
  Q_ASSERT(!controllers.empty());

  if (!toolManager) return;

  if (m_Initialized)
  {
    // remove old observers
    this->Uninitialize();
  }

  m_ToolManager = toolManager;

  // set enabled only if a segmentation is selected
  mitk::DataNode* node = m_ToolManager->GetWorkingData(0);
  QWidget::setEnabled( node != NULL );

  // react whenever the set of selected segmentation changes
  m_ToolManager->WorkingDataChanged += mitk::MessageDelegate<QmitkSlicesInterpolator>( this, &QmitkSlicesInterpolator::OnToolManagerWorkingDataModified );
  m_ToolManager->ReferenceDataChanged += mitk::MessageDelegate<QmitkSlicesInterpolator>( this, &QmitkSlicesInterpolator::OnToolManagerReferenceDataModified );

  // connect to the slice navigation controller. after each change, call the interpolator
  foreach(mitk::SliceNavigationController* slicer, controllers)
  {
    //Has to be initialized
    m_LastSNC = slicer;

    m_TimeStep.insert(slicer, slicer->GetTime()->GetPos());

    itk::MemberCommand<QmitkSlicesInterpolator>::Pointer deleteCommand = itk::MemberCommand<QmitkSlicesInterpolator>::New();
    deleteCommand->SetCallbackFunction( this, &QmitkSlicesInterpolator::OnSliceNavigationControllerDeleted);
    m_ControllerToDeleteObserverTag.insert(slicer, slicer->AddObserver(itk::DeleteEvent(), deleteCommand));

    itk::MemberCommand<QmitkSlicesInterpolator>::Pointer timeChangedCommand = itk::MemberCommand<QmitkSlicesInterpolator>::New();
    timeChangedCommand->SetCallbackFunction( this, &QmitkSlicesInterpolator::OnTimeChanged);
    m_ControllerToTimeObserverTag.insert(slicer, slicer->AddObserver(mitk::SliceNavigationController::TimeSlicedGeometryEvent(NULL,0), timeChangedCommand));

    itk::MemberCommand<QmitkSlicesInterpolator>::Pointer sliceChangedCommand = itk::MemberCommand<QmitkSlicesInterpolator>::New();
    sliceChangedCommand->SetCallbackFunction( this, &QmitkSlicesInterpolator::OnSliceChanged);
    m_ControllerToSliceObserverTag.insert(slicer, slicer->AddObserver(mitk::SliceNavigationController::GeometrySliceEvent(NULL,0), sliceChangedCommand));
  }
  m_ActionToSliceDimensionMap = this->CreateActionToSliceDimension();

  m_Initialized = true;
}

void QmitkSlicesInterpolator::Uninitialize()
{
  if (m_ToolManager.IsNotNull())
  {
    m_ToolManager->WorkingDataChanged -= mitk::MessageDelegate<QmitkSlicesInterpolator>(this, &QmitkSlicesInterpolator::OnToolManagerWorkingDataModified);
    m_ToolManager->ReferenceDataChanged -= mitk::MessageDelegate<QmitkSlicesInterpolator>(this, &QmitkSlicesInterpolator::OnToolManagerReferenceDataModified);
  }

  foreach(mitk::SliceNavigationController* slicer, m_ControllerToSliceObserverTag.keys())
  {
    slicer->RemoveObserver(m_ControllerToDeleteObserverTag.take(slicer));
    slicer->RemoveObserver(m_ControllerToTimeObserverTag.take(slicer));
    slicer->RemoveObserver(m_ControllerToSliceObserverTag.take(slicer));
  }

  m_ActionToSliceDimensionMap.clear();

  m_ToolManager = NULL;

  m_Initialized = false;
}

QmitkSlicesInterpolator::~QmitkSlicesInterpolator()
{
  if (m_Initialized)
  {
    // remove old observers
    this->Uninitialize();
  }

  if ( m_DataStorage.IsNotNull() )
  {
    if(m_DataStorage->Exists(m_3DContourNode))
        m_DataStorage->Remove(m_3DContourNode);

    if(m_DataStorage->Exists(m_InterpolatedSurfaceNode))
        m_DataStorage->Remove(m_InterpolatedSurfaceNode);
  }

  // remove observer
  m_SliceInterpolatorController->RemoveObserver( m_InterpolationInfoChangedObserverTag );
  m_SurfaceInterpolator->RemoveObserver( m_SurfaceInterpolationInfoChangedObserverTag );

  delete m_Timer;
}

void QmitkSlicesInterpolator::HideAllInterpolationControls()
{
  this->Show2DInterpolationControls(false);
  this->Show3DInterpolationControls(false);
}

void QmitkSlicesInterpolator::Show2DInterpolationControls(bool show)
{
  m_BtnApply2D->setVisible(show);
  m_BtnApplyForAllSlices2D->setVisible(show);
}

void QmitkSlicesInterpolator::Show3DInterpolationControls(bool show)
{
  m_BtnApply3D->setVisible(show);
  m_ChkShowPositionNodes->setVisible(show);
}

void QmitkSlicesInterpolator::EnableInterpolation(bool enabled)
{
  m_GroupBoxEnableExclusiveInterpolationMode->setChecked(enabled);
}

void QmitkSlicesInterpolator::ActivateInterpolation(bool enabled)
{
  if (enabled)
  {
    if (m_3DInterpolationEnabled)
    {
      this->Show3DInterpolationControls(true);
      this->Show3DInterpolationResult(false);
    }
    else
    {
      this->Show2DInterpolationControls(true);
      this->Activate2DInterpolation(true);
    }
  }
  else
  {
    mitk::UndoController::GetCurrentUndoModel()->Clear();
    this->HideAllInterpolationControls();
    this->Activate2DInterpolation(false);
    this->Activate3DInterpolation(false);
    this->Show3DInterpolationResult(false);
  }
}

void QmitkSlicesInterpolator::OnInterpolationMethodChanged(int index)
{
  switch(index)
  {
    default:
    case 0: // 2D
      m_GroupBoxEnableExclusiveInterpolationMode->setTitle("2D Interpolation (Enabled)");
      this->HideAllInterpolationControls();
      this->Show2DInterpolationControls(true);
      this->Activate2DInterpolation(true);
      this->Activate3DInterpolation(false);
      break;

    case 1: // 3D
      m_GroupBoxEnableExclusiveInterpolationMode->setTitle("3D Interpolation (Enabled)");
      this->HideAllInterpolationControls();
      this->Show3DInterpolationControls(true);
      this->Activate2DInterpolation(false);
      this->Activate3DInterpolation(true);
      break;
  }
}

void QmitkSlicesInterpolator::OnShowMarkers(bool state)
{
  mitk::DataStorage::SetOfObjects::ConstPointer allContourMarkers =
    m_DataStorage->GetSubset(mitk::NodePredicateProperty::New("isContourMarker", mitk::BoolProperty::New(true)));

  for (mitk::DataStorage::SetOfObjects::ConstIterator it = allContourMarkers->Begin(); it != allContourMarkers->End(); ++it)
  {
    it->Value()->SetProperty("helper object", mitk::BoolProperty::New(!state));
  }
}

void QmitkSlicesInterpolator::OnWorkingImageModified(const itk::EventObject&)
{
  this->m_SliceInterpolatorController->BuildLabelCount();
  this->UpdateVisibleSuggestion();
}

void QmitkSlicesInterpolator::OnToolManagerWorkingDataModified()
{
  mitk::DataNode* workingNode = this->m_ToolManager->GetWorkingData(0);
  if (!workingNode) return;

  mitk::LabelSetImage::Pointer newImage = dynamic_cast< mitk::LabelSetImage* >( workingNode->GetData() );
  if ( newImage.IsNull() ) return;

  if (newImage->GetDimension() > 4 || newImage->GetDimension() < 3)
  {
    MITK_ERROR << "slices interpolator needs a 3D or 3D+t segmentation, not 2D.";
    return;
  }

  if (m_WorkingImage != newImage)
  {
    if (m_WorkingImage.IsNotNull())
      m_WorkingImage->RemoveObserver( m_WorkingImageObserverID );

    m_WorkingImage = newImage;

    // observe Modified() event of image
    itk::ReceptorMemberCommand<QmitkSlicesInterpolator>::Pointer command = itk::ReceptorMemberCommand<QmitkSlicesInterpolator>::New();
    command->SetCallbackFunction( this, &QmitkSlicesInterpolator::OnWorkingImageModified );
    m_WorkingImageObserverID = m_WorkingImage->AddObserver( itk::ModifiedEvent(), command );
  }

  m_SliceInterpolatorController->SetWorkingImage( m_WorkingImage );
  this->UpdateVisibleSuggestion();
}

void QmitkSlicesInterpolator::OnToolManagerReferenceDataModified()
{
/*
  if (m_2DInterpolationEnabled)
  {
    this->Activate2DInterpolation( true ); // re-initialize if needed
  }
  if (m_3DInterpolationEnabled)
  {
    this->Show3DInterpolationResult(false);
  }
*/
}

void QmitkSlicesInterpolator::OnTimeChanged(itk::Object* sender, const itk::EventObject& e)
{
  //Check if we really have a GeometryTimeEvent
  if (!dynamic_cast<const mitk::SliceNavigationController::GeometryTimeEvent*>(&e))
    return;

  mitk::SliceNavigationController* slicer = dynamic_cast<mitk::SliceNavigationController*>(sender);
  Q_ASSERT(slicer);

  m_TimeStep[slicer]/* = event.GetPos()*/;

  //TODO Macht das hier wirklich Sinn????
  if (m_LastSNC == slicer)
  {
    slicer->SendSlice();//will trigger a new interpolation
  }
}

void QmitkSlicesInterpolator::OnSliceChanged(itk::Object *sender, const itk::EventObject &e)
{
  //Check whether we really have a GeometrySliceEvent
  if (!dynamic_cast<const mitk::SliceNavigationController::GeometrySliceEvent*>(&e))
    return;

  mitk::SliceNavigationController* slicer = dynamic_cast<mitk::SliceNavigationController*>(sender);

  if (this->TranslateAndInterpolateChangedSlice(e, slicer))
  {
    slicer->GetRenderer()->RequestUpdate();
  }
}

bool QmitkSlicesInterpolator::TranslateAndInterpolateChangedSlice(const itk::EventObject& e, mitk::SliceNavigationController* slicer)
{
  if (!m_2DInterpolationEnabled) return false;

  const mitk::SliceNavigationController::GeometrySliceEvent& event = dynamic_cast<const mitk::SliceNavigationController::GeometrySliceEvent&>(e);

  mitk::TimeSlicedGeometry* tsg = event.GetTimeSlicedGeometry();
  if (tsg && m_TimeStep.contains(slicer))
  {
    mitk::SlicedGeometry3D* slicedGeometry = dynamic_cast<mitk::SlicedGeometry3D*>(tsg->GetGeometry3D(m_TimeStep[slicer]));
    if (slicedGeometry)
    {
      m_LastSNC = slicer;
      mitk::PlaneGeometry* plane = dynamic_cast<mitk::PlaneGeometry*>(slicedGeometry->GetGeometry2D( event.GetPos() ));
      if (plane != NULL)
      {
        this->Interpolate( plane, m_TimeStep[slicer], slicer );
        return true;
      }
    }
  }

  return false;
}

void QmitkSlicesInterpolator::Interpolate( mitk::PlaneGeometry* plane, unsigned int timeStep, mitk::SliceNavigationController* slicer )
{
  int clickedSliceDimension(-1);
  int clickedSliceIndex(-1);

  // calculate real slice position, i.e. slice of the image and not slice of the TimeSlicedGeometry
  // see if timestep is needed here
  mitk::SegTool2D::DetermineAffectedImageSlice( m_WorkingImage, plane, clickedSliceDimension, clickedSliceIndex );

  mitk::Image::Pointer interpolation = m_SliceInterpolatorController->Interpolate( clickedSliceDimension, clickedSliceIndex, plane, timeStep );

  if (interpolation.IsNotNull())
  {
    mitk::ImageToContourModelFilter::Pointer contourExtractor = mitk::ImageToContourModelFilter::New();
    contourExtractor->SetInput(interpolation);
    contourExtractor->SetUseProgressBar(false);

    try
    {
      contourExtractor->Update();
    }
    catch ( itk::ExceptionObject & excep )
    {
      MITK_ERROR << "Exception caught: " << excep.GetDescription();
      return;
    }

    int numberOfContours = contourExtractor->GetNumberOfIndexedOutputs();

    m_FeedbackContour = contourExtractor->GetOutput(0);
    m_FeedbackContour->DisconnectPipeline();

    m_FeedbackNode->SetData( m_FeedbackContour );

    const mitk::Color& color = m_WorkingImage->GetActiveLabelColor();
    m_FeedbackNode->SetProperty("contour.color", mitk::ColorProperty::New(color));
  }
  else
  {
    m_FeedbackContour->Clear(timeStep);
  }

  m_LastSNC = slicer;
  m_LastSliceIndex = clickedSliceIndex;
}

void QmitkSlicesInterpolator::OnSurfaceInterpolationFinished()
{
  mitk::Surface::Pointer interpolatedSurface = m_SurfaceInterpolator->GetInterpolationResult();
  mitk::DataNode* workingNode = m_ToolManager->GetWorkingData(0);

  if (interpolatedSurface.IsNotNull() && workingNode &&
     workingNode->IsVisible(mitk::BaseRenderer::GetInstance( mitk::BaseRenderer::GetRenderWindowByName("Coronal"))))
  {
    m_BtnApply3D->setEnabled(true);
    m_InterpolatedSurfaceNode->SetData(interpolatedSurface);
    m_3DContourNode->SetData(m_SurfaceInterpolator->GetContoursAsSurface());

    this->Show3DInterpolationResult(true);

    if( !m_DataStorage->Exists(m_InterpolatedSurfaceNode) && !m_DataStorage->Exists(m_3DContourNode))
    {
      m_DataStorage->Add(m_3DContourNode);
      m_DataStorage->Add(m_InterpolatedSurfaceNode);
    }
  }
  else if (interpolatedSurface.IsNull())
  {
    m_BtnApply3D->setEnabled(false);

    if (m_DataStorage->Exists(m_InterpolatedSurfaceNode))
    {
      this->Show3DInterpolationResult(false);
    }
  }

  foreach (mitk::SliceNavigationController* slicer, m_ControllerToTimeObserverTag.keys())
  {
    slicer->GetRenderer()->RequestUpdate();
  }
}

mitk::Image::Pointer QmitkSlicesInterpolator::GetWorkingSlice()
{
  unsigned int timeStep = m_LastSNC->GetTime()->GetPos();

  const mitk::PlaneGeometry* planeGeometry = m_LastSNC->GetCurrentPlaneGeometry();
  if (!planeGeometry) return NULL;

  //Make sure that for reslicing and overwriting the same alogrithm is used. We can specify the mode of the vtk reslicer
  vtkSmartPointer<mitkVtkImageOverwrite> reslice = vtkSmartPointer<mitkVtkImageOverwrite>::New();
  //set to false to extract a slice
  reslice->SetOverwriteMode(false);
  reslice->Modified();

  //use ExtractSliceFilter with our specific vtkImageReslice for overwriting and extracting
  mitk::ExtractSliceFilter::Pointer extractor =  mitk::ExtractSliceFilter::New(reslice);
  extractor->SetInput( m_WorkingImage );
  extractor->SetTimeStep( timeStep );
  extractor->SetWorldGeometry( planeGeometry );
  extractor->SetVtkOutputRequest(false);
  extractor->SetResliceTransformByGeometry( m_WorkingImage->GetTimeSlicedGeometry()->GetGeometry3D( timeStep ) );

  extractor->Modified();

  try
  {
    extractor->Update();
  }
  catch ( itk::ExceptionObject & excep )
  {
    MITK_ERROR << "Exception caught: " << excep.GetDescription();
    return NULL;
  }

  mitk::Image::Pointer slice = extractor->GetOutput();

  //specify the undo operation with the non edited slice
  m_undoOperation = new mitk::DiffSliceOperation(
    m_WorkingImage, extractor->GetVtkOutput(), slice->GetGeometry(), timeStep, const_cast<mitk::PlaneGeometry*>(planeGeometry));

  slice->DisconnectPipeline();

  return slice;
}

void QmitkSlicesInterpolator::OnAcceptInterpolationClicked()
{
  if (m_WorkingImage.IsNotNull() && m_FeedbackNode->GetData())
  {
    mitk::Image::Pointer slice = this->GetWorkingSlice();
    if (slice.IsNull()) return;

    unsigned int timeStep = m_LastSNC->GetTime()->GetPos();

    const mitk::PlaneGeometry* planeGeometry = m_LastSNC->GetCurrentPlaneGeometry();
    if (!planeGeometry) return;

    mitk::ContourModel::Pointer projectedContour = mitk::ContourModel::New();
    mitk::ContourUtils::ProjectContourTo2DSlice( slice, m_FeedbackContour, projectedContour, timeStep );
    if (projectedContour.IsNull()) return;

    mitk::ContourUtils::FillContourInSlice( projectedContour, slice, m_WorkingImage->GetActiveLabelIndex(), timeStep );

    //Make sure that for reslicing and overwriting the same alogrithm is used. We can specify the mode of the vtk reslicer
    vtkSmartPointer<mitkVtkImageOverwrite> overwrite = vtkSmartPointer<mitkVtkImageOverwrite>::New();
    overwrite->SetInputSlice(slice->GetVtkImageData());
    //set overwrite mode to true to write back to the image volume
    overwrite->SetOverwriteMode(true);
    overwrite->Modified();

    mitk::ExtractSliceFilter::Pointer extractor =  mitk::ExtractSliceFilter::New(overwrite);
    extractor->SetInput( m_WorkingImage );
    extractor->SetTimeStep( timeStep );
    extractor->SetWorldGeometry( planeGeometry );
    extractor->SetVtkOutputRequest(true);
    extractor->SetResliceTransformByGeometry( m_WorkingImage->GetTimeSlicedGeometry()->GetGeometry3D( timeStep ) );

    extractor->Modified();

    try
    {
      extractor->Update();
    }
    catch ( itk::ExceptionObject & excep )
    {
      MITK_ERROR << "Exception caught: " << excep.GetDescription();
      return;
    }

    //the image was modified within the pipeline, but not marked so
    m_WorkingImage->Modified();
    m_WorkingImage->GetVtkImageData()->Modified();

    //specify the undo operation with the edited slice
    m_doOperation = new mitk::DiffSliceOperation(
      m_WorkingImage, extractor->GetVtkOutput(),slice->GetGeometry(), timeStep, const_cast<mitk::PlaneGeometry*>(planeGeometry));

    //create an operation event for the undo stack
    mitk::OperationEvent* undoStackItem = new mitk::OperationEvent(
      mitk::DiffSliceOperationApplier::GetInstance(), m_doOperation, m_undoOperation, "Segmentation" );

    //add it to the undo controller
    mitk::UndoController::GetCurrentUndoModel()->SetOperationEvent( undoStackItem );

    //clear the pointers as the operation are stored in the undocontroller and also deleted from there
    m_undoOperation = NULL;
    m_doOperation = NULL;

    m_FeedbackContour->Clear();

    mitk::RenderingManager::GetInstance()->RequestUpdateAll();
  }
}

void QmitkSlicesInterpolator::AcceptAllInterpolations(mitk::SliceNavigationController* slicer)
{
  /*
   * What exactly is done here:
   * 1. We create an empty diff image for the current segmentation
   * 2. All interpolated slices are written into the diff image
   * 3. Then the diffimage is applied to the original segmentation
   */
  if (m_WorkingImage.IsNotNull())
  {
/*
    //making interpolation separately undoable
    mitk::UndoStackItem::IncCurrObjectEventId();
    mitk::UndoStackItem::IncCurrGroupEventId();
    mitk::UndoStackItem::ExecuteIncrement();
*/
    /*
    mitk::Image::Pointer image3D = m_WorkingImage;
    unsigned int timeStep = slicer->GetTime()->GetPos();
    if (m_WorkingImage->GetDimension() == 4)
    {
      mitk::ImageTimeSelector::Pointer timeSelector = mitk::ImageTimeSelector::New();
      timeSelector->SetInput( m_WorkingImage );
      timeSelector->SetTimeNr( timeStep );
      timeSelector->Update();
      image3D = timeSelector->GetOutput();
    }
    */
    /*
    // create a empty diff image for the undo operation
    mitk::Image::Pointer diffImage = mitk::Image::New();
    diffImage->Initialize( image3D );

    // Set all pixels to zero
    mitk::PixelType pixelType( mitk::MakeScalarPixelType<unsigned char>()  );
    memset( diffImage->GetData(), 0, (pixelType.GetBpe() >> 3) * diffImage->GetDimension(0) * diffImage->GetDimension(1) * diffImage->GetDimension(2) );
*/
    // Since we need to shift the plane it must be clone so that the original plane isn't altered
    mitk::PlaneGeometry::Pointer reslicePlane = slicer->GetCurrentPlaneGeometry()->Clone();
    unsigned int timeStep = slicer->GetTime()->GetPos();
    int sliceDimension(-1);
    int sliceIndex(-1);
    mitk::SegTool2D::DetermineAffectedImageSlice( m_WorkingImage, reslicePlane, sliceDimension, sliceIndex );

    unsigned int zslices = m_WorkingImage->GetDimension( sliceDimension );
    mitk::ProgressBar::GetInstance()->AddStepsToDo(zslices);

    mitk::Point3D origin = reslicePlane->GetOrigin();
    unsigned int totalChangedSlices(0);

    for (unsigned int sliceIndex = 0; sliceIndex < zslices; ++sliceIndex)
    {
      // Transforming the current origin of the reslice plane
      // so that it matches the one of the next slice
      m_WorkingImage->GetSlicedGeometry()->WorldToIndex(origin, origin);
      origin[sliceDimension] = sliceIndex;
      m_WorkingImage->GetSlicedGeometry()->IndexToWorld(origin, origin);
      reslicePlane->SetOrigin(origin);
      //Set the slice as 'input'
      mitk::Image::Pointer interpolation = m_SliceInterpolatorController->Interpolate( sliceDimension, sliceIndex, reslicePlane, timeStep );

      if (interpolation.IsNotNull()) // we don't check if interpolation is necessary/sensible - but m_InterpolatorController does
      {
/*
        //Setting up the reslicing pipeline which allows us to write the interpolation results back into
        //the image volume
        vtkSmartPointer<mitkVtkImageOverwrite> reslice = vtkSmartPointer<mitkVtkImageOverwrite>::New();

        //set overwrite mode to true to write back to the image volume
        reslice->SetInputSlice(interpolation->GetVtkImageData());
        reslice->SetOverwriteMode(true);
        reslice->Modified();

        mitk::ExtractSliceFilter::Pointer diffslicewriter =  mitk::ExtractSliceFilter::New(reslice);
        diffslicewriter->SetInput( m_WorkingImage ); //diffImage );
        diffslicewriter->SetTimeStep( timeStep );
        diffslicewriter->SetWorldGeometry(reslicePlane);
        diffslicewriter->SetVtkOutputRequest(true);
        diffslicewriter->SetResliceTransformByGeometry( m_WorkingImage->GetTimeSlicedGeometry()->GetGeometry3D( timeStep ) );

        diffslicewriter->Modified();
        diffslicewriter->Update();
*/
        //Make sure that for reslicing and overwriting the same alogrithm is used. We can specify the mode of the vtk reslicer
        vtkSmartPointer<mitkVtkImageOverwrite> overwrite = vtkSmartPointer<mitkVtkImageOverwrite>::New();
        overwrite->SetInputSlice(interpolation->GetVtkImageData());
        //set overwrite mode to true to write back to the image volume
        overwrite->SetOverwriteMode(true);
        overwrite->Modified();

        mitk::ExtractSliceFilter::Pointer extractor =  mitk::ExtractSliceFilter::New(overwrite);
        extractor->SetInput( m_WorkingImage );
        extractor->SetTimeStep( timeStep );
        extractor->SetWorldGeometry( reslicePlane );
        extractor->SetVtkOutputRequest(true);
        extractor->SetResliceTransformByGeometry( m_WorkingImage->GetTimeSlicedGeometry()->GetGeometry3D( timeStep ) );

        extractor->Modified();

        try
        {
          extractor->Update();
        }
        catch ( itk::ExceptionObject & excep )
        {
          MITK_ERROR << "Exception caught: " << excep.GetDescription();
          return;
        }

        //the image was modified within the pipeline, but not marked so
        m_WorkingImage->Modified();
        m_WorkingImage->GetVtkImageData()->Modified();

        ++totalChangedSlices;

        mitk::RenderingManager::GetInstance()->RequestUpdateAll(mitk::RenderingManager::REQUEST_UPDATE_2DWINDOWS);
      }
      mitk::ProgressBar::GetInstance()->Progress();
    }
/*
    if (totalChangedSlices > 0)
    {
        // create do/undo operations
        mitk::ApplyDiffImageOperation* doOp = new mitk::ApplyDiffImageOperation( mitk::OpTEST, m_WorkingImage, diffImage, timeStep );
        mitk::ApplyDiffImageOperation* undoOp = new mitk::ApplyDiffImageOperation( mitk::OpTEST, m_WorkingImage, diffImage, timeStep );
        undoOp->SetFactor( -1.0 );
        std::stringstream comment;
        comment << "Accept all interpolations (" << totalChangedSlices << ")";
        mitk::OperationEvent* undoStackItem = new mitk::OperationEvent( mitk::DiffImageApplier::GetInstanceForUndo(), doOp, undoOp, comment.str() );
        mitk::UndoController::GetCurrentUndoModel()->SetOperationEvent( undoStackItem );

        // acutally apply the changes here to the original image
        mitk::DiffImageApplier::GetInstanceForUndo()->ExecuteOperation( doOp );
    }
*/
//    m_FeedbackNode->SetData(NULL);
    mitk::RenderingManager::GetInstance()->RequestUpdateAll();
  }
}

void QmitkSlicesInterpolator::FinishInterpolation(mitk::SliceNavigationController* slicer)
{
  //this redirect is for calling from outside

  if (slicer == NULL)
    this->OnAcceptAllInterpolationsClicked();
  else
    this->AcceptAllInterpolations( slicer );
}

void QmitkSlicesInterpolator::OnAcceptAllInterpolationsClicked()
{
  QMenu orientationPopup(this);
  std::map<QAction*, mitk::SliceNavigationController*>::const_iterator it;
  for(it = m_ActionToSliceDimensionMap.begin(); it != m_ActionToSliceDimensionMap.end(); it++)
    orientationPopup.addAction(it->first);

  connect( &orientationPopup, SIGNAL(triggered(QAction*)), this, SLOT(OnAcceptAllPopupActivated(QAction*)) );

  orientationPopup.exec( QCursor::pos() );
}

void QmitkSlicesInterpolator::OnAccept3DInterpolationClicked()
{
  if (m_InterpolatedSurfaceNode.IsNotNull() && m_InterpolatedSurfaceNode->GetData())
  {
    mitk::SurfaceToImageFilter::Pointer s2iFilter = mitk::SurfaceToImageFilter::New();
    s2iFilter->MakeOutputBinaryOn();
    s2iFilter->SetInput(dynamic_cast<mitk::Surface*>(m_InterpolatedSurfaceNode->GetData()));

    // check if ToolManager holds valid ReferenceData
    if (m_ToolManager->GetReferenceData(0) == NULL || m_ToolManager->GetWorkingData(0) == NULL)
    {
        return;
    }
    s2iFilter->SetImage(dynamic_cast<mitk::Image*>(m_ToolManager->GetReferenceData(0)->GetData()));
    s2iFilter->Update();

    mitk::DataNode* segmentationNode = m_ToolManager->GetWorkingData(0);
    segmentationNode->SetData(s2iFilter->GetOutput());
    m_CmbInterpolation->setCurrentIndex(0);
    mitk::RenderingManager::GetInstance()->RequestUpdateAll();
    this->Show3DInterpolationResult(false);
  }
}

void QmitkSlicesInterpolator::OnAcceptAllPopupActivated(QAction* action)
{
  try
  {
    ActionToSliceDimensionMapType::const_iterator iter = m_ActionToSliceDimensionMap.find( action );
    if (iter != m_ActionToSliceDimensionMap.end())
    {
      mitk::SliceNavigationController* slicer = iter->second;
      this->AcceptAllInterpolations( slicer );
    }
  }
  catch(...)
  {
    /* Showing message box with possible memory error */
    QMessageBox errorInfo;
    errorInfo.setWindowTitle("Interpolation Process");
    errorInfo.setIcon(QMessageBox::Critical);
    errorInfo.setText("An error occurred during interpolation. Possible cause: Not enough memory!");
    errorInfo.exec();

    //additional error message
    MITK_ERROR << "Ill construction in " __FILE__ " l. " << __LINE__ ;
  }
}

void QmitkSlicesInterpolator::Activate2DInterpolation(bool on)
{
  if (m_WorkingImage.IsNull()) return;

  m_2DInterpolationEnabled = on;

  if ( m_DataStorage.IsNotNull() )
  {
    if (on && !m_DataStorage->Exists(m_FeedbackNode))
    {
      m_DataStorage->Add( m_FeedbackNode );
    }
    else if (!on && m_DataStorage->Exists(m_FeedbackNode))
    {
      m_DataStorage->Remove( m_FeedbackNode );
    }
  }

  if (on)
  {
    this->UpdateVisibleSuggestion();
    return;
  }

/*
  mitk::DataNode* referenceNode = m_ToolManager->GetReferenceData(0);
  if (referenceNode)
  {
    mitk::Image* referenceImage = dynamic_cast<mitk::Image*>(referenceNode->GetData());
    m_SliceInterpolatorController->SetReferenceImage( referenceImage ); // may be NULL
  }
*/
}

void QmitkSlicesInterpolator::Run3DInterpolation()
{
  m_SurfaceInterpolator->Interpolate();
}

void QmitkSlicesInterpolator::StartUpdateInterpolationTimer()
{
  m_Timer->start(500);
}

void QmitkSlicesInterpolator::StopUpdateInterpolationTimer()
{
  m_Timer->stop();
  m_InterpolatedSurfaceNode->SetProperty("color", mitk::ColorProperty::New(255.0,255.0,0.0));
  mitk::RenderingManager::GetInstance()->RequestUpdate(mitk::BaseRenderer::GetInstance( mitk::BaseRenderer::GetRenderWindowByName("3D"))->GetRenderWindow());
}

void QmitkSlicesInterpolator::ChangeSurfaceColor()
{
  float currentColor[3];
  m_InterpolatedSurfaceNode->GetColor(currentColor);

  float yellow[3] = {255.0,255.0,0.0};

  if( currentColor[2] == yellow[2])
  {
    m_InterpolatedSurfaceNode->SetProperty("color", mitk::ColorProperty::New(255.0,255.0,255.0));
  }
  else
  {
    m_InterpolatedSurfaceNode->SetProperty("color", mitk::ColorProperty::New(yellow));
  }
  m_InterpolatedSurfaceNode->Update();
  mitk::RenderingManager::GetInstance()->RequestUpdate(mitk::BaseRenderer::GetInstance( mitk::BaseRenderer::GetRenderWindowByName("3D"))->GetRenderWindow());
}

void QmitkSlicesInterpolator::Activate3DInterpolation(bool on)
{
  m_3DInterpolationEnabled = on;

  try
  {
    if ( m_DataStorage.IsNotNull() && m_ToolManager && m_3DInterpolationEnabled)
    {
      mitk::DataNode* workingNode = m_ToolManager->GetWorkingData(0);

      if (workingNode)
      {
        bool isInterpolationResult(false);
        workingNode->GetBoolProperty("3DInterpolationResult",isInterpolationResult);

        if ((workingNode->IsSelected() &&
             workingNode->IsVisible(mitk::BaseRenderer::GetInstance( mitk::BaseRenderer::GetRenderWindowByName("Coronal")))) &&
            !isInterpolationResult && m_3DInterpolationEnabled)
        {
          int ret = QMessageBox::Yes;

          if (m_SurfaceInterpolator->EstimatePortionOfNeededMemory() > 0.5)
          {
            QMessageBox msgBox;
            msgBox.setText("Due to short handed system memory the 3D interpolation may be very slow!");
            msgBox.setInformativeText("Are you sure you want to activate the 3D interpolation?");
            msgBox.setStandardButtons(QMessageBox::No | QMessageBox::Yes);
            ret = msgBox.exec();
          }

          if (m_Watcher.isRunning())
            m_Watcher.waitForFinished();

          if (ret == QMessageBox::Yes)
          {
            m_Future = QtConcurrent::run(this, &QmitkSlicesInterpolator::Run3DInterpolation);
            m_Watcher.setFuture(m_Future);
          }
          else
          {
            m_CmbInterpolation->setCurrentIndex(0);
          }
        }
        else if (!m_3DInterpolationEnabled)
        {
          this->Show3DInterpolationResult(false);
          m_BtnApply3D->setEnabled(m_3DInterpolationEnabled);
        }
      }
      else
      {
        QWidget::setEnabled( false );
        m_ChkShowPositionNodes->setEnabled(m_3DInterpolationEnabled);
      }
    }
    if (!m_3DInterpolationEnabled)
    {
       this->Show3DInterpolationResult(false);
       m_BtnApply3D->setEnabled(m_3DInterpolationEnabled);
    }
  }
  catch(...)
  {
    MITK_ERROR<<"Error with 3D surface interpolation!";
  }
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}


void QmitkSlicesInterpolator::UpdateVisibleSuggestion()
{
  if (m_2DInterpolationEnabled && m_LastSNC)
  {
    // determine which one is the current view, try to do an initial interpolation
    mitk::BaseRenderer* renderer = m_LastSNC->GetRenderer();
    if (renderer && renderer->GetMapperID() == mitk::BaseRenderer::Standard2D)
    {
      const mitk::TimeSlicedGeometry* timeSlicedGeometry = dynamic_cast<const mitk::TimeSlicedGeometry*>( renderer->GetWorldGeometry() );
      if (timeSlicedGeometry)
      {
        mitk::SliceNavigationController::GeometrySliceEvent event( const_cast<mitk::TimeSlicedGeometry*>(timeSlicedGeometry), renderer->GetSlice() );
        this->TranslateAndInterpolateChangedSlice(event, m_LastSNC);
      }
    }
  }

  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}

void QmitkSlicesInterpolator::OnSliceInterpolationInfoChanged(const itk::EventObject& /*e*/)
{
  // something (e.g. undo) changed the interpolation info, we should refresh our display
  this->UpdateVisibleSuggestion();
}

void QmitkSlicesInterpolator::OnSurfaceInterpolationInfoChanged(const itk::EventObject& /*e*/)
{
  if(m_3DInterpolationEnabled)
  {
    if (m_Watcher.isRunning())
      m_Watcher.waitForFinished();
    m_Future = QtConcurrent::run(this, &QmitkSlicesInterpolator::Run3DInterpolation);
    m_Watcher.setFuture(m_Future);
  }
}

void QmitkSlicesInterpolator::SetCurrentContourListID()
{
  // New ContourList = hide current interpolation
  Show3DInterpolationResult(false);

  if ( m_DataStorage.IsNotNull() && m_ToolManager && m_LastSNC )
  {
    mitk::DataNode* workingNode = m_ToolManager->GetWorkingData(0);

    if (workingNode)
    {
      bool isInterpolationResult(false);
      workingNode->GetBoolProperty("3DInterpolationResult",isInterpolationResult);

      bool isVisible (workingNode->IsVisible(m_LastSNC->GetRenderer()));
      if (isVisible && !isInterpolationResult)
      {
        QWidget::setEnabled( true );

        //TODO Aufruf hier pruefen!
        mitk::Vector3D spacing = workingNode->GetData()->GetGeometry( m_LastSNC->GetTime()->GetPos() )->GetSpacing();
        double minSpacing (100);
        double maxSpacing (0);
        for (int i =0; i < 3; i++)
        {
          if (spacing[i] < minSpacing)
          {
            minSpacing = spacing[i];
          }
          else if (spacing[i] > maxSpacing)
          {
            maxSpacing = spacing[i];
          }
        }

        m_SurfaceInterpolator->SetSegmentationImage(dynamic_cast<mitk::Image*>(workingNode->GetData()));
        m_SurfaceInterpolator->SetMaxSpacing(maxSpacing);
        m_SurfaceInterpolator->SetMinSpacing(minSpacing);
        m_SurfaceInterpolator->SetDistanceImageVolume(50000);

        m_SurfaceInterpolator->SetCurrentSegmentationInterpolationList(dynamic_cast<mitk::Image*>(workingNode->GetData()));

        if (m_3DInterpolationEnabled)
        {
          if (m_Watcher.isRunning())
            m_Watcher.waitForFinished();
          m_Future = QtConcurrent::run(this, &QmitkSlicesInterpolator::Run3DInterpolation);
          m_Watcher.setFuture(m_Future);
        }
      }
    }
    else
    {
      QWidget::setEnabled(false);
    }
  }
}

void QmitkSlicesInterpolator::Show3DInterpolationResult(bool status)
{
   if (m_InterpolatedSurfaceNode.IsNotNull())
      m_InterpolatedSurfaceNode->SetVisibility(status);

   if (m_3DContourNode.IsNotNull())
      m_3DContourNode->SetVisibility(status, mitk::BaseRenderer::GetInstance( mitk::BaseRenderer::GetRenderWindowByName("stdmulti.widget4")));

   mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}

void QmitkSlicesInterpolator::OnSliceNavigationControllerDeleted(const itk::Object *sender, const itk::EventObject& /*e*/)
{
  //Don't know how to avoid const_cast here?!
  mitk::SliceNavigationController* slicer = dynamic_cast<mitk::SliceNavigationController*>(const_cast<itk::Object*>(sender));
  if (slicer)
  {
    m_ControllerToTimeObserverTag.remove(slicer);
    m_ControllerToSliceObserverTag.remove(slicer);
    m_ControllerToDeleteObserverTag.remove(slicer);
  }
}
