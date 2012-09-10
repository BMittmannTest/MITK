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
// Blueberry
#include <mitkIRenderingManager.h>
#include <mitkIRenderWindowPart.h>
#include <mitkILinkedRenderWindowPart.h>


// Qmitk
#include "QmitkToFUtilView.h"
#include <QmitkStdMultiWidget.h>
#include <QmitkTextOverlay.h>

// Qt
#include <QMessageBox>
#include <QString>

// MITK
#include <mitkBaseRenderer.h>
#include <mitkGlobalInteraction.h>
#include <mitkLookupTableProperty.h>
#include <mitkToFDistanceImageToPointSetFilter.h>
#include <mitkTransferFunction.h>
#include <mitkTransferFunctionProperty.h>

// VTK
#include <vtkCamera.h>

// ITK
#include <itkCommand.h>

const std::string QmitkToFUtilView::VIEW_ID = "org.mitk.views.tofutil";

QmitkToFUtilView::QmitkToFUtilView()
: QmitkAbstractView()
, m_Controls(NULL), m_MultiWidget( NULL )
, m_MitkDistanceImage(NULL), m_MitkAmplitudeImage(NULL), m_MitkIntensityImage(NULL), m_Surface(NULL)
, m_DistanceImageNode(NULL), m_AmplitudeImageNode(NULL), m_IntensityImageNode(NULL), m_RGBImageNode(NULL), m_SurfaceNode(NULL)
, m_ToFImageRecorder(NULL), m_ToFImageGrabber(NULL), m_ToFDistanceImageToSurfaceFilter(NULL), m_ToFCompositeFilter(NULL)
, m_SurfaceDisplayCount(0), m_2DDisplayCount(0)
, m_RealTimeClock(NULL)
, m_StepsForFramerate(100)
, m_2DTimeBefore(0.0)
, m_2DTimeAfter(0.0)
, m_VideoEnabled(false)
{
  this->m_Frametimer = new QTimer(this);

  this->m_ToFDistanceImageToSurfaceFilter = mitk::ToFDistanceImageToSurfaceFilter::New();
  this->m_ToFCompositeFilter = mitk::ToFCompositeFilter::New();
  this->m_ToFImageRecorder = mitk::ToFImageRecorder::New();
  this->m_ToFSurfaceVtkMapper3D = mitk::ToFSurfaceVtkMapper3D::New();
}

QmitkToFUtilView::~QmitkToFUtilView()
{
        OnToFCameraStopped();
        OnToFCameraDisconnected();
}

void QmitkToFUtilView::SetFocus()
{
  m_Controls->m_ToFConnectionWidget->setFocus();
}

void QmitkToFUtilView::CreateQtPartControl( QWidget *parent )
{
  // build up qt view, unless already done
  if ( !m_Controls )
  {
    // create GUI widgets from the Qt Designer's .ui file
    m_Controls = new Ui::QmitkToFUtilViewControls;
    m_Controls->setupUi( parent );

    connect(m_Frametimer, SIGNAL(timeout()), this, SLOT(OnUpdateCamera()));
    connect( (QObject*)(m_Controls->m_ToFConnectionWidget), SIGNAL(ToFCameraConnected()), this, SLOT(OnToFCameraConnected()) );
    connect( (QObject*)(m_Controls->m_ToFConnectionWidget), SIGNAL(ToFCameraDisconnected()), this, SLOT(OnToFCameraDisconnected()) );
    connect( (QObject*)(m_Controls->m_ToFConnectionWidget), SIGNAL(ToFCameraSelected(const QString)), this, SLOT(OnToFCameraSelected(const QString)) );
    connect( (QObject*)(m_Controls->m_ToFRecorderWidget), SIGNAL(ToFCameraStarted()), this, SLOT(OnToFCameraStarted()) );
    connect( (QObject*)(m_Controls->m_ToFRecorderWidget), SIGNAL(ToFCameraStopped()), this, SLOT(OnToFCameraStopped()) );
    connect( (QObject*)(m_Controls->m_ToFRecorderWidget), SIGNAL(RecordingStarted()), this, SLOT(OnToFCameraStopped()) );
    connect( (QObject*)(m_Controls->m_ToFRecorderWidget), SIGNAL(RecordingStopped()), this, SLOT(OnToFCameraStarted()) );
    connect( (QObject*)(m_Controls->m_TextureCheckBox), SIGNAL(toggled(bool)), this, SLOT(OnTextureCheckBoxChecked(bool)) );
    connect( (QObject*)(m_Controls->m_VideoTextureCheckBox), SIGNAL(toggled(bool)), this, SLOT(OnVideoTextureCheckBoxChecked(bool)) );
  }
}

void QmitkToFUtilView::Activated()
{
  //get the current RenderWindowPart or open a new one if there is none
  if(this->GetRenderWindowPart(OPEN))
  {
    mitk::ILinkedRenderWindowPart* linkedRenderWindowPart = dynamic_cast<mitk::ILinkedRenderWindowPart*>(this->GetRenderWindowPart());
    if(linkedRenderWindowPart == 0)
    {
      MITK_ERROR << "No linked StdMultiWidget avaiable!!!";
    }
    else
    {
      linkedRenderWindowPart->EnableSlicingPlanes(false);
    }
    GetRenderWindowPart()->GetRenderWindow("axial")->GetSliceNavigationController()->SetDefaultViewDirection(mitk::SliceNavigationController::Axial);
    GetRenderWindowPart()->GetRenderWindow("axial")->GetSliceNavigationController()->SliceLockedOn();
    GetRenderWindowPart()->GetRenderWindow("sagittal")->GetSliceNavigationController()->SetDefaultViewDirection(mitk::SliceNavigationController::Axial);
    GetRenderWindowPart()->GetRenderWindow("sagittal")->GetSliceNavigationController()->SliceLockedOn();
    GetRenderWindowPart()->GetRenderWindow("coronal")->GetSliceNavigationController()->SetDefaultViewDirection(mitk::SliceNavigationController::Axial);
    GetRenderWindowPart()->GetRenderWindow("coronal")->GetSliceNavigationController()->SliceLockedOn();

    this->GetRenderWindowPart()->GetRenderingManager()->InitializeViews();

    this->UseToFVisibilitySettings(true);

    m_Controls->m_ToFCompositeFilterWidget->SetToFCompositeFilter(this->m_ToFCompositeFilter);
    m_Controls->m_ToFCompositeFilterWidget->SetDataStorage(this->GetDataStorage());

    if (this->m_ToFImageGrabber.IsNull())
    {
      m_Controls->m_ToFRecorderWidget->setEnabled(false);
      m_Controls->m_ToFVisualisationSettingsWidget->setEnabled(false);
      m_Controls->m_ToFCompositeFilterWidget->setEnabled(false);
      m_Controls->tofMeasurementWidget->setEnabled(false);
      m_Controls->SurfacePropertiesBox->setEnabled(false);
    }
  }
}

void QmitkToFUtilView::ActivatedZombieView(berry::IWorkbenchPartReference::Pointer /*zombieView*/)
{
  ResetGUIToDefault();
}

void QmitkToFUtilView::Deactivated()
{

}

void QmitkToFUtilView::Visible()
{
}

void QmitkToFUtilView::Hidden()
{

    ResetGUIToDefault();
}

void QmitkToFUtilView::OnToFCameraConnected()
{
  this->m_SurfaceDisplayCount = 0;
  this->m_2DDisplayCount = 0;
  this->m_ToFImageGrabber = m_Controls->m_ToFConnectionWidget->GetToFImageGrabber();

  this->m_ToFImageRecorder->SetCameraDevice(this->m_ToFImageGrabber->GetCameraDevice());
  m_Controls->m_ToFRecorderWidget->SetParameter(this->m_ToFImageGrabber, this->m_ToFImageRecorder);
  m_Controls->m_ToFRecorderWidget->setEnabled(true);
  m_Controls->m_ToFRecorderWidget->ResetGUIToInitial();

  // initialize measurement widget
  m_Controls->tofMeasurementWidget->InitializeWidget(this->GetRenderWindowPart()->GetRenderWindows(),this->GetDataStorage());

  //TODO
  this->m_RealTimeClock = mitk::RealTimeClock::New();
  this->m_2DTimeBefore = this->m_RealTimeClock->GetCurrentStamp();

  try
  {
    this->m_VideoSource = mitk::OpenCVVideoSource::New();

    this->m_VideoSource->SetVideoCameraInput(0, false);
    this->m_VideoSource->StartCapturing();
    if(!this->m_VideoSource->IsCapturingEnabled())
    {
      MITK_INFO << "unable to initialize video grabbing/playback";
      this->m_VideoEnabled = false;
      m_Controls->m_VideoTextureCheckBox->setEnabled(false);
    }
    else
    {
      this->m_VideoEnabled = true;
      m_Controls->m_VideoTextureCheckBox->setEnabled(true);
    }

    if (this->m_VideoEnabled)
    {

      this->m_VideoSource->FetchFrame();
      this->m_VideoCaptureHeight = this->m_VideoSource->GetImageHeight();
      this->m_VideoCaptureWidth = this->m_VideoSource->GetImageWidth();
      this->m_VideoTexture = this->m_VideoSource->GetVideoTexture();

      this->m_ToFDistanceImageToSurfaceFilter->SetTextureImageWidth(this->m_VideoCaptureWidth);
      this->m_ToFDistanceImageToSurfaceFilter->SetTextureImageHeight(this->m_VideoCaptureHeight);


      this->m_ToFSurfaceVtkMapper3D->SetTextureWidth(this->m_VideoCaptureWidth);
      this->m_ToFSurfaceVtkMapper3D->SetTextureHeight(this->m_VideoCaptureHeight);
    }
  }
  catch (std::logic_error& e)
  {
    QMessageBox::warning(NULL, "Warning", QString(e.what()));
    MITK_ERROR << e.what();
    return;
  }

  this->RequestRenderWindowUpdate();
}

void QmitkToFUtilView::ResetGUIToDefault()
{
  if(this->GetRenderWindowPart())
  {
    mitk::ILinkedRenderWindowPart* linkedRenderWindowPart = dynamic_cast<mitk::ILinkedRenderWindowPart*>(this->GetRenderWindowPart());
    if(linkedRenderWindowPart == 0)
    {
      MITK_ERROR << "No linked StdMultiWidget avaiable!!!";
    }
    else
    {
      linkedRenderWindowPart->EnableSlicingPlanes(true);
    }
    GetRenderWindowPart()->GetRenderWindow("axial")->GetSliceNavigationController()->SetDefaultViewDirection(mitk::SliceNavigationController::Axial);
    GetRenderWindowPart()->GetRenderWindow("axial")->GetSliceNavigationController()->SliceLockedOff();
    GetRenderWindowPart()->GetRenderWindow("sagittal")->GetSliceNavigationController()->SetDefaultViewDirection(mitk::SliceNavigationController::Sagittal);
    GetRenderWindowPart()->GetRenderWindow("sagittal")->GetSliceNavigationController()->SliceLockedOff();
    GetRenderWindowPart()->GetRenderWindow("coronal")->GetSliceNavigationController()->SetDefaultViewDirection(mitk::SliceNavigationController::Frontal);
    GetRenderWindowPart()->GetRenderWindow("coronal")->GetSliceNavigationController()->SliceLockedOff();

    this->UseToFVisibilitySettings(false);

    //global reinit
    this->GetRenderWindowPart()->GetRenderingManager()->InitializeViews(/*this->GetDataStorage()->ComputeBoundingGeometry3D(this->GetDataStorage()->GetAll())*/);
    this->RequestRenderWindowUpdate();
  }
}

void QmitkToFUtilView::OnToFCameraDisconnected()
{
  m_Controls->m_ToFRecorderWidget->OnStop();
  m_Controls->m_ToFRecorderWidget->setEnabled(false);
  m_Controls->m_ToFVisualisationSettingsWidget->setEnabled(false);
  m_Controls->tofMeasurementWidget->setEnabled(false);
  m_Controls->SurfacePropertiesBox->setEnabled(false);
  //clean up measurement widget
  m_Controls->tofMeasurementWidget->CleanUpWidget();

  if(this->m_VideoSource)
  {
    this->m_VideoSource->StopCapturing();
    this->m_VideoSource = NULL;
  }
}

void QmitkToFUtilView::OnToFCameraStarted()
{
  if (m_ToFImageGrabber.IsNotNull())
  {
    // initial update of image grabber
    this->m_ToFImageGrabber->Update();

    this->m_ToFCompositeFilter->SetInput(0,this->m_ToFImageGrabber->GetOutput(0));
    this->m_ToFCompositeFilter->SetInput(1,this->m_ToFImageGrabber->GetOutput(1));
    this->m_ToFCompositeFilter->SetInput(2,this->m_ToFImageGrabber->GetOutput(2));

    // initial update of composite filter
    this->m_ToFCompositeFilter->Update();
    this->m_MitkDistanceImage = m_ToFCompositeFilter->GetOutput(0);
    this->m_DistanceImageNode = ReplaceNodeData("Distance image",m_MitkDistanceImage);
    this->m_MitkAmplitudeImage = m_ToFCompositeFilter->GetOutput(1);
    this->m_AmplitudeImageNode = ReplaceNodeData("Amplitude image",m_MitkAmplitudeImage);
    this->m_MitkIntensityImage = m_ToFCompositeFilter->GetOutput(2);
    this->m_IntensityImageNode = ReplaceNodeData("Intensity image",m_MitkIntensityImage);

    std::string rgbFileName;
    m_ToFImageGrabber->GetCameraDevice()->GetStringProperty("RGBImageFileName",rgbFileName);
    if ((m_SelectedCamera=="Microsoft Kinect")||(rgbFileName!=""))
    {
      this->m_RGBImageNode = ReplaceNodeData("RGB image",this->m_ToFImageGrabber->GetOutput(3));
    }
    else
    {
      this->m_RGBImageNode = NULL;
    }

    this->m_ToFDistanceImageToSurfaceFilter->SetInput(0,m_MitkDistanceImage);
    this->m_ToFDistanceImageToSurfaceFilter->SetInput(1,m_MitkAmplitudeImage);
    this->m_ToFDistanceImageToSurfaceFilter->SetInput(2,m_MitkIntensityImage);
    this->m_Surface = this->m_ToFDistanceImageToSurfaceFilter->GetOutput(0);
    this->m_SurfaceNode = ReplaceNodeData("Surface",m_Surface);

    this->UseToFVisibilitySettings(true);

    m_Controls->m_ToFCompositeFilterWidget->UpdateFilterParameter();
    // initialize visualization widget
    m_Controls->m_ToFVisualisationSettingsWidget->Initialize(this->m_DistanceImageNode, this->m_AmplitudeImageNode, this->m_IntensityImageNode);
    // set distance image to measurement widget
    m_Controls->tofMeasurementWidget->SetDistanceImage(m_MitkDistanceImage);

    this->m_Frametimer->start(0);

    m_Controls->m_ToFVisualisationSettingsWidget->setEnabled(true);
    m_Controls->m_ToFCompositeFilterWidget->setEnabled(true);
    m_Controls->tofMeasurementWidget->setEnabled(true);
    m_Controls->SurfacePropertiesBox->setEnabled(true);

    if (m_Controls->m_TextureCheckBox->isChecked())
    {
      OnTextureCheckBoxChecked(true);
    }
    if (m_Controls->m_VideoTextureCheckBox->isChecked())
    {
      OnVideoTextureCheckBoxChecked(true);
    }
  }
  m_Controls->m_TextureCheckBox->setEnabled(true);
}

void QmitkToFUtilView::OnToFCameraStopped()
{
  m_Controls->m_ToFVisualisationSettingsWidget->setEnabled(false);
  m_Controls->m_ToFCompositeFilterWidget->setEnabled(false);
  m_Controls->tofMeasurementWidget->setEnabled(false);
  m_Controls->SurfacePropertiesBox->setEnabled(false);

  this->m_Frametimer->stop();
}

void QmitkToFUtilView::OnToFCameraSelected(const QString selected)
{
  m_SelectedCamera = selected;
  if ((selected=="PMD CamBoard")||(selected=="PMD O3D"))
  {
    MITK_INFO<<"Surface representation currently not available for CamBoard and O3. Intrinsic parameters missing.";
    this->m_Controls->m_SurfaceCheckBox->setEnabled(false);
    this->m_Controls->m_TextureCheckBox->setEnabled(false);
    this->m_Controls->m_VideoTextureCheckBox->setEnabled(false);
    this->m_Controls->m_SurfaceCheckBox->setChecked(false);
    this->m_Controls->m_TextureCheckBox->setChecked(false);
    this->m_Controls->m_VideoTextureCheckBox->setChecked(false);
  }
  else
  {
    this->m_Controls->m_SurfaceCheckBox->setEnabled(true);
    this->m_Controls->m_TextureCheckBox->setEnabled(true); // TODO enable when bug 8106 is solved
    this->m_Controls->m_VideoTextureCheckBox->setEnabled(true);
  }
}

void QmitkToFUtilView::OnUpdateCamera()
{
  if (m_Controls->m_VideoTextureCheckBox->isChecked() && this->m_VideoEnabled && this->m_VideoSource)
  {
    this->m_VideoTexture = this->m_VideoSource->GetVideoTexture();
    ProcessVideoTransform();
  }

  if (m_Controls->m_SurfaceCheckBox->isChecked())
  {
    // update surface
    m_ToFDistanceImageToSurfaceFilter->SetTextureIndex(m_Controls->m_ToFVisualisationSettingsWidget->GetSelectedImageIndex());
    this->m_Surface->Update();

    vtkColorTransferFunction* colorTransferFunction = m_Controls->m_ToFVisualisationSettingsWidget->GetSelectedColorTransferFunction();

    this->m_ToFSurfaceVtkMapper3D->SetVtkScalarsToColors(colorTransferFunction);

    if (this->m_SurfaceDisplayCount<2)
    {
      this->m_SurfaceNode->SetData(this->m_Surface);
      this->m_SurfaceNode->SetMapper(mitk::BaseRenderer::Standard3D, m_ToFSurfaceVtkMapper3D);

      this->GetRenderWindowPart()->GetRenderingManager()->InitializeViews(
        this->m_Surface->GetTimeSlicedGeometry(), mitk::RenderingManager::REQUEST_UPDATE_3DWINDOWS, true);

      mitk::Point3D surfaceCenter= this->m_Surface->GetGeometry()->GetCenter();
      vtkCamera* camera3d = GetRenderWindowPart()->GetRenderWindow("3d")->GetRenderer()->GetVtkRenderer()->GetActiveCamera();
      camera3d->SetPosition(0,0,-50);
      camera3d->SetViewUp(0,-1,0);
      camera3d->SetFocalPoint(0,0,surfaceCenter[2]);
      camera3d->SetViewAngle(40);
      camera3d->SetClippingRange(1, 10000);
    }
    this->m_SurfaceDisplayCount++;

  }
  else
  {
    // update pipeline
    this->m_MitkDistanceImage->Update();
  }

  this->RequestRenderWindowUpdate();

  this->m_2DDisplayCount++;
  if ((this->m_2DDisplayCount % this->m_StepsForFramerate) == 0)
  {
    this->m_2DTimeAfter = this->m_RealTimeClock->GetCurrentStamp() - this->m_2DTimeBefore;
    MITK_INFO << " 2D-Display-framerate (fps): " << this->m_StepsForFramerate / (this->m_2DTimeAfter/1000);
    this->m_2DTimeBefore = this->m_RealTimeClock->GetCurrentStamp();
  }
}

void QmitkToFUtilView::ProcessVideoTransform()
{
  IplImage *src, *dst;
  src = cvCreateImageHeader(cvSize(this->m_VideoCaptureWidth, this->m_VideoCaptureHeight), IPL_DEPTH_8U, 3);
  src->imageData = (char*)this->m_VideoTexture;

  CvPoint2D32f srcTri[3], dstTri[3];
  CvMat* rot_mat = cvCreateMat(2,3,CV_32FC1);
  CvMat* warp_mat = cvCreateMat(2,3,CV_32FC1);
  dst = cvCloneImage(src);
  dst->origin = src->origin;
  cvZero( dst );

  int xOffset = 0;//m_Controls->m_XOffsetSpinBox->value();
  int yOffset = 0;//m_Controls->m_YOffsetSpinBox->value();
  int zoom = 0;//m_Controls->m_ZoomSpinBox->value();

  // Compute warp matrix
  srcTri[0].x = 0 + zoom;
  srcTri[0].y = 0 + zoom;
  srcTri[1].x = src->width - 1 - zoom;
  srcTri[1].y = 0 + zoom;
  srcTri[2].x = 0 + zoom;
  srcTri[2].y = src->height - 1 - zoom;

  dstTri[0].x = 0;
  dstTri[0].y = 0;
  dstTri[1].x = src->width - 1;
  dstTri[1].y = 0;
  dstTri[2].x = 0;
  dstTri[2].y = src->height - 1;

  cvGetAffineTransform( srcTri, dstTri, warp_mat );
  cvWarpAffine( src, dst, warp_mat );
  cvCopy ( dst, src );


  // Compute warp matrix
  srcTri[0].x = 0;
  srcTri[0].y = 0;
  srcTri[1].x = src->width - 1;
  srcTri[1].y = 0;
  srcTri[2].x = 0;
  srcTri[2].y = src->height - 1;

  dstTri[0].x = srcTri[0].x + xOffset;
  dstTri[0].y = srcTri[0].y + yOffset;
  dstTri[1].x = srcTri[1].x + xOffset;
  dstTri[1].y = srcTri[1].y + yOffset;
  dstTri[2].x = srcTri[2].x + xOffset;
  dstTri[2].y = srcTri[2].y + yOffset;

  cvGetAffineTransform( srcTri, dstTri, warp_mat );
  cvWarpAffine( src, dst, warp_mat );
  cvCopy ( dst, src );

  src->imageData = NULL;
  cvReleaseImage( &src );
  cvReleaseImage( &dst );
  cvReleaseMat( &rot_mat );
  cvReleaseMat( &warp_mat );

}

void QmitkToFUtilView::OnTextureCheckBoxChecked(bool checked)
{
  if(m_SurfaceNode.IsNotNull())
  {
    if (checked)
    {
      this->m_SurfaceNode->SetBoolProperty("scalar visibility", true);
    }
    else
    {
      this->m_SurfaceNode->SetBoolProperty("scalar visibility", false);
    }
  }
}

void QmitkToFUtilView::OnVideoTextureCheckBoxChecked(bool checked)
{
  if (checked)
  {
    if (this->m_VideoEnabled)
    {
      this->m_ToFSurfaceVtkMapper3D->SetTexture(this->m_VideoTexture);
    }
    else
    {
      this->m_ToFSurfaceVtkMapper3D->SetTexture(NULL);
    }
  }
}

void QmitkToFUtilView::OnChangeCoronalWindowOutput(int index)
{
  this->OnToFCameraStopped();
  if(index == 0)
  {
    if(this->m_IntensityImageNode.IsNotNull())
      this->m_IntensityImageNode->SetVisibility(false);
    if(this->m_RGBImageNode.IsNotNull())
      this->m_RGBImageNode->SetVisibility(true);
  }
  else if(index == 1)
  {
    if(this->m_IntensityImageNode.IsNotNull())
      this->m_IntensityImageNode->SetVisibility(true);
    if(this->m_RGBImageNode.IsNotNull())
      this->m_RGBImageNode->SetVisibility(false);
  }
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
  this->OnToFCameraStarted();
}

mitk::DataNode::Pointer QmitkToFUtilView::ReplaceNodeData( std::string nodeName, mitk::BaseData* data )
{

  mitk::DataNode::Pointer node = this->GetDataStorage()->GetNamedNode(nodeName);
  if (node.IsNull())
  {
    node = mitk::DataNode::New();
    node->SetData(data);
    node->SetName(nodeName);
    node->SetBoolProperty("binary",false);
    this->GetDataStorage()->Add(node);
  }
  else
  {
    node->SetData(data);
  }
  return node;
}

void QmitkToFUtilView::UseToFVisibilitySettings(bool useToF)
{
  // set node properties
  if (m_DistanceImageNode.IsNotNull())
  {
    this->m_DistanceImageNode->SetProperty( "visible" , mitk::BoolProperty::New( true ));
    this->m_DistanceImageNode->SetVisibility( !useToF, mitk::BaseRenderer::GetInstance(GetRenderWindowPart()->GetRenderWindow("sagittal")->GetRenderWindow() ) );
    this->m_DistanceImageNode->SetVisibility( !useToF, mitk::BaseRenderer::GetInstance(GetRenderWindowPart()->GetRenderWindow("coronal")->GetRenderWindow() ) );
    this->m_DistanceImageNode->SetVisibility( !useToF, mitk::BaseRenderer::GetInstance(GetRenderWindowPart()->GetRenderWindow("3d")->GetRenderWindow() ) );
    this->m_DistanceImageNode->SetBoolProperty("use color",!useToF);
    this->m_DistanceImageNode->GetPropertyList()->DeleteProperty("LookupTable");
  }
  if (m_AmplitudeImageNode.IsNotNull())
  {
    this->m_AmplitudeImageNode->SetProperty( "visible" , mitk::BoolProperty::New( true ));
    this->m_AmplitudeImageNode->SetVisibility( !useToF, mitk::BaseRenderer::GetInstance(GetRenderWindowPart()->GetRenderWindow("axial")->GetRenderWindow() ) );
    this->m_AmplitudeImageNode->SetVisibility( !useToF, mitk::BaseRenderer::GetInstance(GetRenderWindowPart()->GetRenderWindow("coronal")->GetRenderWindow() ) );
    this->m_AmplitudeImageNode->SetVisibility( !useToF, mitk::BaseRenderer::GetInstance(GetRenderWindowPart()->GetRenderWindow("3d")->GetRenderWindow() ) );
    this->m_AmplitudeImageNode->SetBoolProperty("use color",!useToF);
    this->m_AmplitudeImageNode->GetPropertyList()->DeleteProperty("LookupTable");
  }
  if (m_IntensityImageNode.IsNotNull())
  {
    this->m_IntensityImageNode->SetProperty( "visible" , mitk::BoolProperty::New( true ));
    this->m_IntensityImageNode->SetVisibility( !useToF, mitk::BaseRenderer::GetInstance(GetRenderWindowPart()->GetRenderWindow("axial")->GetRenderWindow() ) );
    this->m_IntensityImageNode->SetVisibility( !useToF, mitk::BaseRenderer::GetInstance(GetRenderWindowPart()->GetRenderWindow("sagittal")->GetRenderWindow() ) );
    this->m_IntensityImageNode->SetVisibility( !useToF, mitk::BaseRenderer::GetInstance(GetRenderWindowPart()->GetRenderWindow("3d")->GetRenderWindow() ) );
    this->m_IntensityImageNode->SetBoolProperty("use color",!useToF);
    this->m_IntensityImageNode->GetPropertyList()->DeleteProperty("LookupTable");
  }
  if ((m_RGBImageNode.IsNotNull()))
  {
    this->m_RGBImageNode->SetProperty( "visible" , mitk::BoolProperty::New( true ));
    this->m_RGBImageNode->SetVisibility( !useToF, mitk::BaseRenderer::GetInstance(GetRenderWindowPart()->GetRenderWindow("axial")->GetRenderWindow() ) );
    this->m_RGBImageNode->SetVisibility( !useToF, mitk::BaseRenderer::GetInstance(GetRenderWindowPart()->GetRenderWindow("sagittal")->GetRenderWindow() ) );
    this->m_RGBImageNode->SetVisibility( !useToF, mitk::BaseRenderer::GetInstance(GetRenderWindowPart()->GetRenderWindow("3d")->GetRenderWindow() ) );
  }
  // initialize images
  if (m_MitkDistanceImage.IsNotNull())
  {
    this->GetRenderWindowPart()->GetRenderingManager()->InitializeViews(
      this->m_MitkDistanceImage->GetTimeSlicedGeometry(), mitk::RenderingManager::REQUEST_UPDATE_2DWINDOWS, true);
  }
  if(this->m_SurfaceNode.IsNotNull())
  {
    QHash<QString, QmitkRenderWindow*> renderWindowHashMap = this->GetRenderWindowPart()->GetRenderWindows();
    QHashIterator<QString, QmitkRenderWindow*> i(renderWindowHashMap);
    while (i.hasNext()){
      i.next();
      this->m_SurfaceNode->SetVisibility( false, mitk::BaseRenderer::GetInstance(i.value()->GetRenderWindow()) );
    }
    this->m_SurfaceNode->SetVisibility( true, mitk::BaseRenderer::GetInstance(GetRenderWindowPart()->GetRenderWindow("3d")->GetRenderWindow() ) );
  }
  //disable/enable gradient background
  this->GetRenderWindowPart()->EnableDecorations(!useToF, QStringList(QString("background")));
}
