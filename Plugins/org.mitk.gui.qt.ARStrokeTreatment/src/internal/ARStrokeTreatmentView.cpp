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
#include <berryISelectionService.h>
#include <berryIWorkbenchWindow.h>

// Qmitk
#include "ARStrokeTreatmentView.h"

// Qt
#include <QMessageBox>
#include <QmitkDataStorageComboBox.h>

// mitk image
#include <QTimer>
#include <mitkIRenderingManager.h>
#include <mitkImage.h>
#include <mitkImageGenerator.h>
#include <mitkOpenCVToMitkImageFilter.h>

const std::string ARStrokeTreatmentView::VIEW_ID = "org.mitk.views.arstroketreatment";

void ARStrokeTreatmentView::SetFocus()
{
  // m_Controls.buttonPerformImageProcessing->setFocus();
}

void ARStrokeTreatmentView::CreateQtPartControl(QWidget *parent)
{
  // create GUI widgets from the Qt Designer's .ui file
  m_Controls.setupUi(parent);
  // connect(m_Controls.buttonPerformImageProcessing, &QPushButton::clicked, this,
  // &ARStrokeTreatmentView::DoImageProcessing);
  CreateConnections();
  // create timer for tracking and video grabbing
  m_UpdateTimer = new QTimer(this);
  m_UpdateTimer->setInterval(100);
}

void ARStrokeTreatmentView::OnSelectionChanged(berry::IWorkbenchPart::Pointer /*source*/,
                                               const QList<mitk::DataNode::Pointer> &nodes)
{
  // iterate all selected objects, adjust warning visibility
  foreach (mitk::DataNode::Pointer node, nodes)
  {
    if (node.IsNotNull() && dynamic_cast<mitk::Image *>(node->GetData()))
    {
      // m_Controls.labelWarning->setVisible(false);
      return;
    }
  }

  // m_Controls.labelWarning->setVisible(true);
}

void ARStrokeTreatmentView::CreateConnections()
{
  // connect(m_Controls.m_TrackingDeviceSelectionWidget,
  //        SIGNAL(NavigationDataSourceSelected(mitk::NavigationDataSource::Pointer)),
  //        this,
  //        SLOT(OnSetupNavigation()));
  connect(m_Controls.m_TrackerGrabbingPushButton, SIGNAL(clicked()), this, SLOT(OnTrackingGrabberPushed()));
  connect(m_Controls.m_VideoGrabbingPushButton, SIGNAL(clicked()), this, SLOT(OnVideoGrabberPushed()));
  return;
}

void ARStrokeTreatmentView::OnTrackingGrabberPushed()
{
  if (!m_TrackingActive)
  {
    m_TrackingActive = true;
    m_Controls.m_TrackerGrabbingPushButton->setText("Stop Tracking");
    connect(m_UpdateTimer, SIGNAL(timeout()), this, SLOT(UpdateTrackingData()));
    if (!m_UpdateTimer->isActive())
    {
      m_UpdateTimer->start();
    }
  }
  if (m_TrackingActive)
  {
    m_TrackingActive = false;
    m_Controls.m_TrackerGrabbingPushButton->setText("Start Tracking");
    disconnect(m_UpdateTimer, SIGNAL(timeout()), this, SLOT(UpdateTrackingData()));
    // m_TrackingSource = NULL; // check if necessary
    if (!m_TrackingActive && !m_VideoGrabbingActive) // Check if timer can be deactiated
    {
      m_UpdateTimer = NULL;
    }
  }
  return;
}

void ARStrokeTreatmentView::UpdateTrackingData()
{
  m_TrackingData = m_Controls.m_TrackingDeviceSelectionWidget->GetSelectedNavigationDataSource()->GetOutput(
    m_Controls.m_TrackingDeviceSelectionWidget->GetSelectedToolID());
  mitk::NavigationData::Pointer navData =
    m_Controls.m_TrackingDeviceSelectionWidget->GetSelectedNavigationDataSource()->GetOutput(0);
  MITK_INFO << navData->GetPosition();
  return;
}

void ARStrokeTreatmentView::OnVideoGrabberPushed()
{
  if (!m_VideoGrabbingActive) // If not grabbing video data
  {
    m_VideoGrabbingActive = true;
    m_Controls.m_VideoGrabbingPushButton->setText("Stop Video");
    connect(m_UpdateTimer, SIGNAL(timeout()), this, SLOT(UpdateImageData()));
    if (!m_UpdateTimer->isActive())
    {
      m_UpdateTimer->start();
    }
  }
  if (m_VideoGrabbingActive)
  {
    m_VideoGrabbingActive = false;
    m_Controls.m_VideoGrabbingPushButton->setText("Start Video");
    disconnect(m_UpdateTimer, SIGNAL(timeout()), this, SLOT(UpdateImageData()));
    // m_Controls.m_StartGrabbing->setText("Start Video Grabbing");
    if (!m_TrackingActive && !m_VideoGrabbingActive) // Check if timer can be deactiated
    {
      m_UpdateTimer = NULL;
    }
  }
   else if (m_Controls.m_MITKImage->isChecked())
  {
     /* mitk::DataNode::Pointer imageNode = mitk::DataNode::New(); // CHECK IF NEEDED
    imageNode->SetName("Open CV Example Image Stream");
    imageNode->SetData(mitk::Image::New());
    m_ConversionFilter = mitk::OpenCVToMitkImageFilter::New();
    this->GetDataStorage()->Add(imageNode);
    OnUpdateImage();
    // Initialize view on Image
    mitk::IRenderWindowPart *renderWindow = this->GetRenderWindowPart();
    if (renderWindow != NULL)
      renderWindow->GetRenderingManager()->InitializeViews(
        imageNode->GetData()->GetGeometry(), mitk::RenderingManager::REQUEST_UPDATE_ALL, true);

    m_UpdateTimerTracking->setInterval(20);
    m_UpdateTimerTracking->start();*/
  }
}

void ARStrokeTreatmentView::UpdateImageData()
{
  m_VideoCapture = new cv::VideoCapture(); // open the default camera
  if (!m_VideoCapture->isOpened())
  {
    return;
  }         // check if we succeeded
  if (true) // m_Controls.m_SeparateWindow->isChecked()
  {
    mitk::DataNode::Pointer imageNode = mitk::DataNode::New();
    std::stringstream nodeName;
    nodeName << "Live Image Stream";
    imageNode->SetName(nodeName.str());
    // create a dummy image (gray values 0..255) for correct initialization of level window, etc.
    mitk::Image::Pointer dummyImage = mitk::ImageGenerator::GenerateRandomImage<float>(100, 100, 1, 1, 1, 1, 1, 255, 0);
    imageNode->SetData(dummyImage);
    while (m_VideoGrabbingActive)
      this->GetDataStorage()->Add(imageNode);
    cv::Mat frame;
    {
      *m_VideoCapture >> frame; // get a new frame from camera
      m_ConversionFilter->SetOpenCVMat(frame);
      m_ConversionFilter->Update();
      m_ConversionFilter->GetOutput();
    }
  }
}

void ARStrokeTreatmentView::DoImageProcessing()
{
  QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
  if (nodes.empty())
    return;

  mitk::DataNode *node = nodes.front();

  if (!node)
  {
    // Nothing selected. Inform the user and return
    QMessageBox::information(nullptr, "Template", "Please load and select an image before starting image processing.");
    return;
  }

  // here we have a valid mitk::DataNode

  // a node itself is not very useful, we need its data item (the image)
  mitk::BaseData *data = node->GetData();
  if (data)
  {
    // test if this data item is an image or not (could also be a surface or something totally different)
    mitk::Image *image = dynamic_cast<mitk::Image *>(data);
    if (image)
    {
      std::stringstream message;
      std::string name;
      message << "Performing image processing for image ";
      if (node->GetName(name))
      {
        // a property called "name" was found for this DataNode
        message << "'" << name << "'";
      }
      message << ".";
      MITK_INFO << message.str();

      // actually do something here...
    }
  }
}

void ARStrokeTreatmentView::TestText()
{
  MITK_INFO << "TestText succesfully printed! Yay!";
  return;
}
