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
#include "ARStrokeTreatmentRegistration.h"
#include "ARStrokeTreatmentView.h"

// Qt
#include <QMessageBox>
#include <QSettings>
#include <QThread>
#include <QTimer>
#include <QmitkDataStorageComboBox.h>
#include <qfiledialog.h>

// MITK
#include "internal/org_mitk_gui_qt_ARStrokeTreatment_Activator.h"
#include <mitkCone.h>
#include <mitkIOUtil.h>
#include <mitkIRenderingManager.h>
#include <mitkImage.h>
#include <mitkImageGenerator.h>
#include <mitkLog.h>
#include <mitkNavigationToolStorageDeserializer.h>
#include <mitkNavigationToolStorageSerializer.h>
#include <mitkNeedleProjectionFilter.h>
#include <mitkNodePredicateDataType.h>
#include <mitkNodePredicateNot.h>
#include <mitkNodePredicateOr.h>
#include <mitkNodePredicateProperty.h>
#include <mitkOpenCVToMitkImageFilter.h>
#include <mitkTrackingDeviceTypeCollection.h>
#include <mitkTrackingVolumeGenerator.h>
#include <mitkUnspecifiedTrackingTypeInformation.h>
// for exceptions
#include <mitkIGTException.h>
#include <mitkIGTIOException.h>

const std::string ARStrokeTreatmentView::VIEW_ID = "org.mitk.views.arstroketreatment";

ARStrokeTreatmentView::ARStrokeTreatmentView()
  : QmitkAbstractView(), m_Controls(nullptr), m_DeviceTypeCollection(nullptr), m_ToolProjectionNode(nullptr)
{
  m_TrackingLoggingTimer = new QTimer(this);
  m_TrackingRenderTimer = new QTimer(this);
  m_TimeoutTimer = new QTimer(this);
  m_tracking = false;
  m_connected = false;
  m_logging = false;
  m_ShowHideToolAxis = false;
  m_loggedFrames = 0;
  m_SimpleModeEnabled = false;
  m_NeedleProjectionFilter = mitk::NeedleProjectionFilter::New();

  // create filename for autosaving of tool storage
  QString loggingPathWithoutFilename = QString(mitk::LoggingBackend::GetLogFile().c_str());
  if (!loggingPathWithoutFilename.isEmpty()) // if there already is a path for the MITK logging file use this one
  {
    // extract path from path+filename (if someone knows a better way to do this feel free to change it)
    int lengthOfFilename = QFileInfo(QString::fromStdString(mitk::LoggingBackend::GetLogFile())).fileName().size();
    loggingPathWithoutFilename.resize(loggingPathWithoutFilename.size() - lengthOfFilename);
    m_AutoSaveFilename = loggingPathWithoutFilename + "TrackingToolboxAutoSave.IGTToolStorage";
  }
  else // if not: use a temporary path from IOUtil
  {
    m_AutoSaveFilename = QString(mitk::IOUtil::GetTempPath().c_str()) + "TrackingToolboxAutoSave.IGTToolStorage";
  }
  MITK_INFO("IGT Tracking Toolbox") << "Filename for auto saving of IGT ToolStorages: "
                                    << m_AutoSaveFilename.toStdString();

  //! [Thread 1]
  // initialize worker thread
  m_WorkerThread = new QThread();
  m_Worker = new ARStrokeTreatmentTracking();
  //! [Thread 1]

  ctkPluginContext *pluginContext = mitk::PluginActivator::GetContext();
  if (pluginContext)
  {
    QString interfaceName = QString::fromStdString(us_service_interface_iid<mitk::TrackingDeviceTypeCollection>());
    QList<ctkServiceReference> serviceReference = pluginContext->getServiceReferences(interfaceName);

    if (serviceReference.size() > 0)
    {
      m_DeviceTypeServiceReference = serviceReference.at(0);
      const ctkServiceReference &r = serviceReference.at(0);
      m_DeviceTypeCollection = pluginContext->getService<mitk::TrackingDeviceTypeCollection>(r);
    }
    else
    {
      MITK_INFO << "No Tracking Device Collection!";
    }
  }
}
ARStrokeTreatmentView::~ARStrokeTreatmentView()
{
  this->StoreUISettings();
  m_TrackingLoggingTimer->stop();
  m_TrackingRenderTimer->stop();
  m_TimeoutTimer->stop();
  delete m_TrackingLoggingTimer;
  delete m_TrackingRenderTimer;
  delete m_TimeoutTimer;
  try
  {
    //! [Thread 2]
    // wait for thread to finish
    m_WorkerThread->terminate();
    m_WorkerThread->wait();
    // clean up worker thread
    if (m_WorkerThread)
    {
      delete m_WorkerThread;
    }
    if (m_Worker)
    {
      delete m_Worker;
    }
    //! [Thread 2]

    // remove the tracking volume
    this->GetDataStorage()->Remove(m_TrackingVolumeNode);
    // unregister microservices
    if (m_toolStorage)
    {
      m_toolStorage->UnRegisterMicroservice();
    }

    if (m_IGTLMessageProvider.IsNotNull())
    {
      m_IGTLMessageProvider->UnRegisterMicroservice();
    }
  }
  catch (std::exception &e)
  {
    MITK_WARN << "Unexpected exception during clean up of tracking toolbox view: " << e.what();
  }
  catch (...)
  {
    MITK_WARN << "Unexpected unknown error during clean up of tracking toolbox view!";
  }
  // store tool storage and UI settings for persistence
  this->AutoSaveToolStorage();
  this->StoreUISettings();

  m_DeviceTypeCollection = nullptr;
  mitk::PluginActivator::GetContext()->ungetService(m_DeviceTypeServiceReference);
};

void ARStrokeTreatmentView::SetFocus()
{
  // m_Controls.buttonPerformImageProcessing->setFocus();
}

void ARStrokeTreatmentView::CreateQtPartControl(QWidget *parent)
{
  // create GUI widgets from the Qt Designer's .ui file
  m_Controls = new Ui::ARStrokeTreatmentControls;
  m_Controls->setupUi(parent);
  // connect(m_Controls.buttonPerformImageProcessing, &QPushButton::clicked, this,
  // &ARStrokeTreatmentView::DoImageProcessing);
  // create timer for tracking and video grabbing
  m_UpdateTimer = new QTimer(this);
  m_UpdateTimer->setInterval(100);
  m_UpdateTimer->start();

  CreateConnections();

  m_Controls->m_RegistrationWidget->setDataStorage(this->GetDataStorage());
  m_Controls->m_RegistrationWidget->HideStaticRegistrationRadioButton(true);
  m_Controls->m_RegistrationWidget->HideContinousRegistrationRadioButton(true);
  m_Controls->m_RegistrationWidget->HideUseICPRegistrationCheckbox(true);

  m_Controls->m_DataStorageComboBox->SetDataStorage(this->GetDataStorage());
  m_Controls->m_DataStorageComboBox->SetAutoSelectNewItems(false);

  m_Controls->m_VideoPausePushButton->setDisabled(true);
}

void ARStrokeTreatmentView::OnSelectionChanged(berry::IWorkbenchPart::Pointer /*source*/,
                                               const QList<mitk::DataNode::Pointer> &nodes)
{
  m_Controls->m_AutomaticRegistrationWidget->Initialize(this->GetDataStorage());
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
  connect(m_Controls->m_ScalingPushButton, SIGNAL(clicked()), this, SLOT(OnScalingChanged()));
  connect(m_Controls->m_ScalingComboBox, SIGNAL(currentIndexChanged(int)), this, SLOT(OnScalingComboBoxChanged()));
  connect(m_Controls->m_Transform, SIGNAL(clicked()), this, SLOT(OnTransformClicked()));
  connect(m_Controls->m_CalculateTREPushButton, SIGNAL(clicked()), this, SLOT(OnCalculateTREClicked()));
  connect(m_Controls->m_ChangeDisplayPushButton, SIGNAL(clicked()), this, SLOT(OnChangeDisplayStyle()));
  // connect(m_Controls.m_TrackingDeviceSelectionWidget,
  //        SIGNAL(NavigationDataSourceSelected(mitk::NavigationDataSource::Pointer)),
  //        this,
  //        SLOT(OnSetupNavigation()));
  connect(m_Controls->m_VideoGrabbingPushButton, SIGNAL(clicked()), this, SLOT(OnVideoGrabberPushed()));
  connect(m_Controls->m_VideoPausePushButton, SIGNAL(clicked()), this, SLOT(OnVideoPausePushButton()));
  connect(m_UpdateTimer, SIGNAL(timeout()), this, SLOT(UpdateLiveData()));

  // create connections for Registration Widget
  connect(m_Controls->m_ChooseSelectedPointer, SIGNAL(clicked()), this, SLOT(PointerSelectionChanged()));
  connect(m_Controls->m_ChooseSelectedImage, SIGNAL(clicked()), this, SLOT(ImageSelectionChanged()));
  // create connections
  connect(m_Controls->m_LoadTools, SIGNAL(clicked()), this, SLOT(OnLoadTools()));
  connect(m_Controls->m_ConnectDisconnectButton, SIGNAL(clicked()), this, SLOT(OnConnectDisconnect()));
  connect(m_Controls->m_StartStopTrackingButton, SIGNAL(clicked()), this, SLOT(OnStartStopTracking()));
  connect(m_Controls->m_ConnectSimpleMode, SIGNAL(clicked()), this, SLOT(OnConnectDisconnect()));
  connect(m_Controls->m_StartTrackingSimpleMode, SIGNAL(clicked()), this, SLOT(OnStartStopTracking()));
  connect(m_Controls->m_FreezeUnfreezeTrackingButton, SIGNAL(clicked()), this, SLOT(OnFreezeUnfreezeTracking()));
  connect(m_TrackingLoggingTimer, SIGNAL(timeout()), this, SLOT(UpdateLoggingTrackingTimer()));
  connect(m_TrackingRenderTimer, SIGNAL(timeout()), this, SLOT(UpdateRenderTrackingTimer()));
  connect(m_TimeoutTimer, SIGNAL(timeout()), this, SLOT(OnTimeOut()));
  connect(m_Controls->m_ChooseFile, SIGNAL(clicked()), this, SLOT(OnChooseFileClicked()));
  connect(m_Controls->m_StartLogging, SIGNAL(clicked()), this, SLOT(StartLogging()));
  connect(m_Controls->m_StopLogging, SIGNAL(clicked()), this, SLOT(StopLogging()));
  connect(m_Controls->m_VolumeSelectionBox,
          SIGNAL(currentIndexChanged(QString)),
          this,
          SLOT(OnTrackingVolumeChanged(QString)));
  connect(m_Controls->m_ShowTrackingVolume, SIGNAL(clicked()), this, SLOT(OnShowTrackingVolumeChanged()));
  connect(m_Controls->m_AutoDetectTools, SIGNAL(clicked()), this, SLOT(OnAutoDetectTools()));
  connect(m_Controls->m_ResetTools, SIGNAL(clicked()), this, SLOT(OnResetTools()));
  connect(m_Controls->m_AddSingleTool, SIGNAL(clicked()), this, SLOT(OnAddSingleTool()));
  connect(m_Controls->m_NavigationToolCreationWidget,
          SIGNAL(NavigationToolFinished()),
          this,
          SLOT(OnAddSingleToolFinished()));
  connect(m_Controls->m_NavigationToolCreationWidget, SIGNAL(Canceled()), this, SLOT(OnAddSingleToolCanceled()));
  connect(m_Controls->m_CsvFormat, SIGNAL(clicked()), this, SLOT(OnToggleFileExtension()));
  connect(m_Controls->m_XmlFormat, SIGNAL(clicked()), this, SLOT(OnToggleFileExtension()));
  connect(m_Controls->m_UseDifferentUpdateRates, SIGNAL(clicked()), this, SLOT(OnToggleDifferentUpdateRates()));
  connect(m_Controls->m_RenderUpdateRate, SIGNAL(valueChanged(int)), this, SLOT(OnChangeRenderUpdateRate()));
  connect(m_Controls->m_DisableAllTimers, SIGNAL(stateChanged(int)), this, SLOT(EnableDisableTimerButtons(int)));
  connect(m_Controls->m_advancedUI, SIGNAL(clicked()), this, SLOT(OnToggleAdvancedSimpleMode()));
  connect(m_Controls->m_SimpleUI, SIGNAL(clicked()), this, SLOT(OnToggleAdvancedSimpleMode()));
  connect(m_Controls->showHideToolProjectionCheckBox, SIGNAL(clicked()), this, SLOT(OnShowHideToolProjectionClicked()));
  connect(m_Controls->showHideToolAxisCheckBox, SIGNAL(clicked()), this, SLOT(OnShowHideToolAxisClicked()));

  connect(m_Controls->m_toolselector, SIGNAL(currentIndexChanged(int)), this, SLOT(SelectToolProjection(int)));

  // connections for the tracking device configuration widget
  connect(
    m_Controls->m_ConfigurationWidget, SIGNAL(TrackingDeviceSelectionChanged()), this, SLOT(OnTrackingDeviceChanged()));

  //! [Thread 3]
  // connect worker thread
  connect(
    m_Worker, SIGNAL(AutoDetectToolsFinished(bool, QString)), this, SLOT(OnAutoDetectToolsFinished(bool, QString)));
  connect(m_Worker, SIGNAL(ConnectDeviceFinished(bool, QString)), this, SLOT(OnConnectFinished(bool, QString)));
  connect(m_Worker, SIGNAL(StartTrackingFinished(bool, QString)), this, SLOT(OnStartTrackingFinished(bool, QString)));
  connect(m_Worker, SIGNAL(StopTrackingFinished(bool, QString)), this, SLOT(OnStopTrackingFinished(bool, QString)));
  connect(m_Worker, SIGNAL(DisconnectDeviceFinished(bool, QString)), this, SLOT(OnDisconnectFinished(bool, QString)));
  connect(m_WorkerThread, SIGNAL(started()), m_Worker, SLOT(ThreadFunc()));

  connect(
    m_Worker, SIGNAL(ConnectDeviceFinished(bool, QString)), m_Controls->m_ConfigurationWidget, SLOT(OnConnected(bool)));
  connect(m_Worker,
          SIGNAL(DisconnectDeviceFinished(bool, QString)),
          m_Controls->m_ConfigurationWidget,
          SLOT(OnDisconnected(bool)));
  connect(m_Worker,
          SIGNAL(StartTrackingFinished(bool, QString)),
          m_Controls->m_ConfigurationWidget,
          SLOT(OnStartTracking(bool)));
  connect(m_Worker,
          SIGNAL(StopTrackingFinished(bool, QString)),
          m_Controls->m_ConfigurationWidget,
          SLOT(OnStopTracking(bool)));

  // Add Listener, so that we know when the toolStorage changed.
  std::string m_Filter =
    "(" + us::ServiceConstants::OBJECTCLASS() + "=" + "org.mitk.services.NavigationToolStorage" + ")";
  mitk::PluginActivator::GetContext()->connectServiceListener(this, "OnToolStorageChanged", QString(m_Filter.c_str()));

  // move the worker to the thread
  m_Worker->moveToThread(m_WorkerThread);
  //! [Thread 3]

  // initialize widgets
  m_Controls->m_TrackingToolsStatusWidget->SetShowPositions(true);
  m_Controls->m_TrackingToolsStatusWidget->SetTextAlignment(Qt::AlignLeft);
  m_Controls->m_simpleWidget->setVisible(false);

  // initialize tracking volume node
  m_TrackingVolumeNode = mitk::DataNode::New();
  m_TrackingVolumeNode->SetName("TrackingVolume");
  m_TrackingVolumeNode->SetBoolProperty("Backface Culling", true);
  mitk::Color red;
  red.SetRed(1);
  m_TrackingVolumeNode->SetColor(red);

  // initialize buttons
  m_Controls->m_AutoDetectTools->setVisible(false); // only visible if supported by tracking device
  m_Controls->m_StartStopTrackingButton->setEnabled(false);
  m_Controls->m_StartTrackingSimpleMode->setEnabled(false);
  m_Controls->m_FreezeUnfreezeTrackingButton->setEnabled(false);

  // initialize warning labels
  m_Controls->m_RenderWarningLabel->setVisible(false);
  m_Controls->m_TrackingFrozenLabel->setVisible(false);

  // initialize projection buttons
  // first it is disabled when the tool is connected check box for projection and axis is enabled
  m_Controls->showHideToolAxisCheckBox->setEnabled(false);
  m_Controls->showHideToolProjectionCheckBox->setEnabled(false);
  m_Controls->m_toolselector->setEnabled(false);

  // Update List of available models for selected tool.
  std::vector<mitk::TrackingDeviceData> Compatibles;
  if ((m_Controls == nullptr) || // check all these stuff for NULL, latterly this causes crashes from time to time
      (m_Controls->m_ConfigurationWidget == nullptr) ||
      (m_Controls->m_ConfigurationWidget->GetTrackingDevice().IsNull()))
  {
    MITK_ERROR << "Couldn't get current tracking device or an object is nullptr, something went wrong!";
    return;
  }
  else
  {
    Compatibles =
      m_DeviceTypeCollection->GetDeviceDataForLine(m_Controls->m_ConfigurationWidget->GetTrackingDevice()->GetType());
  }
  m_Controls->m_VolumeSelectionBox->clear();
  for (std::size_t i = 0; i < Compatibles.size(); i++)
  {
    m_Controls->m_VolumeSelectionBox->addItem(Compatibles[i].Model.c_str());
  }

  // initialize tool storage
  m_toolStorage = mitk::NavigationToolStorage::New(GetDataStorage());
  m_toolStorage->SetName("TrackingToolbox Default Storage");
  m_toolStorage->RegisterAsMicroservice();

  // set home directory as default path for logfile
  m_Controls->m_LoggingFileName->setText(QDir::toNativeSeparators(QDir::homePath()) + QDir::separator() +
                                         "logfile.csv");

  // tracking device may be changed already by the persistence of the
  // QmitkTrackingDeciveConfigurationWidget
  this->OnTrackingDeviceChanged();

  this->LoadUISettings();

  // add tracking volume node only to data storage
  this->GetDataStorage()->Add(m_TrackingVolumeNode);
  if (!m_Controls->m_ShowTrackingVolume->isChecked())
    m_TrackingVolumeNode->SetOpacity(0.0);
  else
    m_TrackingVolumeNode->SetOpacity(0.25);

  // Update List of available models for selected tool.
  m_Controls->m_VolumeSelectionBox->clear();
  for (std::size_t i = 0; i < Compatibles.size(); i++)
  {
    m_Controls->m_VolumeSelectionBox->addItem(Compatibles[i].Model.c_str());
  }
  return;
}

void ARStrokeTreatmentView::OnTransformClicked()
{
  if (m_Controls->m_TrackingDeviceSelectionWidget->GetSelectedToolID() == -1)
  {
    MITK_INFO << "Could not find a tool, please select a tool in the TrackingDeviceSelectionWidget! ;-)";
    return;
  }
  MITK_INFO << m_Controls->m_TrackingDeviceSelectionWidget->GetSelectedToolID();
  mitk::NavigationData::Pointer transformSensorCSToTracking =
    m_Controls->m_TrackingDeviceSelectionWidget->GetSelectedNavigationDataSource()->GetOutput(0); // get first tool
  mitk::AffineTransform3D::Pointer totalTransformation = mitk::AffineTransform3D::New();
  totalTransformation->SetIdentity();
  mitk::AffineTransform3D::Pointer T1 = m_Controls->m_AutomaticRegistrationWidget->GetTransformMarkerCSToImageCS();
  totalTransformation->Compose(m_Controls->m_AutomaticRegistrationWidget->GetInverseTransform(T1));
  mitk::AffineTransform3D::Pointer T2 = m_Controls->m_AutomaticRegistrationWidget->GetTransformMarkerCSToSensorCS();
  totalTransformation->Compose(T2);
  totalTransformation->Compose(transformSensorCSToTracking->GetAffineTransform3D());
  totalTransformation->Modified();
  // first we have to store the original ct image transform to compose it with the new transform later
  MITK_INFO << "bananarama";
  mitk::AffineTransform3D::Pointer imageTransformNew = mitk::AffineTransform3D::New();
  // create new image transform... setting the composed directly leads to an error
  itk::Matrix<mitk::ScalarType, 3, 3> rotationFloatNew = totalTransformation->GetMatrix();
  itk::Vector<mitk::ScalarType, 3> translationFloatNew = totalTransformation->GetOffset();
  imageTransformNew->SetMatrix(rotationFloatNew);
  imageTransformNew->SetOffset(translationFloatNew);
  m_Controls->m_AutomaticRegistrationWidget;
  m_Controls->m_AutomaticRegistrationWidget;
  m_Controls->m_AutomaticRegistrationWidget->GetSurfaceNode()->GetData()->GetGeometry()->SetIndexToWorldTransform(
    imageTransformNew);
  m_Controls->m_AutomaticRegistrationWidget->GetImageNode()->GetData()->GetGeometry()->SetIndexToWorldTransform(
    imageTransformNew);

  for (int counter = 0; counter < m_Controls->m_AutomaticRegistrationWidget->GetPointSetsDataNodes().size(); ++counter)
  {
    m_Controls->m_AutomaticRegistrationWidget->GetPointSetsDataNodes()
      .at(counter)
      ->GetData()
      ->GetGeometry()
      ->SetIndexToWorldTransform(imageTransformNew);
    m_Controls->m_AutomaticRegistrationWidget->GetPointSetsDataNodes().at(counter)->Modified();
  }
  m_AffineTransform = imageTransformNew;
  m_TransformationSet = true;
  GlobalReinit();
}

void ARStrokeTreatmentView::OnCalculateTREClicked()
{
  mitk::DataNode::Pointer pointSetNode = m_Controls->m_AutomaticRegistrationWidget->GetPointSetNode();
  if (pointSetNode.IsNull())
  {
    MITK_WARN << "Cannot calculate TRE. The pointSetsToChooseComboBox node returned a nullptr.";
    return;
  }

  mitk::PointSet::Pointer pointSet = dynamic_cast<mitk::PointSet *>(pointSetNode->GetData());
  if (pointSet.IsNull())
  {
    m_Controls->m_LabelTRE->setText(QString("Unknown"));
    return;
  }
  mitk::Point3D positionEMTrackedGuideWire =
    m_Controls->m_TrackingDeviceSelectionWidget->GetSelectedNavigationDataSource()->GetOutput(1)->GetPosition();
  double distance = pointSet->GetPoint(0).EuclideanDistanceTo(positionEMTrackedGuideWire);
  m_Controls->m_LabelTRE->setText(QString("%1").arg(distance));
}

void ARStrokeTreatmentView::OnChangeDisplayStyle()
{
  const char *img = "Image";
  mitk::NodePredicateDataType::Pointer npImage = mitk::NodePredicateDataType::New(img);
  itk::SmartPointer<const mitk::DataStorage::SetOfObjects> dnPointer = this->GetDataStorage()->GetAll();
  dnPointer->Size();
  for (size_t i = 0; i < dnPointer->Size(); i++)
  {
    if (npImage->CheckNode(dnPointer->ElementAt(i)))
    {
      dnPointer->ElementAt(i)->SetOpacity(0.5);
    }
  }

  QmitkRenderWindow *renderWindow =
    this->GetRenderWindowPart()->GetQmitkRenderWindow(mitk::BaseRenderer::ViewDirection(2));
  this->GetRenderWindowPart()
    ->GetQmitkRenderWindow(mitk::BaseRenderer::ViewDirection(2))
    ->GetCameraRotationController();

  MITK_INFO << int(renderWindow->GetLayoutIndex());
  renderWindow->SetLayoutIndex(mitk::BaseRenderer::ViewDirection(3));
  MITK_INFO << int(renderWindow->GetLayoutIndex());
  renderWindow->updateBehavior();

  this->GlobalReinit();
  return;
}

void ARStrokeTreatmentView::OnVideoGrabberPushed()
{
  if (!m_VideoGrabbingActive)
  {
    m_Controls->m_VideoGrabbingPushButton->setText("Stop Video");
    // Initialize new video grabber
    m_imageNode = mitk::DataNode::New();
    std::stringstream nodeName;
    nodeName << "Live Image Stream";
    m_imageNode->SetName(nodeName.str());
    // create a dummy image (gray values 0..255) for correct initialization of level window, etc.
    mitk::Image::Pointer dummyImage = mitk::ImageGenerator::GenerateRandomImage<float>(100, 100, 1, 1, 1, 1, 1, 255, 0);
    m_imageNode->SetData(dummyImage);
    this->GetDataStorage()->Add(m_imageNode);

    // select video source
    // m_VideoCapture = new cv::VideoCapture("C:/Tools/7.avi");
    m_VideoCapture = new cv::VideoCapture(0);
    mitk::IRenderWindowPart *renderWindow = this->GetRenderWindowPart();
    renderWindow->GetRenderingManager()->InitializeViews(
      m_imageNode->GetData()->GetGeometry(), mitk::RenderingManager::REQUEST_UPDATE_ALL, false);
    m_Controls->m_VideoPausePushButton->setDisabled(false);
    this->GlobalReinit();
  }
  if (m_VideoGrabbingActive)
  {
    ARStrokeTreatmentView::DisableVideoData();
  }
  m_VideoGrabbingActive = !m_VideoGrabbingActive;
}

void ARStrokeTreatmentView::PointerSelectionChanged()
{
  ARStrokeTreatmentView::InitializeRegistration();
  int toolID = m_Controls->m_TrackingDeviceSelectionWidget->GetSelectedToolID();
  m_TrackingData = m_Controls->m_TrackingDeviceSelectionWidget->GetSelectedNavigationDataSource()->GetOutput(toolID);
  m_Controls->m_RegistrationWidget->setTrackerNavigationData(m_TrackingData);
  m_Controls->m_PointerLabel->setText(
    m_Controls->m_TrackingDeviceSelectionWidget->GetSelectedNavigationTool()->GetToolName().c_str());
}

void ARStrokeTreatmentView::ImageSelectionChanged()
{
  ARStrokeTreatmentView::InitializeRegistration();
  m_Controls->m_ImageLabel->setText(m_Controls->m_DataStorageComboBox->GetSelectedNode()->GetName().c_str());
  m_Controls->m_RegistrationWidget->setImageNode(m_Controls->m_DataStorageComboBox->GetSelectedNode());
}

void ARStrokeTreatmentView::InitializeRegistration()
{
  foreach (QmitkRenderWindow *renderWindow, this->GetRenderWindowPart()->GetQmitkRenderWindows().values())
  {
    this->m_Controls->m_RegistrationWidget->AddSliceNavigationController(renderWindow->GetSliceNavigationController());
  }
}

void ARStrokeTreatmentView::UpdateTrackingData()
{
  m_TrackingData = m_Controls->m_TrackingDeviceSelectionWidget->GetSelectedNavigationDataSource()->GetOutput(
    m_Controls->m_TrackingDeviceSelectionWidget->GetSelectedToolID());
  return;
}

void ARStrokeTreatmentView::OnScalingComboBoxChanged()
{
  if (m_Controls->m_ScalingComboBox->currentIndex() == 2)
  {
    m_Controls->m_ScalingDoubleSpinBox->setEnabled(true);
    return;
  }
  else
  {
    m_Controls->m_ScalingDoubleSpinBox->setEnabled(false);
    return;
  }
  return;
}

void ARStrokeTreatmentView::OnScalingChanged()
{
  double m_ScalingFactor;
  m_ScalingChanged = true;
  if (m_Controls->m_ScalingDoubleSpinBox->isEnabled())
  {
    m_ScalingFactor = m_Controls->m_ScalingDoubleSpinBox->value();
  }
  if (m_Controls->m_ScalingComboBox->currentIndex() == 0)
  {
    m_ScalingFactor = 1;
  }
  if (m_Controls->m_ScalingComboBox->currentIndex() == 1)
  {
    m_ScalingFactor = 0.48;
  }
  m_SetSpacing[0] = m_ScalingFactor; // left-right
  m_SetSpacing[1] = m_ScalingFactor; // up
  m_SetSpacing[2] = m_ScalingFactor; // height
  return;
}

void ARStrokeTreatmentView::UpdateLiveData()
{
  m_Controls->m_TrackingDeviceSelectionWidget->update();
  if (m_TrackingActive)
  {
    UpdateTrackingData();
  }
  if (m_VideoGrabbingActive && m_UpdateVideoData)
  {
    cv::Mat frame;
    if (m_VideoCapture->read(frame))
    {
      m_AffineTransform = m_imageNode->GetData()->GetGeometry()->GetIndexToWorldTransform();
      // m_VideoCapture->read(frame);
      m_ConversionFilter->SetOpenCVMat(frame);
      m_ConversionFilter->Update();
      m_imageNode->SetData(m_ConversionFilter->GetOutput());
      std::stringstream nodeName;
      nodeName << "Live Image Stream";
      m_imageNode->SetName(nodeName.str());
      m_imageNode->GetData()->GetGeometry()->SetIndexToWorldTransform(m_AffineTransform);
      m_imageNode->GetData()->GetGeometry()->SetSpacing(m_SetSpacing);
      m_imageNode->Update();
      this->RequestRenderWindowUpdate(mitk::RenderingManager::REQUEST_UPDATE_ALL);
    }
    else
    {
      MITK_ERROR << "No Image could be read. Video Source not found or finished!";
      ARStrokeTreatmentView::DisableVideoData();
      m_TransformationSet = false;
    }
  }
  return;
}

void ARStrokeTreatmentView::DisableVideoData()
{
  m_VideoGrabbingActive = false;
  m_VideoCapture = new cv::VideoCapture;
  m_Controls->m_VideoGrabbingPushButton->setText("Start Video");
  m_Controls->m_VideoPausePushButton->setDisabled(true);
  this->GetDataStorage()->Remove(m_imageNode);
  return;
}

void ARStrokeTreatmentView::OnVideoPausePushButton()
{
  m_UpdateVideoData = !m_UpdateVideoData;
  if (m_UpdateVideoData)
  {
    m_Controls->m_VideoPausePushButton->setText("Pause Video");
  }
  if (!m_UpdateVideoData)
  {
    m_Controls->m_VideoPausePushButton->setText("Continue Video");
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

void ARStrokeTreatmentView::OnLoadTools()
{
  // read in filename
  QString filename =
    QFileDialog::getOpenFileName(nullptr, tr("Open Tool Storage"), "/", tr("Tool Storage Files (*.IGTToolStorage)"));
  if (filename.isNull())
    return;

  // read tool storage from disk
  std::string errorMessage = "";
  mitk::NavigationToolStorageDeserializer::Pointer myDeserializer =
    mitk::NavigationToolStorageDeserializer::New(GetDataStorage());
  // try-catch block for exceptions
  try
  {
    this->ReplaceCurrentToolStorage(myDeserializer->Deserialize(filename.toStdString()), filename.toStdString());
  }
  catch (mitk::IGTException &)
  {
    std::string errormessage = "Error during loading the tool storage file. Please only load tool storage files "
                               "created with the NavigationToolManager view.";
    QMessageBox::warning(nullptr, "Tool Storage Loading Error", errormessage.c_str());
    return;
  }

  if (m_toolStorage->isEmpty())
  {
    errorMessage = myDeserializer->GetErrorMessage();
    MessageBox(errorMessage);
    return;
  }

  // update label
  UpdateToolStorageLabel(filename);

  // update tool preview
  m_Controls->m_TrackingToolsStatusWidget->RemoveStatusLabels();
  m_Controls->m_TrackingToolsStatusWidget->PreShowTools(m_toolStorage);

  // save filename for persistent storage
  m_ToolStorageFilename = filename;
}

void ARStrokeTreatmentView::OnResetTools()
{
  // remove data nodes of surfaces from data storage to clean up
  for (unsigned int i = 0; i < m_toolStorage->GetToolCount(); i++)
  {
    this->GetDataStorage()->Remove(m_toolStorage->GetTool(i)->GetDataNode());
  }
  this->ReplaceCurrentToolStorage(mitk::NavigationToolStorage::New(GetDataStorage()),
                                  "TrackingToolbox Default Storage");
  m_Controls->m_TrackingToolsStatusWidget->RemoveStatusLabels();
  QString toolLabel = QString("<none>");
  m_Controls->m_ToolLabel->setText(toolLabel);
  m_ToolStorageFilename = "";

  RemoveAllToolProjections();
}

void ARStrokeTreatmentView::OnStartStopTracking()
{
  if (!m_connected)
  {
    MITK_WARN << "Can't start tracking if no device is connected. Aborting";
    return;
  }
  if (m_tracking)
  {
    OnStopTracking();
  }
  else
  {
    OnStartTracking();
  }
}

void ARStrokeTreatmentView::OnFreezeUnfreezeTracking()
{
  if (m_Controls->m_FreezeUnfreezeTrackingButton->text() == "Freeze Tracking")
  {
    m_Worker->GetTrackingDeviceSource()->Freeze();
    m_Controls->m_FreezeUnfreezeTrackingButton->setText("Unfreeze Tracking");
    m_Controls->m_TrackingFrozenLabel->setVisible(true);
  }
  else if (m_Controls->m_FreezeUnfreezeTrackingButton->text() == "Unfreeze Tracking")
  {
    m_Worker->GetTrackingDeviceSource()->UnFreeze();
    m_Controls->m_FreezeUnfreezeTrackingButton->setText("Freeze Tracking");
    m_Controls->m_TrackingFrozenLabel->setVisible(false);
  }
}

void ARStrokeTreatmentView::ShowToolProjection(int index)
{
  mitk::DataNode::Pointer toolnode = m_toolStorage->GetTool(index)->GetDataNode();
  QString ToolProjectionName = "ToolProjection" + QString::number(index);
  m_ToolProjectionNode = this->GetDataStorage()->GetNamedNode(ToolProjectionName.toStdString());
  // If node does not exist, create the node for the Pointset
  if (m_ToolProjectionNode.IsNull())
  {
    m_ToolProjectionNode = mitk::DataNode::New();
    m_ToolProjectionNode->SetName(ToolProjectionName.toStdString());
    if (index < static_cast<int>(m_NeedleProjectionFilter->GetNumberOfInputs()))
    {
      m_NeedleProjectionFilter->SelectInput(index);
      m_NeedleProjectionFilter->Update();
      m_ToolProjectionNode->SetData(m_NeedleProjectionFilter->GetProjection());

      m_ToolProjectionNode->SetBoolProperty("show contour", true);
      this->GetDataStorage()->Add(m_ToolProjectionNode, toolnode);
    }
    //  this->FireNodeSelected(node);
  }
  else
  {
    m_ToolProjectionNode->SetBoolProperty("show contour", true);
  }
}

void ARStrokeTreatmentView::RemoveAllToolProjections()
{
  for (size_t i = 0; i < m_toolStorage->GetToolCount(); i++)
  {
    QString toolProjectionName = "ToolProjection" + QString::number(i);

    mitk::DataNode::Pointer node = this->GetDataStorage()->GetNamedNode(toolProjectionName.toStdString());

    // Deactivate and hide the tool projection
    if (!node.IsNull())
    {
      this->GetDataStorage()->Remove(node);
    }
  }
}

void ARStrokeTreatmentView::SelectToolProjection(int idx)
{
  if (m_Controls->showHideToolProjectionCheckBox->isChecked())
  {
    // Deactivate and hide the tool projection
    if (!m_ToolProjectionNode.IsNull())
    {
      this->GetDataStorage()->Remove(m_ToolProjectionNode);
    }

    if (m_NeedleProjectionFilter.IsNotNull())
    {
      m_NeedleProjectionFilter->Update();
    }
    // Refresh the view and the status widget
    mitk::RenderingManager::GetInstance()->RequestUpdateAll();
    // Show the tool projection for the currently selected tool
    ShowToolProjection(idx);
  }
}

void ARStrokeTreatmentView::OnShowHideToolProjectionClicked()
{
  int index = m_Controls->m_toolselector->currentIndex();
  // Activate and show the tool projection
  if (m_Controls->showHideToolProjectionCheckBox->isChecked())
  {
    ShowToolProjection(index);
    m_Controls->showHideToolAxisCheckBox->setEnabled(true);
  }
  else
  {
    RemoveAllToolProjections();
    m_Controls->showHideToolAxisCheckBox->setEnabled(false);
  }
  if (m_NeedleProjectionFilter->GetNumberOfInputs())
  {
    m_NeedleProjectionFilter->Update();
  }
  // Refresh the view and the status widget
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
  //  m_Controls->m_TrackingToolsStatusWidget->Refresh();
}

void ARStrokeTreatmentView::OnShowHideToolAxisClicked()
{
  if (!m_ShowHideToolAxis)
  {
    // Activate and show the tool axis
    m_NeedleProjectionFilter->ShowToolAxis(true);
    m_ShowHideToolAxis = true;
  }
  else
  {
    // Deactivate and hide the tool axis
    m_NeedleProjectionFilter->ShowToolAxis(false);
    m_NeedleProjectionFilter->GetProjection()->RemovePointIfExists(2);
    m_ShowHideToolAxis = false;
  }
  // Update the filter
  if (m_NeedleProjectionFilter->GetNumberOfInputs())
  {
    m_NeedleProjectionFilter->Update();
  }
  // Refresh the view and the status widget
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
  // m_Controls->m_TrackingToolsStatusWidget->Refresh();
}

void ARStrokeTreatmentView::OnConnectDisconnect()
{
  if (m_connected)
  {
    OnDisconnect();
  }
  else
  {
    OnConnect();
  }
}

void ARStrokeTreatmentView::OnConnect()
{
  MITK_DEBUG << "Connect Clicked";
  // check if everything is ready to start tracking
  if (this->m_toolStorage.IsNull())
  {
    MessageBox("Error: No Tools Loaded Yet!");
    return;
  }
  else if (this->m_toolStorage->GetToolCount() == 0)
  {
    MessageBox("Error: No Way To Track Without Tools!");
    return;
  }

  // parse tracking device data
  mitk::TrackingDeviceData data = mitk::UnspecifiedTrackingTypeInformation::GetDeviceDataUnspecified();
  QString qstr = m_Controls->m_VolumeSelectionBox->currentText();
  if ((!qstr.isNull()) || (!qstr.isEmpty()))
  {
    std::string str = qstr.toStdString();
    data = m_DeviceTypeCollection->GetDeviceDataByName(str); // Data will be set later, after device generation
  }

  //! [Thread 4]
  // initialize worker thread
  m_Worker->SetWorkerMethod(ARStrokeTreatmentTracking::eConnectDevice);
  m_Worker->SetTrackingDevice(this->m_Controls->m_ConfigurationWidget->GetTrackingDevice());
  m_Worker->SetInverseMode(m_Controls->m_InverseMode->isChecked());
  m_Worker->SetNavigationToolStorage(this->m_toolStorage);
  m_Worker->SetTrackingDeviceData(data);
  // start worker thread
  m_WorkerThread->start();
  //! [Thread 4]

  // enable checkboxes for projection and tool axis
  m_Controls->showHideToolAxisCheckBox->setEnabled(true);
  m_Controls->showHideToolProjectionCheckBox->setEnabled(true);
  m_Controls->m_toolselector->setEnabled(true);

  // disable buttons
  this->m_Controls->m_MainWidget->setEnabled(false);
}

void ARStrokeTreatmentView::EnableDisableTimerButtons(int enable)
{
  bool enableBool = enable;
  m_Controls->m_UpdateRateOptionsGroupBox->setEnabled(!enableBool);
  m_Controls->m_RenderWarningLabel->setVisible(enableBool);
}

void ARStrokeTreatmentView::OnConnectFinished(bool success, QString errorMessage)
{
  m_WorkerThread->quit();
  m_WorkerThread->wait();

  // enable buttons
  this->m_Controls->m_MainWidget->setEnabled(true);

  if (!success)
  {
    MITK_WARN << errorMessage.toStdString();
    MessageBox(errorMessage.toStdString());
    return;
  }

  //! [Thread 6]
  // get data from worker thread
  m_TrackingDeviceData = m_Worker->GetTrackingDeviceData();
  m_ToolVisualizationFilter = m_Worker->GetToolVisualizationFilter();
  if (m_ToolVisualizationFilter.IsNotNull())
  {
    // Connect the NeedleProjectionFilter to the ToolVisualizationFilter as third filter of the IGT pipeline
    m_NeedleProjectionFilter->ConnectTo(m_ToolVisualizationFilter);
    if (m_Controls->showHideToolProjectionCheckBox->isChecked())
    {
      ShowToolProjection(m_Controls->m_toolselector->currentIndex());
    }
  }

  //! [Thread 6]

  // enable/disable Buttons
  DisableOptionsButtons();
  DisableTrackingConfigurationButtons();

  m_Controls->m_TrackingControlLabel->setText("Status: connected");
  m_Controls->m_ConnectDisconnectButton->setText("Disconnect");
  m_Controls->m_ConnectSimpleMode->setText("Disconnect");
  m_Controls->m_StartStopTrackingButton->setEnabled(true);
  m_Controls->m_StartTrackingSimpleMode->setEnabled(true);
  m_connected = true;

  // During connection, thi sourceID of the tool storage changed. However, Microservice can't be updated on a different
  // thread. UpdateMicroservice is necessary to use filter to get the right storage belonging to a source. Don't do it
  // before m_connected is true, as we don't want to call content of OnToolStorageChanged.
  m_toolStorage->UpdateMicroservice();
}

void ARStrokeTreatmentView::OnDisconnect()
{
  m_Worker->SetWorkerMethod(ARStrokeTreatmentTracking::eDisconnectDevice);
  m_WorkerThread->start();
  m_Controls->m_MainWidget->setEnabled(false);
}

void ARStrokeTreatmentView::OnDisconnectFinished(bool success, QString errorMessage)
{
  m_WorkerThread->quit();
  m_WorkerThread->wait();
  m_Controls->m_MainWidget->setEnabled(true);

  if (!success)
  {
    MITK_WARN << errorMessage.toStdString();
    MessageBox(errorMessage.toStdString());
    return;
  }

  // enable/disable Buttons
  m_Controls->m_StartStopTrackingButton->setEnabled(false);
  m_Controls->m_StartTrackingSimpleMode->setEnabled(false);
  EnableOptionsButtons();
  EnableTrackingConfigurationButtons();
  m_Controls->m_TrackingControlLabel->setText("Status: disconnected");
  m_Controls->m_ConnectDisconnectButton->setText("Connect");
  m_Controls->m_ConnectSimpleMode->setText("Connect");
  m_Controls->m_FreezeUnfreezeTrackingButton->setText("Freeze Tracking");
  m_Controls->m_TrackingFrozenLabel->setVisible(false);
  m_connected = false;
}

void ARStrokeTreatmentView::OnStartTracking()
{
  // show tracking volume
  this->OnTrackingVolumeChanged(m_Controls->m_VolumeSelectionBox->currentText());
  // Reset the view to a defined start. Do it here and not in OnStartTrackingFinished, to give other tracking devices
  // the chance to reset the view to a different direction.
  this->GlobalReinit();
  m_Worker->SetWorkerMethod(ARStrokeTreatmentTracking::eStartTracking);
  m_WorkerThread->start();
  this->m_Controls->m_MainWidget->setEnabled(false);
}

void ARStrokeTreatmentView::OnStartTrackingFinished(bool success, QString errorMessage)
{
  //! [Thread 5]
  m_WorkerThread->quit();
  m_WorkerThread->wait();
  //! [Thread 5]
  this->m_Controls->m_MainWidget->setEnabled(true);

  if (!success)
  {
    MessageBox(errorMessage.toStdString());
    MITK_WARN << errorMessage.toStdString();
    return;
  }

  if (!(m_Controls->m_DisableAllTimers->isChecked()))
  {
    if (m_Controls->m_UseDifferentUpdateRates->isChecked())
    {
      if (m_Controls->m_RenderUpdateRate->value() != 0)
        m_TrackingRenderTimer->start(1000 / (m_Controls->m_RenderUpdateRate->value()));
      m_TrackingLoggingTimer->start(1000 / (m_Controls->m_LogUpdateRate->value()));
    }
    else
    {
      m_TrackingRenderTimer->start(1000 / (m_Controls->m_UpdateRate->value()));
      m_TrackingLoggingTimer->start(1000 / (m_Controls->m_UpdateRate->value()));
    }
  }

  m_Controls->m_TrackingControlLabel->setText("Status: tracking");

  // connect the tool visualization widget
  for (std::size_t i = 0; i < m_Worker->GetTrackingDeviceSource()->GetNumberOfOutputs(); i++)
  {
    m_Controls->m_TrackingToolsStatusWidget->AddNavigationData(m_Worker->GetTrackingDeviceSource()->GetOutput(i));
  }
  m_Controls->m_TrackingToolsStatusWidget->ShowStatusLabels();
  if (m_Controls->m_ShowToolQuaternions->isChecked())
  {
    m_Controls->m_TrackingToolsStatusWidget->SetShowQuaternions(true);
  }
  else
  {
    m_Controls->m_TrackingToolsStatusWidget->SetShowQuaternions(false);
  }

  // if activated enable open IGT link microservice
  if (m_Controls->m_EnableOpenIGTLinkMicroService->isChecked())
  {
    // create convertion filter
    m_IGTLConversionFilter = mitk::NavigationDataToIGTLMessageFilter::New();
    m_IGTLConversionFilter->SetName("IGT Tracking Toolbox");
    QString dataModeSelection = this->m_Controls->m_OpenIGTLinkDataFormat->currentText();
    if (dataModeSelection == "TDATA")
    {
      m_IGTLConversionFilter->SetOperationMode(mitk::NavigationDataToIGTLMessageFilter::ModeSendTDataMsg);
    }
    else if (dataModeSelection == "TRANSFORM")
    {
      m_IGTLConversionFilter->SetOperationMode(mitk::NavigationDataToIGTLMessageFilter::ModeSendTransMsg);
    }
    else if (dataModeSelection == "QTDATA")
    {
      m_IGTLConversionFilter->SetOperationMode(mitk::NavigationDataToIGTLMessageFilter::ModeSendQTDataMsg);
    }
    else if (dataModeSelection == "POSITION")
    {
      m_IGTLConversionFilter->SetOperationMode(mitk::NavigationDataToIGTLMessageFilter::ModeSendQTransMsg);
    }
    m_IGTLConversionFilter->ConnectTo(m_ToolVisualizationFilter);
    m_IGTLConversionFilter->RegisterAsMicroservice();

    // create server and message provider
    m_IGTLServer = mitk::IGTLServer::New(false);
    m_IGTLServer->SetName("Tracking Toolbox IGTL Server");
    m_IGTLMessageProvider = mitk::IGTLMessageProvider::New();
    m_IGTLMessageProvider->SetIGTLDevice(m_IGTLServer);
    m_IGTLMessageProvider->RegisterAsMicroservice();
  }

  m_tracking = true;
  m_Controls->m_ConnectDisconnectButton->setEnabled(false);
  m_Controls->m_StartStopTrackingButton->setText("Stop Tracking");
  m_Controls->m_StartTrackingSimpleMode->setText("Stop\nTracking");
  m_Controls->m_FreezeUnfreezeTrackingButton->setEnabled(true);
}

void ARStrokeTreatmentView::OnStopTracking()
{
  if (!m_tracking)
    return;
  for (unsigned int i = 0; i < m_ToolVisualizationFilter->GetNumberOfIndexedOutputs(); i++)
  {
    mitk::NavigationData::Pointer currentTool = m_ToolVisualizationFilter->GetOutput(i);
    if (currentTool->IsDataValid())
    {
      this->m_toolStorage->GetTool(i)->GetDataNode()->SetColor(mitk::IGTColor_INVALID);
    }
  }

  // refresh view and status widget
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();

  m_TrackingRenderTimer->stop();
  m_TrackingLoggingTimer->stop();

  m_Worker->SetWorkerMethod(ARStrokeTreatmentTracking::eStopTracking);
  m_WorkerThread->start();
  m_Controls->m_MainWidget->setEnabled(false);
}

void ARStrokeTreatmentView::OnStopTrackingFinished(bool success, QString errorMessage)
{
  m_WorkerThread->quit();
  m_WorkerThread->wait();
  m_Controls->m_MainWidget->setEnabled(true);
  if (!success)
  {
    MessageBox(errorMessage.toStdString());
    MITK_WARN << errorMessage.toStdString();
    return;
  }

  m_Controls->m_TrackingControlLabel->setText("Status: connected");
  if (m_logging)
    StopLogging();
  m_Controls->m_TrackingToolsStatusWidget->RemoveStatusLabels();
  m_Controls->m_TrackingToolsStatusWidget->PreShowTools(m_toolStorage);
  m_tracking = false;
  m_Controls->m_StartStopTrackingButton->setText("Start Tracking");
  m_Controls->m_StartTrackingSimpleMode->setText("Start\nTracking");
  m_Controls->m_ConnectDisconnectButton->setEnabled(true);
  m_Controls->m_FreezeUnfreezeTrackingButton->setEnabled(false);

  // unregister open IGT link micro service
  if (m_Controls->m_EnableOpenIGTLinkMicroService->isChecked())
  {
    m_IGTLConversionFilter->UnRegisterMicroservice();
    m_IGTLMessageProvider->UnRegisterMicroservice();
  }
}

void ARStrokeTreatmentView::OnTrackingDeviceChanged()
{
  mitk::TrackingDeviceType Type;

  if (m_Controls->m_ConfigurationWidget->GetTrackingDevice().IsNotNull())
  {
    Type = m_Controls->m_ConfigurationWidget->GetTrackingDevice()->GetType();
    // enable controls because device is valid
    m_Controls->m_TrackingToolsFrame->setEnabled(true);
    m_Controls->m_TrackingControlsFrame->setEnabled(true);
  }
  else
  {
    Type = mitk::UnspecifiedTrackingTypeInformation::GetTrackingDeviceName();
    MessageBox("Error: This tracking device is not included in this project. Please make sure that the device is "
               "installed and activated in your MITK build.");
    m_Controls->m_TrackingToolsFrame->setEnabled(false);
    m_Controls->m_TrackingControlsFrame->setEnabled(false);
    return;
  }

  // Code to enable/disable device specific buttons
  if (m_Controls->m_ConfigurationWidget->GetTrackingDevice()->AutoDetectToolsAvailable())
  {
    m_Controls->m_AutoDetectTools->setVisible(true);
  }
  else
  {
    m_Controls->m_AutoDetectTools->setVisible(false);
  }

  m_Controls->m_AddSingleTool->setEnabled(
    this->m_Controls->m_ConfigurationWidget->GetTrackingDevice()->AddSingleToolIsAvailable());

  // Code to select appropriate tracking volume for current type
  std::vector<mitk::TrackingDeviceData> Compatibles = m_DeviceTypeCollection->GetDeviceDataForLine(Type);
  m_Controls->m_VolumeSelectionBox->clear();
  for (std::size_t i = 0; i < Compatibles.size(); i++)
  {
    m_Controls->m_VolumeSelectionBox->addItem(Compatibles[i].Model.c_str());
  }
}

void ARStrokeTreatmentView::OnTrackingVolumeChanged(QString qstr)
{
  if (qstr.isNull())
    return;
  if (qstr.isEmpty())
    return;

  mitk::TrackingVolumeGenerator::Pointer volumeGenerator = mitk::TrackingVolumeGenerator::New();

  std::string str = qstr.toStdString();

  mitk::TrackingDeviceData data = m_DeviceTypeCollection->GetDeviceDataByName(str);
  m_TrackingDeviceData = data;

  volumeGenerator->SetTrackingDeviceData(data);
  volumeGenerator->Update();

  mitk::Surface::Pointer volumeSurface = volumeGenerator->GetOutput();

  m_TrackingVolumeNode->SetData(volumeSurface);

  if (!m_Controls->m_ShowTrackingVolume->isChecked())
    m_TrackingVolumeNode->SetOpacity(0.0);
  else
    m_TrackingVolumeNode->SetOpacity(0.25);

  GlobalReinit();
}

void ARStrokeTreatmentView::OnShowTrackingVolumeChanged()
{
  if (m_Controls->m_ShowTrackingVolume->isChecked())
  {
    OnTrackingVolumeChanged(m_Controls->m_VolumeSelectionBox->currentText());
    m_TrackingVolumeNode->SetOpacity(0.25);
  }
  else
  {
    m_TrackingVolumeNode->SetOpacity(0.0);
  }
}

void ARStrokeTreatmentView::OnAutoDetectTools()
{
  if (m_Controls->m_ConfigurationWidget->GetTrackingDevice()->AutoDetectToolsAvailable())
  {
    DisableTrackingConfigurationButtons();
    m_Worker->SetWorkerMethod(ARStrokeTreatmentTracking::eAutoDetectTools);
    m_Worker->SetTrackingDevice(m_Controls->m_ConfigurationWidget->GetTrackingDevice().GetPointer());
    m_Worker->SetDataStorage(this->GetDataStorage());
    m_WorkerThread->start();
    m_TimeoutTimer->start(30000);
    // disable controls until worker thread is finished
    this->m_Controls->m_MainWidget->setEnabled(false);
  }
}

void ARStrokeTreatmentView::OnAutoDetectToolsFinished(bool success, QString errorMessage)
{
  // Check, if the thread is running. There might have been a timeOut inbetween and this causes crashes...
  if (m_WorkerThread->isRunning())
  {
    m_TimeoutTimer->stop();
    m_WorkerThread->quit();
    m_WorkerThread->wait();
  }

  // enable controls again
  this->m_Controls->m_MainWidget->setEnabled(true);
  EnableTrackingConfigurationButtons();

  if (!success)
  {
    MITK_WARN << errorMessage.toStdString();
    MessageBox(errorMessage.toStdString());
    EnableTrackingConfigurationButtons();
    return;
  }

  mitk::NavigationToolStorage::Pointer autoDetectedStorage = m_Worker->GetNavigationToolStorage();

  // save detected tools
  std::string _autoDetectText;
  _autoDetectText = "Autodetected ";
  _autoDetectText.append(
    this->m_TrackingDeviceData.Line); // This is the device name as string of the current TrackingDevice.
  _autoDetectText.append(" Storage");
  this->ReplaceCurrentToolStorage(autoDetectedStorage, _autoDetectText);
  // auto save the new storage to hard disc (for persistence)
  AutoSaveToolStorage();
  // update label
  QString toolLabel =
    QString("Loaded Tools: ") + QString::number(m_toolStorage->GetToolCount()) + " Tools (Auto Detected)";
  m_Controls->m_ToolLabel->setText(toolLabel);
  // update tool preview
  m_Controls->m_TrackingToolsStatusWidget->RemoveStatusLabels();
  m_Controls->m_TrackingToolsStatusWidget->PreShowTools(m_toolStorage);

  EnableTrackingConfigurationButtons();

  // print a logging message about the detected tools
  switch (m_toolStorage->GetToolCount())
  {
    case 0:
      MITK_INFO("IGT Tracking Toolbox") << "Found no tools. Empty ToolStorage was autosaved to "
                                        << m_ToolStorageFilename.toStdString();
      break;
    case 1:
      MITK_INFO("IGT Tracking Toolbox") << "Found one tool. ToolStorage was autosaved to "
                                        << m_ToolStorageFilename.toStdString();
      break;
    default:
      MITK_INFO("IGT Tracking Toolbox") << "Found " << m_toolStorage->GetToolCount()
                                        << " tools. ToolStorage was autosaved to "
                                        << m_ToolStorageFilename.toStdString();
  }
}

void ARStrokeTreatmentView::MessageBox(std::string s)
{
  QMessageBox msgBox;
  msgBox.setText(s.c_str());
  msgBox.exec();
}

void ARStrokeTreatmentView::UpdateRenderTrackingTimer()
{
  // update filter
  m_ToolVisualizationFilter->Update();
  MITK_DEBUG << "Number of outputs ToolVisualizationFilter: " << m_ToolVisualizationFilter->GetNumberOfIndexedOutputs();
  MITK_DEBUG << "Number of inputs ToolVisualizationFilter: " << m_ToolVisualizationFilter->GetNumberOfIndexedInputs();

  // update tool colors to show tool status
  for (unsigned int i = 0; i < m_ToolVisualizationFilter->GetNumberOfIndexedOutputs(); i++)
  {
    mitk::NavigationData::Pointer currentTool = m_ToolVisualizationFilter->GetOutput(i);
    if (currentTool->IsDataValid())
    {
      this->m_toolStorage->GetTool(i)->GetDataNode()->SetColor(mitk::IGTColor_VALID);
    }
    else
    {
      this->m_toolStorage->GetTool(i)->GetDataNode()->SetColor(mitk::IGTColor_WARNING);
    }
  }

  // Update the NeedleProjectionFilter
  if (m_NeedleProjectionFilter.IsNotNull())
  {
    m_NeedleProjectionFilter->Update();
  }

  // refresh view and status widget
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
  m_Controls->m_TrackingToolsStatusWidget->Refresh();
}

void ARStrokeTreatmentView::UpdateLoggingTrackingTimer()
{
  // update logging
  if (m_logging)
  {
    this->m_loggingFilter->Update();
    m_loggedFrames = this->m_loggingFilter->GetNumberOfRecordedSteps();
    this->m_Controls->m_LoggedFramesLabel->setText("Logged Frames: " + QString::number(m_loggedFrames));
    // check if logging stopped automatically
    if ((m_loggedFrames > 1) && (!m_loggingFilter->GetRecording()))
    {
      StopLogging();
    }
  }
  // refresh status widget
  m_Controls->m_TrackingToolsStatusWidget->Refresh();
}

void ARStrokeTreatmentView::OnChooseFileClicked()
{
  QDir currentPath = QFileInfo(m_Controls->m_LoggingFileName->text()).dir();

  // if no path was selected (QDir would select current working dir then) or the
  // selected path does not exist -> use home directory
  if (currentPath == QDir() || !currentPath.exists())
  {
    currentPath = QDir(QDir::homePath());
  }

  QString filename =
    QFileDialog::getSaveFileName(nullptr, tr("Choose Logging File"), currentPath.absolutePath(), "*.*");
  if (filename == "")
    return;
  this->m_Controls->m_LoggingFileName->setText(filename);
  this->OnToggleFileExtension();
}
// bug-16470: toggle file extension after clicking on radio button
void ARStrokeTreatmentView::OnToggleFileExtension()
{
  QString currentInputText = this->m_Controls->m_LoggingFileName->text();
  QString currentFile = QFileInfo(currentInputText).baseName();
  QDir currentPath = QFileInfo(currentInputText).dir();
  if (currentFile.isEmpty())
  {
    currentFile = "logfile";
  }
  // Setting currentPath to default home path when currentPath is empty or it does not exist
  if (currentPath == QDir() || !currentPath.exists())
  {
    currentPath = QDir::homePath();
  }
  // check if csv radio button is clicked
  if (this->m_Controls->m_CsvFormat->isChecked())
  {
    // you needn't add a seperator to the input text when currentpath is the rootpath
    if (currentPath.isRoot())
    {
      this->m_Controls->m_LoggingFileName->setText(QDir::toNativeSeparators(currentPath.absolutePath()) + currentFile +
                                                   ".csv");
    }

    else
    {
      this->m_Controls->m_LoggingFileName->setText(QDir::toNativeSeparators(currentPath.absolutePath()) +
                                                   QDir::separator() + currentFile + ".csv");
    }
  }
  // check if xml radio button is clicked
  else if (this->m_Controls->m_XmlFormat->isChecked())
  {
    // you needn't add a seperator to the input text when currentpath is the rootpath
    if (currentPath.isRoot())
    {
      this->m_Controls->m_LoggingFileName->setText(QDir::toNativeSeparators(currentPath.absolutePath()) + currentFile +
                                                   ".xml");
    }
    else
    {
      this->m_Controls->m_LoggingFileName->setText(QDir::toNativeSeparators(currentPath.absolutePath()) +
                                                   QDir::separator() + currentFile + ".xml");
    }
  }
}

void ARStrokeTreatmentView::OnToggleAdvancedSimpleMode()
{
  if (m_SimpleModeEnabled)
  {
    m_Controls->m_simpleWidget->setVisible(false);
    m_Controls->m_MainWidget->setVisible(true);
    m_Controls->m_SimpleUI->setChecked(false);
    m_SimpleModeEnabled = false;
  }
  else
  {
    m_Controls->m_simpleWidget->setVisible(true);
    m_Controls->m_MainWidget->setVisible(false);
    m_SimpleModeEnabled = true;
  }
}

void ARStrokeTreatmentView::OnToggleDifferentUpdateRates()
{
  if (m_Controls->m_UseDifferentUpdateRates->isChecked())
  {
    if (m_Controls->m_RenderUpdateRate->value() == 0)
      m_Controls->m_RenderWarningLabel->setVisible(true);
    else
      m_Controls->m_RenderWarningLabel->setVisible(false);

    m_Controls->m_UpdateRate->setEnabled(false);
    m_Controls->m_OptionsUpdateRateLabel->setEnabled(false);

    m_Controls->m_RenderUpdateRate->setEnabled(true);
    m_Controls->m_OptionsRenderUpdateRateLabel->setEnabled(true);

    m_Controls->m_LogUpdateRate->setEnabled(true);
    m_Controls->m_OptionsLogUpdateRateLabel->setEnabled(true);
  }

  else
  {
    m_Controls->m_RenderWarningLabel->setVisible(false);

    m_Controls->m_UpdateRate->setEnabled(true);
    m_Controls->m_OptionsUpdateRateLabel->setEnabled(true);

    m_Controls->m_RenderUpdateRate->setEnabled(false);
    m_Controls->m_OptionsRenderUpdateRateLabel->setEnabled(false);

    m_Controls->m_LogUpdateRate->setEnabled(false);
    m_Controls->m_OptionsLogUpdateRateLabel->setEnabled(false);
  }
}

void ARStrokeTreatmentView::OnChangeRenderUpdateRate()
{
  if (m_Controls->m_RenderUpdateRate->value() == 0)
    m_Controls->m_RenderWarningLabel->setVisible(true);
  else
    m_Controls->m_RenderWarningLabel->setVisible(false);
}

void ARStrokeTreatmentView::StartLogging()
{
  if (m_ToolVisualizationFilter.IsNull())
  {
    MessageBox("Cannot activate logging without a connected device. Configure and connect a tracking device first.");
    return;
  }

  if (!m_logging)
  {
    // initialize logging filter
    m_loggingFilter = mitk::NavigationDataRecorder::New();
    m_loggingFilter->SetRecordOnlyValidData(m_Controls->m_SkipInvalidData->isChecked());

    m_loggingFilter->ConnectTo(m_ToolVisualizationFilter);

    if (m_Controls->m_LoggingLimit->isChecked())
    {
      m_loggingFilter->SetRecordCountLimit(m_Controls->m_LoggedFramesLimit->value());
    }

    // start filter with try-catch block for exceptions
    try
    {
      m_loggingFilter->StartRecording();
    }
    catch (mitk::IGTException &)
    {
      std::string errormessage = "Error during start recording. Recorder already started recording?";
      QMessageBox::warning(nullptr, "IGTPlayer: Error", errormessage.c_str());
      m_loggingFilter->StopRecording();
      return;
    }

    // update labels / logging variables
    this->m_Controls->m_LoggingLabel->setText("Logging ON");
    this->m_Controls->m_LoggedFramesLabel->setText("Logged Frames: 0");
    m_loggedFrames = 0;
    m_logging = true;
    DisableLoggingButtons();
  }
}

void ARStrokeTreatmentView::StopLogging()
{
  if (m_logging)
  {
    // stop logging
    m_loggingFilter->StopRecording();
    m_logging = false;

    // update GUI
    this->m_Controls->m_LoggingLabel->setText("Logging OFF");
    EnableLoggingButtons();

    // write the results to a file
    if (m_Controls->m_CsvFormat->isChecked())
    {
      mitk::IOUtil::Save(m_loggingFilter->GetNavigationDataSet(),
                         this->m_Controls->m_LoggingFileName->text().toStdString());
    }
    else if (m_Controls->m_XmlFormat->isChecked())
    {
      mitk::IOUtil::Save(m_loggingFilter->GetNavigationDataSet(),
                         this->m_Controls->m_LoggingFileName->text().toStdString());
    }
  }
}

void ARStrokeTreatmentView::OnAddSingleTool()
{
  QString Identifier = "Tool#";
  QString Name = "NewTool";
  if (m_toolStorage.IsNotNull())
  {
    Identifier += QString::number(m_toolStorage->GetToolCount());
    Name += QString::number(m_toolStorage->GetToolCount());
  }
  else
  {
    Identifier += "0";
    Name += "0";
  }
  m_Controls->m_NavigationToolCreationWidget->Initialize(
    GetDataStorage(), Identifier.toStdString(), Name.toStdString());
  m_Controls->m_NavigationToolCreationWidget->SetTrackingDeviceType(
    m_Controls->m_ConfigurationWidget->GetTrackingDevice()->GetType(), false);
  m_Controls->m_TrackingToolsWidget->setCurrentIndex(1);

  // disable tracking volume during tool editing
  lastTrackingVolumeState = m_Controls->m_ShowTrackingVolume->isChecked();
  if (lastTrackingVolumeState)
    m_Controls->m_ShowTrackingVolume->click();
  GlobalReinit();
}

void ARStrokeTreatmentView::OnAddSingleToolFinished()
{
  m_Controls->m_TrackingToolsWidget->setCurrentIndex(0);
  if (this->m_toolStorage.IsNull())
  {
    // this shouldn't happen!
    MITK_WARN << "No ToolStorage available, cannot add tool, aborting!";
    return;
  }
  m_toolStorage->AddTool(m_Controls->m_NavigationToolCreationWidget->GetCreatedTool());
  m_Controls->m_TrackingToolsStatusWidget->PreShowTools(m_toolStorage);
  m_Controls->m_ToolLabel->setText("<manually added>");

  // displya in tool selector
  // m_Controls->m_toolselector->addItem(QString::fromStdString(m_Controls->m_NavigationToolCreationWidget->GetCreatedTool()->GetToolName()));

  // auto save current storage for persistence
  MITK_INFO << "Auto saving manually added tools for persistence.";
  AutoSaveToolStorage();

  // enable tracking volume again
  if (lastTrackingVolumeState)
    m_Controls->m_ShowTrackingVolume->click();
  GlobalReinit();
}

void ARStrokeTreatmentView::OnAddSingleToolCanceled()
{
  m_Controls->m_TrackingToolsWidget->setCurrentIndex(0);

  // enable tracking volume again
  if (lastTrackingVolumeState)
    m_Controls->m_ShowTrackingVolume->click();
  GlobalReinit();
}

void ARStrokeTreatmentView::GlobalReinit()
{
  // get all nodes that have not set "includeInBoundingBox" to false
  mitk::NodePredicateNot::Pointer pred = mitk::NodePredicateNot::New(
    mitk::NodePredicateProperty::New("includeInBoundingBox", mitk::BoolProperty::New(false)));

  mitk::DataStorage::SetOfObjects::ConstPointer rs = this->GetDataStorage()->GetSubset(pred);
  // calculate bounding geometry of these nodes
  auto bounds = this->GetDataStorage()->ComputeBoundingGeometry3D(rs, "visible");

  // initialize the views to the bounding geometry
  mitk::RenderingManager::GetInstance()->InitializeViews(bounds);
}

void ARStrokeTreatmentView::DisableLoggingButtons()
{
  m_Controls->m_StartLogging->setEnabled(false);
  m_Controls->m_LoggingFileName->setEnabled(false);
  m_Controls->m_ChooseFile->setEnabled(false);
  m_Controls->m_LoggingLimit->setEnabled(false);
  m_Controls->m_LoggedFramesLimit->setEnabled(false);
  m_Controls->m_CsvFormat->setEnabled(false);
  m_Controls->m_XmlFormat->setEnabled(false);
  m_Controls->m_SkipInvalidData->setEnabled(false);
  m_Controls->m_StopLogging->setEnabled(true);
}

void ARStrokeTreatmentView::EnableLoggingButtons()
{
  m_Controls->m_StartLogging->setEnabled(true);
  m_Controls->m_LoggingFileName->setEnabled(true);
  m_Controls->m_ChooseFile->setEnabled(true);
  m_Controls->m_LoggingLimit->setEnabled(true);
  m_Controls->m_LoggedFramesLimit->setEnabled(true);
  m_Controls->m_CsvFormat->setEnabled(true);
  m_Controls->m_XmlFormat->setEnabled(true);
  m_Controls->m_SkipInvalidData->setEnabled(true);
  m_Controls->m_StopLogging->setEnabled(false);
}

void ARStrokeTreatmentView::DisableOptionsButtons()
{
  m_Controls->m_ShowTrackingVolume->setEnabled(false);
  m_Controls->m_UseDifferentUpdateRates->setEnabled(false);
  m_Controls->m_UpdateRate->setEnabled(false);
  m_Controls->m_OptionsUpdateRateLabel->setEnabled(false);
  m_Controls->m_RenderUpdateRate->setEnabled(false);
  m_Controls->m_OptionsRenderUpdateRateLabel->setEnabled(false);
  m_Controls->m_LogUpdateRate->setEnabled(false);
  m_Controls->m_OptionsLogUpdateRateLabel->setEnabled(false);
  m_Controls->m_DisableAllTimers->setEnabled(false);
  m_Controls->m_OtherOptionsGroupBox->setEnabled(false);
  m_Controls->m_EnableOpenIGTLinkMicroService->setEnabled(false);
  m_Controls->m_OpenIGTLinkDataFormat->setEnabled(false);
}

void ARStrokeTreatmentView::EnableOptionsButtons()
{
  m_Controls->m_ShowTrackingVolume->setEnabled(true);
  m_Controls->m_UseDifferentUpdateRates->setEnabled(true);
  m_Controls->m_DisableAllTimers->setEnabled(true);
  m_Controls->m_OtherOptionsGroupBox->setEnabled(true);
  m_Controls->m_EnableOpenIGTLinkMicroService->setEnabled(true);
  m_Controls->m_OpenIGTLinkDataFormat->setEnabled(true);
  OnToggleDifferentUpdateRates();
}

void ARStrokeTreatmentView::EnableTrackingControls()
{
  m_Controls->m_TrackingControlsFrame->setEnabled(true);
}

void ARStrokeTreatmentView::DisableTrackingControls()
{
  m_Controls->m_TrackingControlsFrame->setEnabled(false);
}

void ARStrokeTreatmentView::EnableTrackingConfigurationButtons()
{
  m_Controls->m_AutoDetectTools->setEnabled(true);
  m_Controls->m_AddSingleTool->setEnabled(
    this->m_Controls->m_ConfigurationWidget->GetTrackingDevice()->AddSingleToolIsAvailable());
  m_Controls->m_LoadTools->setEnabled(true);
  m_Controls->m_ResetTools->setEnabled(true);
}

void ARStrokeTreatmentView::DisableTrackingConfigurationButtons()
{
  m_Controls->m_AutoDetectTools->setEnabled(false);
  m_Controls->m_AddSingleTool->setEnabled(false);
  m_Controls->m_LoadTools->setEnabled(false);
  m_Controls->m_ResetTools->setEnabled(false);
}

void ARStrokeTreatmentView::ReplaceCurrentToolStorage(mitk::NavigationToolStorage::Pointer newStorage,
                                                      std::string newStorageName)
{
  // first: get rid of the old one
  // don't reset if there is no tool storage. BugFix #17793
  if (m_toolStorage.IsNotNull())
  {
    m_toolStorage->UnLockStorage(); // only to be sure...
    m_toolStorage->UnRegisterMicroservice();
    m_toolStorage = nullptr;
  }

  // now: replace by the new one
  m_toolStorage = newStorage;
  m_toolStorage->SetName(newStorageName);
  m_toolStorage->RegisterAsMicroservice();
}

void ARStrokeTreatmentView::OnTimeOut()
{
  MITK_WARN << "TimeOut. Quitting the thread...";
  m_WorkerThread->quit();
  // only if we can't quit use terminate.
  if (!m_WorkerThread->wait(1000))
  {
    MITK_ERROR << "Can't quit the thread. Terminating... Might cause further problems, be careful!";
    m_WorkerThread->terminate();
    m_WorkerThread->wait();
  }
  m_TimeoutTimer->stop();
}

void ARStrokeTreatmentView::OnToolStorageChanged(const ctkServiceEvent event)
{
  // don't listen to any changes during connection, toolStorage is locked anyway, so this are only changes of e.g.
  // sourceID which are not relevant for the widget.
  if (!m_connected && (event.getType() == ctkServiceEvent::MODIFIED))
  {
    m_Controls->m_ConfigurationWidget->OnToolStorageChanged();

    m_Controls->m_toolselector->clear();
    for (size_t i = 0; i < m_toolStorage->GetToolCount(); i++)
    {
      m_Controls->m_toolselector->addItem(QString::fromStdString(m_toolStorage->GetTool(i)->GetToolName()));
    }
  }
}

//! [StoreUISettings]
void ARStrokeTreatmentView::StoreUISettings()
{
  // persistence service does not directly work in plugins for now
  // -> using QSettings
  QSettings settings;

  settings.beginGroup(QString::fromStdString(VIEW_ID));
  MITK_DEBUG << "Store UI settings";
  // set the values of some widgets and attrbutes to the QSettings
  settings.setValue("ShowTrackingVolume", QVariant(m_Controls->m_ShowTrackingVolume->isChecked()));
  settings.setValue("toolStorageFilename", QVariant(m_ToolStorageFilename));
  settings.setValue("VolumeSelectionBox", QVariant(m_Controls->m_VolumeSelectionBox->currentIndex()));
  settings.setValue("SimpleModeEnabled", QVariant(m_SimpleModeEnabled));
  settings.setValue("OpenIGTLinkDataFormat", QVariant(m_Controls->m_OpenIGTLinkDataFormat->currentIndex()));
  settings.setValue("EnableOpenIGTLinkMicroService",
                    QVariant(m_Controls->m_EnableOpenIGTLinkMicroService->isChecked()));
  settings.endGroup();
}
//! [StoreUISettings]

//! [LoadUISettings]
void ARStrokeTreatmentView::LoadUISettings()
{
  // persistence service does not directly work in plugins for now -> using QSettings
  QSettings settings;

  settings.beginGroup(QString::fromStdString(VIEW_ID));

  // set some widgets and attributes by the values from the QSettings
  m_Controls->m_ShowTrackingVolume->setChecked(settings.value("ShowTrackingVolume", false).toBool());
  m_Controls->m_EnableOpenIGTLinkMicroService->setChecked(
    settings.value("EnableOpenIGTLinkMicroService", true).toBool());
  m_Controls->m_VolumeSelectionBox->setCurrentIndex(settings.value("VolumeSelectionBox", 0).toInt());
  m_Controls->m_OpenIGTLinkDataFormat->setCurrentIndex(settings.value("OpenIGTLinkDataFormat", 0).toInt());
  m_ToolStorageFilename = settings.value("toolStorageFilename", QVariant("")).toString();
  if (settings.value("SimpleModeEnabled", false).toBool())
  {
    this->OnToggleAdvancedSimpleMode();
  }
  settings.endGroup();
  //! [LoadUISettings]

  //! [LoadToolStorage]
  // try to deserialize the tool storage from the given tool storage file name
  if (!m_ToolStorageFilename.isEmpty())
  {
    // try-catch block for exceptions
    try
    {
      mitk::NavigationToolStorageDeserializer::Pointer myDeserializer =
        mitk::NavigationToolStorageDeserializer::New(GetDataStorage());
      m_toolStorage->UnRegisterMicroservice();
      m_toolStorage = myDeserializer->Deserialize(m_ToolStorageFilename.toStdString());
      m_toolStorage->RegisterAsMicroservice();

      // update label
      UpdateToolStorageLabel(m_ToolStorageFilename);

      // update tool preview
      m_Controls->m_TrackingToolsStatusWidget->RemoveStatusLabels();
      m_Controls->m_TrackingToolsStatusWidget->PreShowTools(m_toolStorage);
    }
    catch (const mitk::IGTException &e)
    {
      MITK_WARN("QmitkMITKIGTTrackingToolBoxView")
        << "Error during restoring tools. Problems with file (" << m_ToolStorageFilename.toStdString()
        << "), please check the file? Error message: " << e.GetDescription();
      this->OnResetTools(); // if there where errors reset the tool storage to avoid problems later on
    }
  }
  //! [LoadToolStorage]
}

void ARStrokeTreatmentView::UpdateToolStorageLabel(QString pathOfLoadedStorage)
{
  QFileInfo myPath(pathOfLoadedStorage); // use this to seperate filename from path
  QString toolLabel = myPath.fileName();
  if (toolLabel.size() > 45) // if the tool storage name is to long trimm the string
  {
    toolLabel.resize(40);
    toolLabel += "[...]";
  }
  m_Controls->m_ToolLabel->setText(toolLabel);
}

void ARStrokeTreatmentView::AutoSaveToolStorage()
{
  m_ToolStorageFilename = m_AutoSaveFilename;
  mitk::NavigationToolStorageSerializer::Pointer mySerializer = mitk::NavigationToolStorageSerializer::New();
  mySerializer->Serialize(m_ToolStorageFilename.toStdString(), m_toolStorage);
}
