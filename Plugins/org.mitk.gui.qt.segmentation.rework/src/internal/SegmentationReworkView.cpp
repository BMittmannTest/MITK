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
#include "SegmentationReworkView.h"

#include <itkFileTools.h>

// Qt
#include <QMessageBox>
#include <QProgressBar>

#include <mitkIOUtil.h>
#include <QmitkNewSegmentationDialog.h>
#include <mitkSegTool2D.h>
#include <mitkToolManagerProvider.h>
#include <mitkOverwriteSliceImageFilter.h>
#include <mitkSegmentationInterpolationController.h>
#include <mitkImagePixelWriteAccessor.h>
#include <mitkExtractSliceFilter.h>
#include <mitkVtkImageOverwrite.h>

#include <mitkDICOMQIIOMimeTypes.h>
#include <mitkDICOMDCMTKTagScanner.h>

const std::string SegmentationReworkView::VIEW_ID = "org.mitk.views.segmentationreworkview";

void SegmentationReworkView::SetFocus() {}

void SegmentationReworkView::CreateQtPartControl(QWidget *parent)
{
  // create GUI widgets from the Qt Designer's .ui file
  m_Controls.setupUi(parent);
  m_Parent = parent;
  counter = 0;

  qRegisterMetaType< std::vector<std::string> >("std::vector<std::string>");

  //m_Controls.verticalWidget->setVisible(false);
  m_Controls.cleanDicomBtn->setVisible(false);
  m_Controls.individualWidget_2->setVisible(false);

  m_Controls.sliderWidget->setMinimum(1);
  m_Controls.sliderWidget->setMaximum(100);
  m_Controls.sliderWidget->setTickInterval(1);
  m_Controls.sliderWidget->setSingleStep(1);
  m_Controls.radioA->setChecked(true);

  connect(m_Controls.buttonUpload, &QPushButton::clicked, this, &SegmentationReworkView::UploadNewSegmentation);
  connect(m_Controls.buttonNewSeg, &QPushButton::clicked, this, &SegmentationReworkView::CreateNewSegmentationC);
  connect(m_Controls.cleanDicomBtn, &QPushButton::clicked, this, &SegmentationReworkView::CleanDicomFolder);
  connect(m_Controls.restartConnection, &QPushButton::clicked, this, &SegmentationReworkView::RestartConnection);
  connect(m_Controls.checkIndiv, &QCheckBox::stateChanged, this, &SegmentationReworkView::OnIndividualCheckChange);
  connect(m_Controls.sliderWidget, &ctkSliderWidget::valueChanged, this, &SegmentationReworkView::OnSliderWidgetChanged);

  m_DownloadBaseDir = mitk::IOUtil::GetTempPathA() + "segrework";
  MITK_INFO << "using download base dir: " << m_DownloadBaseDir;
  m_UploadBaseDir = mitk::IOUtil::GetTempPathA() + "uploadSeg";

  if (!itksys::SystemTools::FileIsDirectory(m_DownloadBaseDir))
  {
    itk::FileTools::CreateDirectory(m_DownloadBaseDir);
  }

  if (!itksys::SystemTools::FileIsDirectory(m_UploadBaseDir))
  {
    itk::FileTools::CreateDirectory(m_UploadBaseDir);
  }

  utility::string_t port = U("2020");
  utility::string_t address = U("http://127.0.0.1:");
  address.append(port);

  m_HttpHandler = std::unique_ptr<SegmentationReworkREST>(new SegmentationReworkREST(address));

  connect(m_HttpHandler.get(),
    &SegmentationReworkREST::InvokeUpdateChartWidget,
    this,
    &SegmentationReworkView::UpdateChartWidget);
  connect(this, &SegmentationReworkView::InvokeLoadData, this, &SegmentationReworkView::LoadData);

  connect(this, &SegmentationReworkView::InvokeProgress, this, &SegmentationReworkView::AddProgress);

  m_HttpHandler->SetPutCallback(std::bind(&SegmentationReworkView::RESTPutCallback, this, std::placeholders::_1));
  m_HttpHandler->SetGetImageSegCallback(std::bind(&SegmentationReworkView::RESTGetCallback, this, std::placeholders::_1));
  m_HttpHandler->SetGetAddSeriesCallback(std::bind(&SegmentationReworkView::RESTGetCallbackGeneric, this, std::placeholders::_1));
  m_HttpHandler->Open().wait();

  MITK_INFO << "Listening for requests at: " << utility::conversions::to_utf8string(address);

  //utility::string_t pacsURL = U("http://jip-dktk/dcm4chee-arc/aets/DCM4CHEE");
  utility::string_t pacsURL = U("http://193.174.48.78:8090/dcm4chee-arc/aets/DCM4CHEE");
  utility::string_t restURL = U("http://localhost:8000");

  m_DICOMWeb = new mitk::DICOMWeb(pacsURL);
  MITK_INFO << "requests to pacs are sent to: " << utility::conversions::to_utf8string(pacsURL);
  m_Controls.dcm4cheeURL->setText({ (utility::conversions::to_utf8string(pacsURL).c_str()) });

  m_RestService = new mitk::RESTClient(restURL);
}

void SegmentationReworkView::OnSliderWidgetChanged(double value)
{
  std::map<double, double>::iterator it;
  unsigned int count = 0;
  for (it = m_ScoreMap.begin(); it != m_ScoreMap.end(); it++)
  {
    if (it->second < value)
    {
      count++;
    }
  }
  QString labelsToDelete = "slices to delete: " + QString::number(count);
  m_Controls.slicesToDeleteLabel->setText(labelsToDelete);

  std::map<double, double> thresholdMap;

  for (it = m_ScoreMap.begin(); it != m_ScoreMap.end(); it++)
  {
    thresholdMap.insert(std::map<double, double>::value_type(it->first, value));
  }

  m_Controls.chartWidget->RemoveData(m_thresholdLabel);
  m_Controls.chartWidget->AddData2D(thresholdMap, m_thresholdLabel);
  m_Controls.chartWidget->SetChartType(m_thresholdLabel, QmitkChartWidget::ChartType::line);
  m_Controls.chartWidget->Show();
}

void SegmentationReworkView::AddProgress(int progress, QString status)
{
  auto futureValue = m_Controls.progressBar->value() + progress;
  if (futureValue >= 100) 
  {
    m_Controls.progressBar->setValue(0);
    m_Controls.progressBar->setFormat("");
  } else 
  {
    m_Controls.progressBar->setFormat(status.append(" %p%"));
    m_Controls.progressBar->setValue(futureValue);
  }
}

void SegmentationReworkView::RestartConnection()
{
  auto host = m_Controls.dcm4cheeHostValue->text();
  std::string url = host.toStdString() + "/dcm4chee-arc/aets/DCM4CHEE";

  m_DICOMWeb = new mitk::DICOMWeb(utility::conversions::to_string_t(url));
  MITK_INFO << "requests to pacs are sent to: " << url;
  m_Controls.dcm4cheeURL->setText({ (utility::conversions::to_utf8string(url).c_str()) });
}

void SegmentationReworkView::OnIndividualCheckChange(int state)
{
  if(state == Qt::Unchecked)
  {
    m_Controls.individualWidget_2->setVisible(false);
  }
  else if (state == Qt::Checked)
  {
    m_Controls.individualWidget_2->setVisible(true);
  }
}

void SegmentationReworkView::RESTPutCallback(const SegmentationReworkREST::DicomDTO &dto)
{
  emit InvokeProgress(20, { "display graph and query structured report" });
  SetSimilarityGraph(dto.simScoreArray, dto.minSliceStart);

  m_SRUID = dto.srSeriesUID;
  m_GroundTruth = dto.groundTruth;

  typedef std::map<std::string, std::string> ParamMap;
  ParamMap seriesInstancesParams;

  seriesInstancesParams.insert((ParamMap::value_type({"StudyInstanceUID"}, dto.studyUID)));
  seriesInstancesParams.insert((ParamMap::value_type({"SeriesInstanceUID"}, dto.srSeriesUID)));
  seriesInstancesParams.insert((ParamMap::value_type({"includefield"}, {"0040A375"}))); // Current Requested Procedure Evidence Sequence

  try {
    m_DICOMWeb->QuidoRSInstances(seriesInstancesParams).then([=](web::json::value jsonResult) 
    {

      auto firstResult = jsonResult[0];
      auto actualListKey = firstResult.at(U("0040A375")).as_object().at(U("Value")).as_array()[0].as_object().at(U("00081115")).as_object().at(U("Value")).as_array();

      std::string segSeriesUIDA = "";
      std::string segSeriesUIDB = "";
      std::string imageSeriesUID = "";

      for (unsigned int index = 0; index < actualListKey.size(); index++) {
        auto element = actualListKey[index].as_object();
        // get SOP class UID
        auto innerElement = element.at(U("00081199")).as_object().at(U("Value")).as_array()[0];
        auto sopClassUID = innerElement.at(U("00081150")).as_object().at(U("Value")).as_array()[0].as_string();

        auto seriesUID = utility::conversions::to_utf8string(element.at(U("0020000E")).as_object().at(U("Value")).as_array()[0].as_string());

        if (sopClassUID == U("1.2.840.10008.5.1.4.1.1.66.4")) // SEG
        {
          if (segSeriesUIDA.length() == 0)
          {
            segSeriesUIDA = seriesUID;
          }
          else
          {
            segSeriesUIDB = seriesUID;
          }
        }
        else if (sopClassUID == U("1.2.840.10008.5.1.4.1.1.2"))  // CT
        {
          imageSeriesUID = seriesUID;
        }
      }

      emit InvokeProgress(10, {"load composite context of structured report"});
      MITK_INFO << "image series UID " << imageSeriesUID;
      MITK_INFO << "seg A series UID " << segSeriesUIDA;
      MITK_INFO << "seg B series UID " << segSeriesUIDB;

      MITK_INFO << "Load related dicom series ...";

      std::string folderPathSeries = mitk::IOUtil::CreateTemporaryDirectory("XXXXXX", m_DownloadBaseDir) +"/";

      std::string pathSegA = mitk::IOUtil::CreateTemporaryDirectory("XXXXXX", m_DownloadBaseDir) + "/";
      std::string pathSegB = mitk::IOUtil::CreateTemporaryDirectory("XXXXXX", m_DownloadBaseDir) + "/";

      auto folderPathSegA = utility::conversions::to_string_t(pathSegA);
      auto folderPathSegB = utility::conversions::to_string_t(pathSegB);

      m_CurrentStudyUID = dto.studyUID;

      std::vector<pplx::task<std::string>> tasks;
      auto imageSeriesTask = m_DICOMWeb->WadoRS(utility::conversions::to_string_t(folderPathSeries), dto.studyUID, imageSeriesUID);
      auto segATask = m_DICOMWeb->WadoRS(folderPathSegA, dto.studyUID, segSeriesUIDA);
      auto segBTask = m_DICOMWeb->WadoRS(folderPathSegB, dto.studyUID, segSeriesUIDB);
      tasks.push_back(imageSeriesTask);
      tasks.push_back(segATask);
      tasks.push_back(segBTask);

      auto joinTask = pplx::when_all(begin(tasks), end(tasks));
      auto filePathList = joinTask.then([&](std::vector<std::string> filePathList) {
        emit InvokeProgress(50, {"load dicom files from disk"});
        InvokeLoadData(filePathList);
      });
    });

  }
  catch (mitk::Exception&  e) {
    MITK_ERROR << e.what();
  }
}

void SegmentationReworkView::RESTGetCallbackGeneric(const SegmentationReworkREST::DicomDTO &dto)
{
  std::vector<pplx::task<std::string>> tasks;

  if (dto.seriesUIDList.size() > 0) {
    for (std::string segSeriesUID : dto.seriesUIDList)
    {
      std::string folderPathSeries = mitk::IOUtil::CreateTemporaryDirectory("XXXXXX", m_DownloadBaseDir) + "/";
      try {
        auto seriesTask = m_DICOMWeb->WadoRS(utility::conversions::to_string_t(folderPathSeries), dto.studyUID, segSeriesUID);
        tasks.push_back(seriesTask);
      } catch (const mitk::Exception &exception) 
      {
        MITK_INFO << exception.what();
        return;
      }
    }
  }

  try 
  {
    auto joinTask = pplx::when_all(begin(tasks), end(tasks));
    auto filePathList = joinTask.then([&](std::vector<std::string> filePathList) { InvokeLoadData(filePathList); });
  }
  catch (const mitk::Exception& exception) 
  {
    MITK_INFO << exception.what();
    return;
  }
}

void SegmentationReworkView::RESTGetCallback(const SegmentationReworkREST::DicomDTO &dto) 
{
  std::string folderPathSeries = mitk::IOUtil::CreateTemporaryDirectory("XXXXXX", m_DownloadBaseDir) + "/";

  MITK_INFO << folderPathSeries;

  std::string pathSeg = mitk::IOUtil::CreateTemporaryDirectory("XXXXXX", m_DownloadBaseDir) + "/";
  auto folderPathSeg = utility::conversions::to_string_t(pathSeg);

  MITK_INFO << pathSeg;

  try {
    std::vector<pplx::task<std::string>> tasks;
    auto imageSeriesTask = m_DICOMWeb->WadoRS(utility::conversions::to_string_t(folderPathSeries), dto.studyUID, dto.imageSeriesUID);

    tasks.push_back(imageSeriesTask);

    if (dto.seriesUIDList.size() > 0) {
      for (std::string segSeriesUID : dto.seriesUIDList)
      {
        auto segTask = m_DICOMWeb->WadoRS(folderPathSeg, dto.studyUID, segSeriesUID);
        tasks.push_back(segTask);
      }
    }
    else {
      auto segATask = m_DICOMWeb->WadoRS(folderPathSeg, dto.studyUID, dto.segSeriesUIDA);
      tasks.push_back(segATask);
    }


    auto joinTask = pplx::when_all(begin(tasks), end(tasks));
    auto filePathList = joinTask.then([&](std::vector<std::string> filePathList) { InvokeLoadData(filePathList); });
  }
  catch (const mitk::Exception &exception) {
    MITK_INFO << exception.what();
  }
}

std::string SegmentationReworkView::GetAlgorithmOfSegByPath(std::string path)
{
  auto scanner = mitk::DICOMDCMTKTagScanner::New();

  mitk::DICOMTagPath algorithmName;
  algorithmName.AddAnySelection(0x0062, 0x0002).AddElement(0x0062, 0x0009);

  mitk::StringList files;
  files.push_back(path);
  scanner->SetInputFiles(files);
  scanner->AddTagPath(algorithmName);

  scanner->Scan();

  mitk::DICOMDatasetAccessingImageFrameList frames = scanner->GetFrameInfoList();
  auto findings = frames.front()->GetTagValueAsString(algorithmName);
  if (findings.size() != 0)
    MITK_INFO << findings.front().value;
  return findings.front().value;
}

void SegmentationReworkView::LoadData(std::vector<std::string> filePathList)
{
  MITK_INFO << "Loading finished. Pushing data to data storage ...";
  auto ds = GetDataStorage();
  auto dataNodes = mitk::IOUtil::Load(filePathList, *ds);
  // reinit view
  mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(ds);

  // find data nodes
  m_Image = dataNodes->at(0);
  m_Image->SetName("image data");
  m_SegA = dataNodes->at(1);
  m_SegB = dataNodes->at(2);

  auto algorithmNameA = GetAlgorithmOfSegByPath(filePathList[1]);
  auto algorithmNameB = GetAlgorithmOfSegByPath(filePathList[2]);
  m_SegA->SetName(algorithmNameA);
  m_SegB->SetName(algorithmNameB);
  m_Controls.labelSegAValue->setText(algorithmNameA.c_str());
  m_Controls.labelSegBValue->setText(algorithmNameB.c_str());
  m_Controls.labelGroundTruthValue->setText(m_GroundTruth.c_str());
  emit InvokeProgress(20, { "" });
}

void SegmentationReworkView::UpdateChartWidget()
{
  m_Controls.chartWidget->Show();
}

void SegmentationReworkView::SetSimilarityGraph(std::vector<double> simScoreArray, int sliceMinStart)
{
  std::string label = "similarity graph";
  m_thresholdLabel = "threshold";
  //m_Controls.chartWidget->Clear();

  double sliceIndex = sliceMinStart;
  for (double score : simScoreArray)
  {
    m_ScoreMap.insert(std::map<double, double>::value_type(sliceIndex, score));
    sliceIndex++;
  }

  std::map<double, double> thresholdMap;

  m_Controls.chartWidget->AddData2D(m_ScoreMap, label);
  m_Controls.chartWidget->AddData2D(thresholdMap, m_thresholdLabel);
  m_Controls.chartWidget->SetChartType(label, QmitkChartWidget::ChartType::line);
  m_Controls.chartWidget->SetChartType(m_thresholdLabel, QmitkChartWidget::ChartType::line);
  m_Controls.chartWidget->SetXAxisLabel("slice number");
  m_Controls.chartWidget->SetYAxisLabel("similarity in percent");
  m_Controls.chartWidget->SetTitle("Similartiy Score for Segmentation Comparison");
}

void SegmentationReworkView::UploadNewSegmentation()
{
  AddProgress(10, { "save SEG to temp folder" });
  std::string folderPathSeg = mitk::IOUtil::CreateTemporaryDirectory("XXXXXX", m_UploadBaseDir) + "/";

  const std::string savePath = folderPathSeg + m_SegC->GetName() + ".dcm";
  const std::string mimeType = mitk::MitkDICOMQIIOMimeTypes::DICOMSEG_MIMETYPE_NAME();
  mitk::IOUtil::Save(m_SegC->GetData(), mimeType, savePath);

  // get Series Instance UID from new SEG
  auto scanner = mitk::DICOMDCMTKTagScanner::New();
  mitk::DICOMTagPath seriesUID(0x0020, 0x000E);

  mitk::StringList files;
  files.push_back(savePath);
  scanner->SetInputFiles(files);
  scanner->AddTagPath(seriesUID);

  scanner->Scan();

  mitk::DICOMDatasetAccessingImageFrameList frames = scanner->GetFrameInfoList();
  auto findings = frames.front()->GetTagValueAsString(seriesUID);
  auto segSeriesUID = findings.front().value;

  AddProgress(20, {"push SEG to PACS"});
  auto filePath = utility::conversions::to_string_t(savePath);
  try {
    m_DICOMWeb->StowRS(filePath, m_CurrentStudyUID).then([=] {
      emit InvokeProgress(40, {"persist reworked SEG to evaluation database"}); 

      MitkUriBuilder queryBuilder(U("tasks/evaluations/"));
      queryBuilder.append_query(U("srUID"), utility::conversions::to_string_t(m_SRUID));
      //queryBuilder.append_path(U("71/"));
      m_RestService->Get(queryBuilder.to_string()).then([=](web::json::value result) {
        MITK_INFO << "after GET";
        MITK_INFO << utility::conversions::to_utf8string(result.to_string());
        auto updatedContent = result.as_array()[0];
        updatedContent[U("reworkedSegmentationUID")] = web::json::value::string(utility::conversions::to_string_t(segSeriesUID));
     
        auto id = updatedContent.at(U("id")).as_integer();
        MITK_INFO << id;
        auto idParam = std::to_string(id).append("/");

        MitkUriBuilder queryBuilder(U("tasks/evaluations"));
        queryBuilder.append_path(utility::conversions::to_string_t(idParam));

        m_RestService->PUT(queryBuilder.to_string(), updatedContent).then([=](web::json::value result) {
          MITK_INFO << "successfully stored";
          emit InvokeProgress(30, { "successfully stored" });
        });
      });
    });
  }
  catch (const std::exception &exception)
  {
    std::cout << exception.what() << std::endl;
  }
}

std::vector<unsigned int> SegmentationReworkView::CreateSegmentation(mitk::Image::Pointer baseSegmentation, double threshold)
{
  MITK_INFO << "handle individual segmentation creation";
  std::map<double, double>::iterator it;

  std::vector<unsigned int> sliceIndices;

  unsigned int count = 0;
  for (it = m_ScoreMap.begin(); it != m_ScoreMap.end(); it++)
  {
    if (it->second < threshold)
    {
      auto index = it->first;
      try 
      {
        mitk::ImagePixelWriteAccessor<unsigned short, 3> imageAccessor(baseSegmentation);
        for (unsigned int x = 0; x < baseSegmentation->GetDimension(0); x++)
        {
          for (unsigned int y = 0; y < baseSegmentation->GetDimension(1); y++)
          {
            imageAccessor.SetPixelByIndex({ { x, y, int(index) } }, 0);
          }
        }
      }
      catch (mitk::Exception& e) {
        MITK_ERROR << e.what();
      }

      count++;
      sliceIndices.push_back(index);
      MITK_INFO << "slice " << it->first <<  " removed ";
	   }
  }
  MITK_INFO << "slices deleted " << count;
  return sliceIndices;
}

void SegmentationReworkView::CreateNewSegmentationC() 
{
  mitk::ToolManager* toolManager = mitk::ToolManagerProvider::GetInstance()->GetToolManager();
  toolManager->InitializeTools();
  toolManager->SetReferenceData(m_Image);
  
  mitk::Image::Pointer baseImage;
  if (m_Controls.radioA->isChecked()) {
    baseImage = dynamic_cast<mitk::Image*>(m_SegA->GetData())->Clone();
  }
  else if (m_Controls.radioB->isChecked()) {
    baseImage = dynamic_cast<mitk::Image*>(m_SegB->GetData())->Clone();
  }

  if (m_Controls.checkIndiv->isChecked())
  {
    auto sliceIndices = CreateSegmentation(baseImage, m_Controls.sliderWidget->value());
  }

  QmitkNewSegmentationDialog* dialog = new QmitkNewSegmentationDialog(m_Parent); // needs a QWidget as parent, "this" is not QWidget

  int dialogReturnValue = dialog->exec();
  if (dialogReturnValue == QDialog::Rejected)
  {
    // user clicked cancel or pressed Esc or something similar
    return;
  }

  // ask the user about an organ type and name, add this information to the image's (!) propertylist
  // create a new image of the same dimensions and smallest possible pixel type
  mitk::Tool* firstTool = toolManager->GetToolById(0);
  if (firstTool)
  {
    try
    {
      std::string newNodeName = dialog->GetSegmentationName().toStdString();
      if (newNodeName.empty())
      {
        newNodeName = "no_name";
      }

      mitk::DataNode::Pointer newSegmentation = firstTool->CreateSegmentationNode(baseImage, newNodeName, dialog->GetColor());
      // initialize showVolume to false to prevent recalculating the volume while working on the segmentation
      newSegmentation->SetProperty("showVolume", mitk::BoolProperty::New(false));
      if (!newSegmentation)
      {
        return; // could be aborted by user
      }

      if (mitk::ToolManagerProvider::GetInstance()->GetToolManager()->GetWorkingData(0))
      {
        mitk::ToolManagerProvider::GetInstance()->GetToolManager()->GetWorkingData(0)->SetSelected(false);
      }
      newSegmentation->SetSelected(true);
      this->GetDataStorage()->Add(newSegmentation, toolManager->GetReferenceData(0)); // add as a child, because the segmentation "derives" from the original
     
      m_SegC = newSegmentation;

      auto referencedImages = m_Image->GetData()->GetProperty("files");
      m_SegC->GetData()->SetProperty("referenceFiles", referencedImages);

    }
    catch (std::bad_alloc)
    {
      QMessageBox::warning(nullptr, tr("Create new segmentation"), tr("Could not allocate memory for new segmentation"));
    }
  }
  else {
    MITK_INFO << "no tools...";
  }
}

void SegmentationReworkView::CleanDicomFolder() 
{
  if (m_SegA || m_SegB || m_SegC) {
    QMessageBox::warning(nullptr, tr("Clean dicom folder"), tr("Please remove the data in data storage before cleaning the download folder"));
    return;
  }

  //std::experimental::filesystem::remove_all(m_DownloadBaseDir);
  // TODO : use POCO
  //itk::FileTools::CreateDirectory(m_DownloadBaseDir);
}
