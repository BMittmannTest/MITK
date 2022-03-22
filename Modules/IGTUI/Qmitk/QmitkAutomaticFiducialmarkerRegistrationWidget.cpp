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

#include "QmitkAutomaticFiducialmarkerRegistrationWidget.h"

#include <QMessageBox>

#include "mitkProperties.h"
#include <mitkImage.h>
#include <mitkImageCast.h>
#include <mitkLabelSetImage.h>
#include <mitkNodePredicateNot.h>
#include <mitkNodePredicateProperty.h>

#include <mitkRenderingManager.h>
#include <mitkStaticIGTHelperFunctions.h>

#include <vtkLandmarkTransform.h>
#include <vtkPoints.h>

#include <itkGeometryUtilities.h>
#include <itkImageRegionIterator.h>
#include <mitkImagePixelReadAccessor.h>

#include <itkCannyEdgeDetectionImageFilter.h>
#include <itkSimpleContourExtractorImageFilter.h>
#include <itkZeroCrossingImageFilter.h>

QmitkAutomaticFiducialmarkerRegistrationWidget::QmitkAutomaticFiducialmarkerRegistrationWidget(QWidget *parent)
  : QWidget(parent),
    m_Controls(nullptr),
    m_TransformMarkerCSToSensorCS(mitk::AffineTransform3D::New()),
    m_TransformMarkerCSToImageCS(mitk::AffineTransform3D::New()),
    m_DataStorage(nullptr)
{
  this->UnsetFloatingImageGeometry();
  this->DefineDataStorageImageFilter();
  this->CreateQtPartControl(this);
}

QmitkAutomaticFiducialmarkerRegistrationWidget::~QmitkAutomaticFiducialmarkerRegistrationWidget()
{
  delete m_Controls;
}

void QmitkAutomaticFiducialmarkerRegistrationWidget::Initialize(itk::SmartPointer<mitk::DataStorage> dataStorage)
{
  MITK_INFO << "Initialize()";
  m_Controls->selectedImageComboBox->SetDataStorage(dataStorage);
  m_Controls->selectedSurfaceComboBox->SetDataStorage(dataStorage);
  m_Controls->pointSetsToChooseComboBox->SetDataStorage(dataStorage);
  this->m_DataStorage = dataStorage;
}

mitk::AffineTransform3D::Pointer QmitkAutomaticFiducialmarkerRegistrationWidget::GetTransformMarkerCSToImageCS()
{
  if (m_DataStorage.IsNull())
  {
    MITK_WARN << "Returning nullptr, as data storage was not set and transform was not determined yet.";
    return nullptr;
  }
  return this->m_TransformMarkerCSToImageCS;
}

mitk::AffineTransform3D::Pointer QmitkAutomaticFiducialmarkerRegistrationWidget::GetTransformMarkerCSToSensorCS()
{
  // The following calculations are related to the 3mm | 15mm fiducial configuration
  MITK_WARN << "Determining the transform Marker-CS --> Sensor-CS is supported actually only for the 3mm | 15mm "
               "fiducial configuration.";

  if (m_TransformMarkerCSToSensorCS.IsNull())
  {
    m_TransformMarkerCSToSensorCS = mitk::AffineTransform3D::New();
  }

  if (m_Controls->selectedTrackerComboBox->currentIndex() ==
      0) // Use the NDI disc tracker for performing the registration:
  {
    MITK_INFO << "Use NDI disc tracker for performing the CT-to-US-registration";
    mitk::Vector3D translationNDI;
    translationNDI[0] = 15.000;
    translationNDI[1] = 8.000;
    translationNDI[2] = -2.500; // use -2.500 instead of 0.000 for the z translation if working with v2 of the
                                // sensor adapter.

    m_TransformMarkerCSToSensorCS->SetOffset(translationNDI);

    // Quaternion (x, y, z, r) --> n = (0,0,1) --> q(0,0,sin(90°),cos(90°))
    mitk::Quaternion qNDI(
      0,
      0,
      1,
      0); // corresponding to a rotation of 180° around the normal z-axis.
          // .transpose() is needed for changing the rows and the columns of the returned rotation_matrix_transpose
    vnl_matrix_fixed<double, 3, 3> vnl_rotation = qNDI.rotation_matrix_transpose().transpose(); // :-)
    mitk::Matrix3D rotationMatrix;

    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        rotationMatrix[i][j] = vnl_rotation[i][j];
      }
    }

    m_TransformMarkerCSToSensorCS->SetMatrix(rotationMatrix);
  }
  else // Use the polhemus RX2 tracker for performing the registration:
  {
    MITK_INFO << "Use Polhemus RX2 tracker for performing the CT-to-US-registration";
    mitk::Vector3D translationPolhemus;
    translationPolhemus[0] = -18.175;
    translationPolhemus[1] = 15.000;
    translationPolhemus[2] = 10.501; // considering a distance from the base plate of 0.315 inch and not 0.313 inch
                                     // use 10.501 instead of 8.001 for the z translation if working with v2 of the
                                     // sensor adapter.

    m_TransformMarkerCSToSensorCS->SetOffset(translationPolhemus);

    // Quaternion (x, y, z, r) --> n = (1,0,0) --> q(sin(90°),0,0,cos(90°))
    mitk::Quaternion q1(
      1,
      0,
      0,
      0); // corresponding to a rotation of 180° around the normal x-axis.
          // .transpose() is needed for changing the rows and the columns of the returned rotation_matrix_transpose
    vnl_matrix_fixed<double, 3, 3> vnl_rotation = q1.rotation_matrix_transpose().transpose(); // :-)
    mitk::Matrix3D rotationMatrix;

    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        rotationMatrix[i][j] = vnl_rotation[i][j];
      }
    }

    m_TransformMarkerCSToSensorCS->SetMatrix(rotationMatrix);
  }
  // The transformation from the sensor-CS to the marker-CS is calculated now.
  return m_TransformMarkerCSToSensorCS;
}

mitk::AffineTransform3D::Pointer QmitkAutomaticFiducialmarkerRegistrationWidget::GetInverseTransform(
  mitk::AffineTransform3D::Pointer transform)
{
  mitk::AffineTransform3D::Pointer inverseTransform = mitk::AffineTransform3D::New();
  mitk::AffineTransform3D::Pointer inverse =
    dynamic_cast<mitk::AffineTransform3D *>(transform->GetInverseTransform().GetPointer());

  if (inverse.IsNull())
  {
    MITK_WARN << "Could not get inverse transform of mitk::AffineTransform3D. Returning nullptr";
    return nullptr;
  }

  inverseTransform->SetOffset(inverse->GetOffset());
  inverseTransform->SetMatrix(inverse->GetMatrix());

  return inverseTransform;
}

void QmitkAutomaticFiducialmarkerRegistrationWidget::CreateQtPartControl(QWidget *parent)
{
  if (!m_Controls)
  {
    m_Controls = new Ui::QmitkAutomaticFiducialmarkerRegistrationWidget;
  }
  m_Controls->setupUi(parent);
  m_Controls->selectedImageComboBox->SetPredicate(m_IsAPatientImagePredicate);
  m_Controls->selectedSurfaceComboBox->SetPredicate(m_IsASurfacePredicate);
  m_Controls->pointSetsToChooseComboBox->SetPredicate(m_IsAPointSetPredicate);

  // create signal/slot connections
  connect(m_Controls->selectedImageComboBox,
          SIGNAL(OnSelectionChanged(const mitk::DataNode *)),
          this,
          SLOT(OnImageToRegisterComboBoxSelectionChanged(const mitk::DataNode *)));
  connect(
    m_Controls->doRegistrationMarkerToImagePushButton, SIGNAL(clicked()), this, SLOT(OnRegisterMarkerCSToImageCS()));
  connect(m_Controls->localizeFiducialMarkerPushButton, SIGNAL(clicked()), this, SLOT(OnLocalizeFiducials()));
  connect(m_Controls->addPontSetPushButton, SIGNAL(clicked()), this, SLOT(OnAddPointSetClicked()));
  connect(m_Controls->removePointSetPushButton, SIGNAL(clicked()), this, SLOT(OnRemovePointSetClicked()));
}

void QmitkAutomaticFiducialmarkerRegistrationWidget::OnImageToRegisterComboBoxSelectionChanged(
  const mitk::DataNode *node)
{
  MITK_INFO << "OnImageToRegisterComboBoxSelectionChanged()";

  if (m_ImageToRegister.IsNotNull())
  {
    // TODO: Define, what will happen if the imageCT is not null...
  }

  if (node == nullptr)
  {
    this->UnsetFloatingImageGeometry();
    m_ImageToRegister = nullptr;
    return;
  }

  mitk::DataNode *selectedFloatingImage = m_Controls->selectedImageComboBox->GetSelectedNode();
  if (selectedFloatingImage == nullptr)
  {
    this->UnsetFloatingImageGeometry();
    m_ImageToRegister = nullptr;
    return;
  }

  mitk::Image::Pointer floatingImage = dynamic_cast<mitk::Image *>(selectedFloatingImage->GetData());
  if (floatingImage.IsNull())
  {
    MITK_WARN << "Failed to cast selected segmentation node to mitk::Image*";
    this->UnsetFloatingImageGeometry();
    m_ImageToRegister = nullptr;
    return;
  }

  m_ImageToRegister = floatingImage;
  this->SetImageGeometryInformation(floatingImage.GetPointer());
}

void QmitkAutomaticFiducialmarkerRegistrationWidget::OnRegisterMarkerCSToImageCS()
{
  this->CreateMarkerModelCoordinateSystemPointSet();

  // Check for initialization
  if (m_MarkerModelCoordinateSystemPointSet.IsNull() || m_MarkerImageToRegisterCoordinateSystemPointSet.IsNull())
  {
    MITK_WARN << "Fiducial Landmarks are not initialized yet, cannot register";
    return;
  }

  // Retrieve fiducials
  if (m_MarkerImageToRegisterCoordinateSystemPointSet->GetSize() != m_MarkerModelCoordinateSystemPointSet->GetSize())
  {
    MITK_WARN << "Not the same number of fiducials, cannot register";
    return;
  }
  else if (m_MarkerImageToRegisterCoordinateSystemPointSet->GetSize() < 3)
  {
    MITK_WARN << "Need at least 3 fiducials, cannot register";
    return;
  }

  //############### conversion to vtk data types (we will use the vtk landmark based transform)
  //########################## convert point sets to vtk poly data
  vtkSmartPointer<vtkPoints> sourcePoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkPoints> targetPoints = vtkSmartPointer<vtkPoints>::New();
  for (int i = 0; i < m_MarkerModelCoordinateSystemPointSet->GetSize(); i++)
  {
    double point[3] = {m_MarkerModelCoordinateSystemPointSet->GetPoint(i)[0],
                       m_MarkerModelCoordinateSystemPointSet->GetPoint(i)[1],
                       m_MarkerModelCoordinateSystemPointSet->GetPoint(i)[2]};
    sourcePoints->InsertNextPoint(point);

    double point_targets[3] = {m_MarkerImageToRegisterCoordinateSystemPointSet->GetPoint(i)[0],
                               m_MarkerImageToRegisterCoordinateSystemPointSet->GetPoint(i)[1],
                               m_MarkerImageToRegisterCoordinateSystemPointSet->GetPoint(i)[2]};
    targetPoints->InsertNextPoint(point_targets);
  }

  //########################### here, the actual transform is computed ##########################
  // compute transform
  vtkSmartPointer<vtkLandmarkTransform> transform = vtkSmartPointer<vtkLandmarkTransform>::New();
  transform->SetSourceLandmarks(sourcePoints);
  transform->SetTargetLandmarks(targetPoints);
  transform->SetModeToRigidBody();
  transform->Modified();
  transform->Update();
  // compute FRE of transform

  double FRE = mitk::StaticIGTHelperFunctions::ComputeFRE(
    m_MarkerModelCoordinateSystemPointSet, m_MarkerImageToRegisterCoordinateSystemPointSet, transform);
  MITK_INFO << "FRE: " << FRE << " mm";

  //#############################################################################################

  //############### conversion back to itk/mitk data types ##########################
  // convert from vtk to itk data types
  itk::Matrix<float, 3, 3> rotationFloat = itk::Matrix<float, 3, 3>();
  itk::Vector<float, 3> translationFloat = itk::Vector<float, 3>();
  itk::Matrix<double, 3, 3> rotationDouble = itk::Matrix<double, 3, 3>();
  itk::Vector<double, 3> translationDouble = itk::Vector<double, 3>();

  vtkSmartPointer<vtkMatrix4x4> m = transform->GetMatrix();
  for (int k = 0; k < 3; k++)
    for (int l = 0; l < 3; l++)
    {
      rotationFloat[k][l] = m->GetElement(k, l);
      rotationDouble[k][l] = m->GetElement(k, l);
    }
  for (int k = 0; k < 3; k++)
  {
    translationFloat[k] = m->GetElement(k, 3);
    translationDouble[k] = m->GetElement(k, 3);
  }
  // create mitk affine transform 3D and save it to the class member
  if (m_TransformMarkerCSToImageCS.IsNull())
  {
    m_TransformMarkerCSToImageCS = mitk::AffineTransform3D::New();
  }
  m_TransformMarkerCSToImageCS->SetMatrix(rotationDouble);
  m_TransformMarkerCSToImageCS->SetOffset(translationDouble);
  MITK_INFO << m_TransformMarkerCSToImageCS;

  mitk::PointSet *pointSet_orig = m_MarkerModelCoordinateSystemPointSet;
  mitk::PointSet::Pointer pointSet_moved = mitk::PointSet::New();

  for (int i = 0; i < pointSet_orig->GetSize(); i++)
  {
    pointSet_moved->InsertPoint(m_TransformMarkerCSToImageCS->TransformPoint(pointSet_orig->GetPoint(i)));
  }

  pointSet_orig->Clear();
  for (int i = 0; i < pointSet_moved->GetSize(); i++)
    pointSet_orig->InsertPoint(pointSet_moved->GetPoint(i));

  // Do a global reinit
  mitk::RenderingManager::GetInstance()->InitializeViewsByBoundingObjects(this->m_DataStorage);
}

mitk::DataNode::Pointer QmitkAutomaticFiducialmarkerRegistrationWidget::GetImageNode()
{
  return m_Controls->selectedImageComboBox->GetSelectedNode();
}

mitk::DataNode::Pointer QmitkAutomaticFiducialmarkerRegistrationWidget::GetSurfaceNode()
{
  return m_Controls->selectedSurfaceComboBox->GetSelectedNode();
}

mitk::DataNode::Pointer QmitkAutomaticFiducialmarkerRegistrationWidget::GetPointSetNode()
{
  return m_Controls->pointSetsToChooseComboBox->GetSelectedNode();
}

std::vector<mitk::DataNode::Pointer> QmitkAutomaticFiducialmarkerRegistrationWidget::GetPointSetsDataNodes()
{
  return m_PointSetsDataNodes;
}

void QmitkAutomaticFiducialmarkerRegistrationWidget::OnLocalizeFiducials()
{
  m_FiducialMarkerCentroids.clear();
  m_FiducialCandidates.clear();
  if (m_MarkerImageToRegisterCoordinateSystemPointSet.IsNotNull())
  {
    m_MarkerImageToRegisterCoordinateSystemPointSet->Clear();
  }

  if (!this->FilterImage())
  {
    QMessageBox msgBox;
    msgBox.setText("Cannot perform filtering of the image. The floating image = nullptr.");
    msgBox.exec();
    return;
  }
  mitk::AffineTransform3D::Pointer transform = m_ImageToRegister->GetGeometry()->GetIndexToWorldTransform();
  MITK_WARN << "IndexToWorldTransform_CTimage = " << transform;

  this->GetCentroidsOfLabeledObjects();

  int numberTrialsNonReducingFiducialCandidates = 0;
  for (int i = 0; (m_FiducialCandidates.size() > 8 &&
                   m_Controls->fiducialMarkerConfigurationComboBox->currentIndex() <= 1 && i < 100) ||
                  (m_FiducialCandidates.size() > 4 &&
                   m_Controls->fiducialMarkerConfigurationComboBox->currentIndex() == 2 && i < 100);
       ++i)
  {
    MITK_INFO << "Size centroids fiducial candidates: " << m_FiducialCandidates.size();
    unsigned int oldNumberFiducialCandidates = m_FiducialCandidates.size();
    if (m_Controls->fiducialMarkerConfigurationComboBox->currentIndex() <= 1)
    {
      if (!this->EliminateFiducialCandidatesByEuclideanDistances() || m_FiducialCandidates.size() < 8)
      {
        QMessageBox msgBox;
        QString text = QString("Have found %1 instead of 8 fiducial candidates.\
        Cannot perform fiducial localization procedure.")
                         .arg(m_FiducialCandidates.size());
        msgBox.setText(text);
        msgBox.exec();
        return;
      }
    }
    else if (m_Controls->fiducialMarkerConfigurationComboBox->currentIndex() == 2)
    {
      if (!this->EliminateFiducialCandidatesByEuclideanDistances() || m_FiducialCandidates.size() < 4)
      {
        QMessageBox msgBox;
        QString text = QString("Have found %1 instead of 4 fiducial candidates.\
        Cannot perform fiducial localization procedure.")
                         .arg(m_FiducialCandidates.size());
        msgBox.setText(text);
        msgBox.exec();
        return;
      }
    }

    if (oldNumberFiducialCandidates - m_FiducialCandidates.size() == 0)
    {
      ++numberTrialsNonReducingFiducialCandidates;
    }
    if (numberTrialsNonReducingFiducialCandidates >= 3)
    {
      MITK_INFO << "Skipping EliminateFiducialCandidatesByEuclideanDistances() as \
                    it does not reduce the fiducial candidates further.";
      break;
    }
  }

  // Find all pairs of fiducial candidates, whose euclidean distance is smaller than
  // a predefined limit value. As one fiducial candidate of each pair does not
  // correspond to a true fiducial, these false candidates can be eliminated, which
  // reduces the total amount of fiducial candidates.
  for (int i = 0; (m_FiducialCandidates.size() > 8 &&
                     m_Controls->fiducialMarkerConfigurationComboBox->currentIndex() <= 1 && (i < 50) ||
                   (m_FiducialCandidates.size() > 4 &&
                    m_Controls->fiducialMarkerConfigurationComboBox->currentIndex() == 2 && (i < 50)));
       ++i)
  {
    this->EliminateNearFiducialCandidatesByMaxDistanceToCentroids();
  }

  // Before calling NumerateFiducialMarks it must be sure,
  // that there rested only 4 or 8 fiducial candidates depending on the fiducial marker configuration.
  this->NumerateFiducialMarks();
}

void QmitkAutomaticFiducialmarkerRegistrationWidget::OnAddPointSetClicked()
{
  mitk::DataNode *selectedPointSet = m_Controls->pointSetsToChooseComboBox->GetSelectedNode();
  if (selectedPointSet == nullptr)
  {
    return;
  }

  mitk::PointSet::Pointer pointSet = dynamic_cast<mitk::PointSet *>(selectedPointSet->GetData());
  if (pointSet.IsNull())
  {
    MITK_WARN << "Failed to cast selected pointSet node to mitk::PointSet*";
    return;
  }
  QString name = QString::fromStdString(selectedPointSet->GetName());
  for (int counter = 0; counter < m_Controls->selectedPointSetsListWidget->count(); ++counter)
  {
    if (m_Controls->selectedPointSetsListWidget->item(counter)->text().compare(name) == 0)
    {
      MITK_INFO << "PointSet already exist in list of selected PointSets. Do not add the PointSet.";
      return;
    }
  }
  m_Controls->selectedPointSetsListWidget->addItem(name);
  m_PointSetsDataNodes.push_back(selectedPointSet);
}

void QmitkAutomaticFiducialmarkerRegistrationWidget::OnRemovePointSetClicked()
{
  int position = m_Controls->selectedPointSetsListWidget->currentRow();
  if (m_Controls->selectedPointSetsListWidget->count() == 0 || position < 0)
  {
    return;
  }

  m_PointSetsDataNodes.erase(m_PointSetsDataNodes.begin() + position);
  QListWidgetItem *item = m_Controls->selectedPointSetsListWidget->currentItem();
  m_Controls->selectedPointSetsListWidget->removeItemWidget(item);
  delete item;
}


void QmitkAutomaticFiducialmarkerRegistrationWidget::UnsetFloatingImageGeometry()
{
  m_ImageDimension[0] = 0;
  m_ImageDimension[1] = 0;
  m_ImageDimension[2] = 0;

  m_ImageSpacing[0] = 1;
  m_ImageSpacing[1] = 1;
  m_ImageSpacing[2] = 1;
}

void QmitkAutomaticFiducialmarkerRegistrationWidget::SetImageGeometryInformation(mitk::Image *image)
{
  m_ImageDimension[0] = image->GetDimension(0);
  m_ImageDimension[1] = image->GetDimension(1);
  m_ImageDimension[2] = image->GetDimension(2);

  m_ImageSpacing[0] = image->GetGeometry()->GetSpacing()[0];
  m_ImageSpacing[1] = image->GetGeometry()->GetSpacing()[1];
  m_ImageSpacing[2] = image->GetGeometry()->GetSpacing()[2];

  MITK_INFO << "GetScalarValueMax = " << image->GetScalarValueMax(0);
  MITK_INFO << "GetScalarValueMin = " << image->GetScalarValueMin(0);
}

double QmitkAutomaticFiducialmarkerRegistrationWidget::GetVoxelVolume()
{
  if (m_ImageToRegister.IsNull())
  {
    return 0.0;
  }

  MITK_INFO << "ImageSpacing = " << m_ImageSpacing;
  return m_ImageSpacing[0] * m_ImageSpacing[1] * m_ImageSpacing[2];
}

double QmitkAutomaticFiducialmarkerRegistrationWidget::GetFiducialVolume(double radius)
{
  return 1.333333333 * 3.141592 * (radius * radius * radius);
}

bool QmitkAutomaticFiducialmarkerRegistrationWidget::FilterImage()
{
  if (m_ImageToRegister.IsNull())
  {
    return false;
  }
  ImageType::Pointer itkImageOld;
  mitk::CastToItkImage(m_ImageToRegister, itkImageOld);

  ImageType::Pointer itkImage1 = ImageType::New();
  mitk::CastToItkImage(m_ImageToRegister, itkImage1);

  this->InitializeImageFilters();

  m_ThresholdFilter->SetInput(itkImage1);
  m_LaplacianFilter1->SetInput(m_ThresholdFilter->GetOutput());
  m_LaplacianFilter2->SetInput(m_LaplacianFilter1->GetOutput());

  // Dynamically adjust lower threshold of binary filter only for MRI images:
  if (m_Controls->selectedModalityComboBox->currentIndex() == 1) // 1 means MRI
  {
    m_LaplacianFilter2->Update();
    mitk::CastToMitkImage(m_LaplacianFilter2->GetOutput(), m_ImageToRegister);

    if (m_ImageToRegister->GetScalarValueMax() < 150)
    {
      m_BinaryThresholdFilter->SetLowerThreshold(2);
    }
    else if (m_ImageToRegister->GetScalarValueMax() < 300)
    {
      m_BinaryThresholdFilter->SetLowerThreshold(10);
    }
    else if (m_ImageToRegister->GetScalarValueMax() < 400)
    {
      m_BinaryThresholdFilter->SetLowerThreshold(15);
    }
    else
    {
      m_BinaryThresholdFilter->SetLowerThreshold(20);
    }
    ImageType::Pointer itkImage2 = ImageType::New();
    mitk::CastToItkImage(m_ImageToRegister, itkImage2);

    m_BinaryThresholdFilter->SetInput(itkImage2);
  }
  else
  {
    m_BinaryThresholdFilter->SetInput(m_LaplacianFilter2->GetOutput());
  }

  m_HoleFillingFilter->SetInput(m_BinaryThresholdFilter->GetOutput());
  m_BinaryImageToShapeLabelMapFilter->SetInput(m_HoleFillingFilter->GetOutput());
  m_BinaryImageToShapeLabelMapFilter->Update();

  ImageType::Pointer binaryImage = ImageType::New();
  binaryImage = m_HoleFillingFilter->GetOutput();

  this->EliminateTooSmallLabeledObjects(binaryImage);
  mitk::CastToMitkImage(binaryImage, m_ImageToRegister);
  // mitk::CastToMitkImage(itkImageOld, m_ImageToRegister);
  return true;
}

void QmitkAutomaticFiducialmarkerRegistrationWidget::InitializeImageFilters()
{
  // Initialize threshold filters
  m_ThresholdFilter = itk::ThresholdImageFilter<ImageType>::New();
  m_ThresholdFilter->SetOutsideValue(0);

  if (m_Controls->selectedModalityComboBox->currentIndex() == 1) // 1 means MRI
  {
    m_ThresholdFilter->SetLower(floor(m_ImageToRegister->GetScalarValueMax() * 0.05));
    m_ThresholdFilter->SetUpper(1800);
  }
  else
  {
    m_ThresholdFilter->SetLower(500);
    if (m_Controls->selectedModalityComboBox->currentIndex() == 2) // 2 means 3D Fluoroskopie:
    {
      m_ThresholdFilter->SetUpper(40000);
    }
    else // Fall CT:
    {
      m_ThresholdFilter->SetUpper(3200);
    }
  }

  // Initialize binary threshold filter 1
  m_BinaryThresholdFilter = BinaryThresholdImageFilterType::New();
  m_BinaryThresholdFilter->SetOutsideValue(0);
  m_BinaryThresholdFilter->SetInsideValue(1);
  if (m_Controls->selectedModalityComboBox->currentIndex() == 1) // 1 means MRI
  {
    // Lower threshold is configured dynamically during image processing
    m_BinaryThresholdFilter->SetUpperThreshold(500);
  }
  else
  {
    m_BinaryThresholdFilter->SetLowerThreshold(350);
    m_BinaryThresholdFilter->SetUpperThreshold(10000);
  }

  // Initialize laplacian recursive gaussian image filter
  m_LaplacianFilter1 = LaplacianRecursiveGaussianImageFilterType::New();
  m_LaplacianFilter2 = LaplacianRecursiveGaussianImageFilterType::New();

  // Initialize binary hole filling filter
  m_HoleFillingFilter = VotingBinaryIterativeHoleFillingImageFilterType::New();
  VotingBinaryIterativeHoleFillingImageFilterType::InputSizeType radius;
  radius.Fill(1);
  m_HoleFillingFilter->SetRadius(radius);
  m_HoleFillingFilter->SetBackgroundValue(0);
  m_HoleFillingFilter->SetForegroundValue(1);
  m_HoleFillingFilter->SetMaximumNumberOfIterations(5);

  // Initialize binary image to shape label map filter
  m_BinaryImageToShapeLabelMapFilter = BinaryImageToShapeLabelMapFilterType::New();
  m_BinaryImageToShapeLabelMapFilter->SetInputForegroundValue(1);
}

double QmitkAutomaticFiducialmarkerRegistrationWidget::GetCharacteristicDistanceAWithUpperMargin()
{
  switch (m_Controls->fiducialMarkerConfigurationComboBox->currentIndex())
  {
      // case 0 is equal to fiducial marker configuration 15mm distance
    case 0:
      return 18.105;

      // case 1 is equal to fiducial marker configuration 20mm distance
    case 1:
      return 24.14;

      // case 2 is equal to fiducial marker configuration I of master project
    case 2:
      return 20.20;
  }
  return 0.0;
}

double QmitkAutomaticFiducialmarkerRegistrationWidget::GetCharacteristicDistanceBWithLowerMargin()
{
  switch (m_Controls->fiducialMarkerConfigurationComboBox->currentIndex())
  {
      // case 0 is equal to fiducial marker configuration 15mm distance
    case 0:
      return 18.105;

      // case 1 is equal to fiducial marker configuration 20mm distance
    case 1:
      return 24.14;

      // case 2 is equal to fiducial marker configuration I of master project
    case 2:
      return 20.20;
  }
  return 0.0;
}

double QmitkAutomaticFiducialmarkerRegistrationWidget::GetCharacteristicDistanceBWithUpperMargin()
{
  switch (m_Controls->fiducialMarkerConfigurationComboBox->currentIndex())
  {
      // case 0 is equal to fiducial marker configuration 15mm distance
    case 0:
      return 23.595;

      // case 1 is equal to fiducial marker configuration 20mm distance
    case 1:
      return 31.46;

      // case 2 is equal to fiducial marker configuration I of master project
    case 2:
      return 28.01;
  }
  return 0.0;
}

double QmitkAutomaticFiducialmarkerRegistrationWidget::GetMinimalFiducialConfigurationDistance()
{
  switch (m_Controls->fiducialMarkerConfigurationComboBox->currentIndex())
  {
      // case 0 is equal to fiducial marker configuration 15mm distance
    case 0:
      return 15.0;

      // case 1 is equal to fiducial marker configuration 20mm distance
    case 1:
      return 20.0;

      // case 2 is equal to fiducial marker configuration I of master project
    case 2:
      return 15.10;
  }
  return 0.0;
}

void QmitkAutomaticFiducialmarkerRegistrationWidget::CreateMarkerModelCoordinateSystemPointSet()
{
  if (m_MarkerModelCoordinateSystemPointSet.IsNull())
  {
    m_MarkerModelCoordinateSystemPointSet = mitk::PointSet::New();
  }
  else
  {
    m_MarkerModelCoordinateSystemPointSet->Clear();
  }

  mitk::Point3D fiducial1;
  mitk::Point3D fiducial2;
  mitk::Point3D fiducial3;
  mitk::Point3D fiducial4;
  mitk::Point3D fiducial5;
  mitk::Point3D fiducial6;
  mitk::Point3D fiducial7;
  mitk::Point3D fiducial8;

  switch (m_Controls->fiducialMarkerConfigurationComboBox->currentIndex())
  {
      // case 0 is equal to fiducial marker configuration 15mm distance
    case 0:
      fiducial1[0] = 0;
      fiducial1[1] = 0;
      fiducial1[2] = 0;

      fiducial2[0] = 0;
      fiducial2[1] = 15;
      fiducial2[2] = 0;

      fiducial3[0] = 15;
      fiducial3[1] = 0;
      fiducial3[2] = 0;

      fiducial4[0] = 30;
      fiducial4[1] = 30;
      fiducial4[2] = 0;

      fiducial5[0] = 0;
      fiducial5[1] = 30;
      fiducial5[2] = 15;

      fiducial6[0] = 15;
      fiducial6[1] = 30;
      fiducial6[2] = 15;

      fiducial7[0] = 30;
      fiducial7[1] = 15;
      fiducial7[2] = 15;

      fiducial8[0] = 30;
      fiducial8[1] = 0;
      fiducial8[2] = 15;
      break;
      // case 1 is equal to fiducial marker configuration 20mm distance
    case 1:
      fiducial1[0] = 0;
      fiducial1[1] = 0;
      fiducial1[2] = 0;

      fiducial2[0] = 0;
      fiducial2[1] = 20;
      fiducial2[2] = 0;

      fiducial3[0] = 20;
      fiducial3[1] = 0;
      fiducial3[2] = 0;

      fiducial4[0] = 40;
      fiducial4[1] = 40;
      fiducial4[2] = 0;

      fiducial5[0] = 0;
      fiducial5[1] = 40;
      fiducial5[2] = 20;

      fiducial6[0] = 20;
      fiducial6[1] = 40;
      fiducial6[2] = 20;

      fiducial7[0] = 40;
      fiducial7[1] = 20;
      fiducial7[2] = 20;

      fiducial8[0] = 40;
      fiducial8[1] = 0;
      fiducial8[2] = 20;
      break;

    case 2:
      MITK_WARN("Is not implemented yet!!! --> TODO!");
      break;
  }

  m_MarkerModelCoordinateSystemPointSet->InsertPoint(0, fiducial1);
  m_MarkerModelCoordinateSystemPointSet->InsertPoint(1, fiducial2);
  m_MarkerModelCoordinateSystemPointSet->InsertPoint(2, fiducial3);
  m_MarkerModelCoordinateSystemPointSet->InsertPoint(3, fiducial4);
  m_MarkerModelCoordinateSystemPointSet->InsertPoint(4, fiducial5);
  m_MarkerModelCoordinateSystemPointSet->InsertPoint(5, fiducial6);
  m_MarkerModelCoordinateSystemPointSet->InsertPoint(6, fiducial7);
  m_MarkerModelCoordinateSystemPointSet->InsertPoint(7, fiducial8);
}

void QmitkAutomaticFiducialmarkerRegistrationWidget::EliminateTooSmallLabeledObjects(ImageType::Pointer binaryImage)
{
  BinaryImageToShapeLabelMapFilterType::OutputImageType::Pointer labelMap =
    m_BinaryImageToShapeLabelMapFilter->GetOutput();
  double voxelVolume = this->GetVoxelVolume();
  double fiducialVolume;
  unsigned int numberOfPixels;

  if (m_Controls->selectedModalityComboBox->currentIndex() == 1) // 1 means MRI
  {
    fiducialVolume = this->GetFiducialVolume(3.0);
    numberOfPixels = ceil(fiducialVolume / voxelVolume);
  }
  else if (m_Controls->fiducialDiameter3mmRadioButton->isChecked())
  {
    fiducialVolume = this->GetFiducialVolume(1.5);
    numberOfPixels = ceil(fiducialVolume / voxelVolume);
  }
  else
  {
    fiducialVolume = this->GetFiducialVolume(2.5);
    numberOfPixels = ceil(fiducialVolume / voxelVolume);
  }

  MITK_INFO << "Voxel Volume = " << voxelVolume << "; Fiducial Volume = " << fiducialVolume;
  MITK_INFO << "Number of pixels = " << numberOfPixels;

  labelMap = m_BinaryImageToShapeLabelMapFilter->GetOutput();
  // The output of this filter is an itk::LabelMap, which contains itk::LabelObject's
  MITK_INFO << "There are " << labelMap->GetNumberOfLabelObjects() << " objects.";

  // Loop over each region
  for (int i = labelMap->GetNumberOfLabelObjects() - 1; i >= 0; --i)
  {
    // Get the ith region
    BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType *labelObject =
      labelMap->GetNthLabelObject(i);

    if (labelObject->Size() < numberOfPixels * 0.6 ||
        (labelObject->Size() > numberOfPixels * 2 && m_Controls->selectedModalityComboBox->currentIndex() == 1))
    {
      /*for (unsigned int pixelId = 0; pixelId < labelObject->Size(); pixelId++)
      {
      binaryImage->SetPixel(labelObject->GetIndex(pixelId), 0);
      }*/
      labelMap->RemoveLabelObject(labelObject);
    }
  }
}

bool QmitkAutomaticFiducialmarkerRegistrationWidget::EliminateFiducialCandidatesByEuclideanDistances()
{
  if ((m_FiducialCandidates.size() < 8 && m_Controls->fiducialMarkerConfigurationComboBox->currentIndex() <= 1) ||
      (m_FiducialCandidates.size() < 4 && m_Controls->fiducialMarkerConfigurationComboBox->currentIndex() == 2))
  {
    return false;
  }

  for (unsigned int counter = 0; counter < m_FiducialCandidates.size(); ++counter)
  {
    int amountOfAcceptedFiducials = 0;
    mitk::Point3D fiducialCentroid(m_FiducialCandidates.at(counter).first);
    // Loop through all fiducial candidates and calculate the distance between the chosen fiducial
    // candidate and the other candidates. For each candidate with a right distance between
    // Configuration A: 7.93mm and 31.0mm (10 mm distance between fiducial centers) or
    // Configuration B: 11.895mm and 45.0mm (15 mm distance between fiducial centers) or
    // Configuration C: 15.86mm and 60.0mm (20 mm distance between fiducial centers)
    //
    // increase the amountOfAcceptedFiducials.
    for (unsigned int position = 0; position < m_FiducialCandidates.size(); ++position)
    {
      if (position == counter)
      {
        continue;
      }
      mitk::Point3D otherCentroid(m_FiducialCandidates.at(position).first);
      double distance = fiducialCentroid.EuclideanDistanceTo(otherCentroid);

      switch (m_Controls->fiducialMarkerConfigurationComboBox->currentIndex())
      {
          // case 0 is equal to fiducial marker configuration 15mm distance
        case 0:
          if (distance > 11.895 && distance < 45.0)
          {
            ++amountOfAcceptedFiducials;
          }
          break;
          // case 1 is equal to fiducial marker configuration 20mm distance
        case 1:
          if (distance > 15.86 && distance < 60.0)
          {
            ++amountOfAcceptedFiducials;
          }
          break;
          // case 2 is equal to fiducial marker configuration I of master project
        case 2:
        {
          if (distance > 11.895 && distance < 53.0)
          {
            ++amountOfAcceptedFiducials;
          }
        }
      }
    }
    // The amountOfAcceptedFiducials must be at least 7 for cases 0-1 and at least 3 for case 2. Otherwise delete the
    // fiducial candidate from the list of candidates.
    if ((amountOfAcceptedFiducials < 7 && m_Controls->fiducialMarkerConfigurationComboBox->currentIndex() <= 1) ||
        (amountOfAcceptedFiducials < 3 && m_Controls->fiducialMarkerConfigurationComboBox->currentIndex() == 2))
    {
      MITK_INFO << "Deleting fiducial candidate at position: " << m_FiducialCandidates.at(counter).first;
      m_FiducialCandidates.erase(m_FiducialCandidates.begin() + counter);
      if ((m_FiducialCandidates.size() < 8 && m_Controls->fiducialMarkerConfigurationComboBox->currentIndex() <= 1) ||
          (m_FiducialCandidates.size() < 4 && m_Controls->fiducialMarkerConfigurationComboBox->currentIndex() == 2))
      {
        return false;
      }
      counter = -1;
    }
  }

  if ((m_FiducialCandidates.size() > 8 && m_Controls->fiducialMarkerConfigurationComboBox->currentIndex() <= 1) ||
      (m_FiducialCandidates.size() > 4 && m_Controls->fiducialMarkerConfigurationComboBox->currentIndex() == 2))
  {
    // Classify the rested fiducial candidates by its characteristic Euclidean distances
    // between the canidates and remove all candidates with a false distance configuration:
    this->ClassifyFiducialCandidates();
  }

  return true;
}

void QmitkAutomaticFiducialmarkerRegistrationWidget::EliminateNearFiducialCandidatesByMaxDistanceToCentroids()
{
  std::vector<std::vector<double>> distanceVectorsFiducials;
  this->CalculateDistancesBetweenFiducials(distanceVectorsFiducials);

  double limit = 10.0;

  std::vector<unsigned int> indices_1;
  std::vector<unsigned int> indices_2;
  // Find all pairs of fiducial candidates, whose euclidean distance is smaller than
  // the above defined limit value. As one fiducial candidate of each pair does not
  // correspond to a true fiducial, these have to be determined in the following
  for (unsigned int i = 0; i < distanceVectorsFiducials.size(); ++i)
  {
    if (distanceVectorsFiducials.at(i).at(0) < 10.0)
    {
      for (unsigned int index_rest_of_vector = i + 1; index_rest_of_vector < distanceVectorsFiducials.size();
           ++index_rest_of_vector)
      {
        if (distanceVectorsFiducials.at(i).at(0) == distanceVectorsFiducials.at(index_rest_of_vector).at(0))
        {
          indices_1.push_back(i);
          indices_2.push_back(index_rest_of_vector);
          continue;
        }
      }
    }
  }

  // Determine the maximum distance to the centroid of the fiducial candidate
  // The true fiducial corresponds to the fiducial candidate with a smaller
  // maximum distance to the centroid.
  std::vector<unsigned int> candidateIndicesToRemove;
  for (unsigned int count = 0; count < indices_1.size(); ++count)
  {
    double maxDistance_1 = this->GetMaxDistanceToCentroidOfFiducialCandidate(indices_1.at(count));
    double maxDistance_2 = this->GetMaxDistanceToCentroidOfFiducialCandidate(indices_2.at(count));
    MITK_INFO << "maxDinstance 1 = " << maxDistance_1 << " | maxDistance 2 = " << maxDistance_2;
    if (maxDistance_1 > maxDistance_2)
    {
      candidateIndicesToRemove.push_back(indices_1.at(count));
    }
    else
    {
      candidateIndicesToRemove.push_back(indices_2.at(count));
    }
  }
  // Sort the distances from low to high indices
  std::sort(candidateIndicesToRemove.begin(), candidateIndicesToRemove.end());

  // Remove false fiducial candidates
  for (int count = candidateIndicesToRemove.size() - 1; count >= 0; --count)
  {
    MITK_INFO << "Removing fiducial candidate at index " << candidateIndicesToRemove.at(count) << " with centroid "
              << m_FiducialCandidates.at(candidateIndicesToRemove.at(count)).first;
    m_FiducialCandidates.erase(m_FiducialCandidates.begin() + candidateIndicesToRemove.at(count));
  }
}

void QmitkAutomaticFiducialmarkerRegistrationWidget::ClassifyFiducialCandidates()
{
  MITK_INFO << "ClassifyFiducialCandidates()";
  std::vector<int> fiducialCandidatesToBeRemoved;
  std::vector<std::vector<double>> distanceVectorsFiducials;
  this->CalculateDistancesBetweenFiducials(distanceVectorsFiducials);

  for (unsigned int counter = 0; counter < distanceVectorsFiducials.size(); ++counter)
  {
    int distanceA = 0; // => 10,00mm distance
    int distanceB = 0; // => 14,14mm distance
    int distanceC = 0; // => 17,32mm distance
    int distanceD = 0; // => 22,36mm distance
    int distanceE = 0; // => 24,49mm distance
    int distanceF = 0; // => 28,28mm distance

    std::vector<double> &distances = distanceVectorsFiducials.at(counter);
    for (unsigned int number = 0; number < distances.size(); ++number)
    {
      double &distance = distances.at(number);
      switch (m_Controls->fiducialMarkerConfigurationComboBox->currentIndex())
      {
          // case 0 is equal to fiducial marker configuration 15mm distance
        case 0:
          if (distance > 11.895 && distance <= 18.105)
          {
            ++distanceA;
          }
          else if (distance > 18.105 && distance <= 23.595)
          {
            ++distanceB;
          }
          else if (distance > 23.595 && distance <= 29.76)
          {
            ++distanceC;
          }
          else if (distance > 29.76 && distance <= 35.1375)
          {
            ++distanceD;
            if (distance > 33.54)
            {
              ++distanceE;
            }
          }
          else if (distance > 35.1375 && distance <= 39.5775)
          {
            ++distanceE;
            if (distance < 36.735)
            {
              ++distanceD;
            }
          }
          else if (distance > 39.5775 && distance <= 45.00)
          {
            ++distanceF;
          }
          break;

          // case 1 is equal to fiducial marker configuration 20mm distance
        case 1:
          if (distance > 15.86 && distance <= 24.14)
          {
            ++distanceA;
          }
          else if (distance > 24.14 && distance <= 31.46)
          {
            ++distanceB;
          }
          else if (distance > 31.46 && distance <= 39.68)
          {
            ++distanceC;
          }
          else if (distance > 39.68 && distance <= 46.85)
          {
            ++distanceD;
            if (distance > 44.00)
            {
              ++distanceE;
            }
          }
          else if (distance > 46.85 && distance <= 52.77)
          {
            ++distanceE;
            if (distance < 49.00)
            {
              ++distanceD;
            }
          }
          else if (distance > 52.77 && distance <= 60.00)
          {
            ++distanceF;
          }
          break;

          // case 2 is equal to fiducial marker configuration I of master project
        case 2:
          if (distance > 11.895 && distance <= 20.20)
          {
            ++distanceA;
          }
          else if (distance > 20.20 && distance <= 28.01)
          {
            ++distanceB;
            if (distance > 24.00)
            {
              ++distanceC;
              ++distanceD;
            }
          }
          else if (distance > 28.01 && distance <= 33.91)
          {
            ++distanceC;
            ++distanceD;
            if (distance < 30.00)
              ++distanceB;
            else if (distance >= 30.00)
              ++distanceE;
          }
          else if (distance > 33.91 && distance <= 41.53)
          {
            ++distanceE;
          }
          else if (distance > 41.53 && distance <= 53.00)
          {
            ++distanceF;
          }
          break;
      }
    } // End for-loop distances-vector

    // Now, having looped through all distances of one fiducial candidate, check
    // if the combination of different distances is known. The >= is due to the
    // possible occurrence of other fiducial candidates that have an distance equal to
    // one of the distances A - E. However, false fiducial candidates outside
    // the fiducial marker does not have the right distance configuration:

    switch (m_Controls->fiducialMarkerConfigurationComboBox->currentIndex())
    {
        // case 0 is equal to fiducial marker configuration 15mm distance
      case 0:
        // case 1 is equal to fiducial marker configuration 20mm distance
      case 1:
        if (((distanceA >= 2 && distanceD >= 2 && distanceE >= 2 && distanceF >= 1) ||
             (distanceA >= 1 && distanceB >= 2 && distanceC >= 1 && distanceD >= 2 && distanceE >= 1) ||
             (distanceB >= 2 && distanceD >= 4 && distanceF >= 1) ||
             (distanceA >= 1 && distanceB >= 1 && distanceD >= 3 && distanceE >= 1 && distanceF >= 1)) == false)
        {
          MITK_INFO << "Detected fiducial candidate with unknown distance configuration.";
          fiducialCandidatesToBeRemoved.push_back(counter);
        }
        break;

        // case 2 is equal to fiducial marker configuration I of master project
      case 2:
        if (((distanceB >= 1 && distanceC >= 1 && distanceD >= 1) ||
             (distanceB >= 1 && distanceE >= 1 && distanceF >= 1) ||
             (distanceA >= 1 && distanceD >= 1 && distanceE >= 1) ||
             (distanceA >= 1 && distanceC >= 1 && distanceF >= 1)) == false)
        {
          MITK_INFO << "Detected fiducial candidate with unknown distance configuration.";
          fiducialCandidatesToBeRemoved.push_back(counter);
        }
        break;
    }
  }
  for (int count = fiducialCandidatesToBeRemoved.size() - 1; count >= 0; --count)
  {
    MITK_INFO << "Removing fiducial candidate " << fiducialCandidatesToBeRemoved.at(count);
    m_FiducialCandidates.erase(m_FiducialCandidates.begin() + fiducialCandidatesToBeRemoved.at(count));
  }
}

void QmitkAutomaticFiducialmarkerRegistrationWidget::GetCentroidsOfLabeledObjects()
{
  MITK_INFO << "GetCentroidsOfLabeledObjects()";
  BinaryImageToShapeLabelMapFilterType::OutputImageType::Pointer labelMap =
    m_BinaryImageToShapeLabelMapFilter->GetOutput();
  for (int i = labelMap->GetNumberOfLabelObjects() - 1; i >= 0; --i)
  {
    // Get the ith region
    ShapeLabelObjectType *labelObject = labelMap->GetNthLabelObject(i);

    mitk::Vector3D centroid;
    centroid[0] = labelObject->GetCentroid()[0];
    centroid[1] = labelObject->GetCentroid()[1];
    centroid[2] = labelObject->GetCentroid()[2];

    m_FiducialCandidates.push_back(FiducialCandidatePairType(centroid, labelObject));

    MITK_INFO << "Position: " << centroid << " | object " << i << " contains " << labelObject->Size() << " pixel";
  }
  // evtl. for later: itk::LabelMapOverlayImageFilter
}

void QmitkAutomaticFiducialmarkerRegistrationWidget::NumerateFiducialMarks()
{
  MITK_INFO << "NumerateFiducialMarks()";
  bool successFiducialNo1;
  bool successFiducialNo4;
  bool successFiducialNo2;
  bool successFiducialNo3;
  bool successFiducialNo2And3;
  bool successFiducialNo5;
  bool successFiducialNo8;
  bool successFiducialNo6;
  bool successFiducialNo7;

  std::vector<std::vector<double>> distanceVectorsFiducials;
  this->CalculateDistancesBetweenFiducials(distanceVectorsFiducials);

  if (m_Controls->fiducialMarkerConfigurationComboBox->currentIndex() <= 1)
  {
    successFiducialNo1 = this->FindFiducialNo1(distanceVectorsFiducials);
    successFiducialNo4 = this->FindFiducialNo4(distanceVectorsFiducials);
    successFiducialNo2And3 = this->FindFiducialNo2And3();
    successFiducialNo5 = this->FindFiducialNo5();
    successFiducialNo8 = this->FindFiducialNo8();
    successFiducialNo6 = this->FindFiducialNo6();
    successFiducialNo7 = this->FindFiducialNo7();

    if (!successFiducialNo1 || !successFiducialNo4 || !successFiducialNo2And3 || !successFiducialNo5 ||
        !successFiducialNo8 || !successFiducialNo6 || !successFiducialNo7)
    {
      QMessageBox msgBox;
      msgBox.setText("Cannot numerate/localize all fiducials successfully.");
      msgBox.exec();
      return;
    }
  }
  else if (m_Controls->fiducialMarkerConfigurationComboBox->currentIndex() == 2)
  {
    double characteristicDistanceAWithUpperMargin = this->GetCharacteristicDistanceAWithUpperMargin();

    // Find Fiducial 4:
    for (unsigned int i = 0; i < distanceVectorsFiducials.size(); ++i)
    {
      std::vector<double> &distances = distanceVectorsFiducials.at(i);

      if (distances.at(0) <= characteristicDistanceAWithUpperMargin &&
          distances.at(2) > 41.53) //  -> characteristic Distance F with lower margin
      {
        MITK_INFO << "Found Fiducial 4 (PointSet number " << i << ")";
        m_FiducialMarkerCentroids.insert(std::pair<int, mitk::Point3D>(4, m_FiducialCandidates.at(i).first));
        distanceVectorsFiducials.erase(distanceVectorsFiducials.begin() + i);
        m_FiducialCandidates.erase(m_FiducialCandidates.begin() + i);
        successFiducialNo4 = true;
        break;
      }
    }

    // Find Fiducial 3:
    for (unsigned int i = 0; i < distanceVectorsFiducials.size(); ++i)
    {
      std::vector<double> &distances = distanceVectorsFiducials.at(i);

      if (distances.at(0) <= characteristicDistanceAWithUpperMargin)
      {
        MITK_INFO << "Found Fiducial 3 (PointSet number " << i << ")";
        m_FiducialMarkerCentroids.insert(std::pair<int, mitk::Point3D>(3, m_FiducialCandidates.at(i).first));
        distanceVectorsFiducials.erase(distanceVectorsFiducials.begin() + i);
        m_FiducialCandidates.erase(m_FiducialCandidates.begin() + i);
        successFiducialNo3 = true;
        break;
      }
    }

    // Find Fiducial 2:
    for (unsigned int i = 0; i < distanceVectorsFiducials.size(); ++i)
    {
      std::vector<double> &distances = distanceVectorsFiducials.at(i);

      if (distances.at(2) > 41.53) //  -> characteristic Distance F with lower margin
      {
        MITK_INFO << "Found Fiducial 2 (PointSet number " << i << ")";
        m_FiducialMarkerCentroids.insert(std::pair<int, mitk::Point3D>(2, m_FiducialCandidates.at(i).first));
        distanceVectorsFiducials.erase(distanceVectorsFiducials.begin() + i);
        m_FiducialCandidates.erase(m_FiducialCandidates.begin() + i);
        successFiducialNo2 = true;
        break;
      }
    }

    // Find Fiducial 1:
    double characteristicDistanceBWithLowerMargin = this->GetCharacteristicDistanceBWithLowerMargin();

    for (unsigned int i = 0; i < distanceVectorsFiducials.size(); ++i)
    {
      std::vector<double> &distances = distanceVectorsFiducials.at(i);

      if (distances.at(0) > characteristicDistanceBWithLowerMargin && distances.at(0) <= 30.0 && distances.at(2) < 37.0)
      {
        MITK_INFO << "Found Fiducial 1 (PointSet number " << i << ")";
        m_FiducialMarkerCentroids.insert(std::pair<int, mitk::Point3D>(1, m_FiducialCandidates.at(i).first));
        distanceVectorsFiducials.erase(distanceVectorsFiducials.begin() + i);
        m_FiducialCandidates.erase(m_FiducialCandidates.begin() + i);
        successFiducialNo1 = true;
        break;
      }
    }

    if (!successFiducialNo4 || !successFiducialNo3 || !successFiducialNo2 || !successFiducialNo1)
    {
      QMessageBox msgBox;
      msgBox.setText("Cannot numerate/localize all fiducials successfully.");
      msgBox.exec();
      return;
    }
  }

  if (m_MarkerImageToRegisterCoordinateSystemPointSet.IsNull())
  {
    m_MarkerImageToRegisterCoordinateSystemPointSet = mitk::PointSet::New();
  }
  else if (m_MarkerImageToRegisterCoordinateSystemPointSet->GetSize() != 0)
  {
    m_MarkerImageToRegisterCoordinateSystemPointSet->Clear();
  }
  MITK_INFO << m_FiducialMarkerCentroids.size();
  for (unsigned int counter = 1; counter <= m_FiducialMarkerCentroids.size(); ++counter)
  {
    m_MarkerImageToRegisterCoordinateSystemPointSet->InsertPoint(counter - 1, m_FiducialMarkerCentroids.at(counter));
  }

  mitk::DataNode::Pointer node = mitk::DataNode::New();
  node->SetData(m_MarkerImageToRegisterCoordinateSystemPointSet);
  node->SetName("MarkerFloatingImageCSPointSet");
  node->SetFloatProperty("pointsize", 5.0);
  this->m_DataStorage->Add(node);
}

void QmitkAutomaticFiducialmarkerRegistrationWidget::CalculateDistancesBetweenFiducials(
  std::vector<std::vector<double>> &distanceVectorsFiducials)
{
  std::vector<double> distancesBetweenFiducials;

  for (unsigned int i = 0; i < m_FiducialCandidates.size(); ++i)
  {
    distancesBetweenFiducials.clear();
    mitk::Point3D fiducialCentroid(m_FiducialCandidates.at(i).first);
    for (unsigned int n = 0; n < m_FiducialCandidates.size(); ++n)
    {
      mitk::Point3D otherCentroid(m_FiducialCandidates.at(n).first);
      distancesBetweenFiducials.push_back(fiducialCentroid.EuclideanDistanceTo(otherCentroid));
    }
    // Sort the distances from low to big numbers
    std::sort(distancesBetweenFiducials.begin(), distancesBetweenFiducials.end());
    // First entry of the distance vector must be 0, so erase it
    if (distancesBetweenFiducials.at(0) == 0.0)
    {
      distancesBetweenFiducials.erase(distancesBetweenFiducials.begin());
    }
    // Add the distance vector to the collecting distances vector
    distanceVectorsFiducials.push_back(distancesBetweenFiducials);
  }

  // Uncomment this block of code for development purposes:
  for (unsigned int i = 0; i < distanceVectorsFiducials.size(); ++i)
  {
    MITK_INFO << "Vector " << i << ":";
    for (unsigned int k = 0; k < distanceVectorsFiducials.at(i).size(); ++k)
    {
      MITK_INFO << distanceVectorsFiducials.at(i).at(k);
    }
  }
}

double QmitkAutomaticFiducialmarkerRegistrationWidget::GetMaxDistanceToCentroidOfFiducialCandidate(unsigned int &index)
{
  mitk::Point3D centroid(m_FiducialCandidates.at(index).first);
  mitk::Point3D pixelPoint;
  ShapeLabelObjectType *labelObject = m_FiducialCandidates.at(index).second;
  double maxDistance = 0.0;
  double distance = 0.0;

  for (unsigned int pixelId = 0; pixelId < labelObject->Size(); pixelId++)
  {
    m_ImageToRegister->GetGeometry()->IndexToWorld(labelObject->GetIndex(pixelId), pixelPoint);
    distance = centroid.EuclideanDistanceTo(pixelPoint);
    maxDistance = maxDistance < distance ? distance : maxDistance;
  }

  return maxDistance;
}

bool QmitkAutomaticFiducialmarkerRegistrationWidget::FindFiducialNo1(
  std::vector<std::vector<double>> &distanceVectorsFiducials)
{
  for (unsigned int i = 0; i < distanceVectorsFiducials.size(); ++i)
  {
    std::vector<double> &distances = distanceVectorsFiducials.at(i);
    if ((distances.size() < 7 && m_Controls->fiducialMarkerConfigurationComboBox->currentIndex() <= 1) ||
        (distances.size() < 3 && m_Controls->fiducialMarkerConfigurationComboBox->currentIndex() == 2))
    {
      MITK_WARN << "Cannot find fiducial 1, there aren't found enough fiducial candidates.";
      return false;
    }
    double characteristicDistanceAWithUpperMargin = this->GetCharacteristicDistanceAWithUpperMargin();

    if (distances.at(0) <= characteristicDistanceAWithUpperMargin &&
        distances.at(1) <= characteristicDistanceAWithUpperMargin)
    {
      MITK_INFO << "Found Fiducial 1 (PointSet number " << i << ")";
      m_FiducialMarkerCentroids.insert(std::pair<int, mitk::Point3D>(1, m_FiducialCandidates.at(i).first));
      distanceVectorsFiducials.erase(distanceVectorsFiducials.begin() + i);
      m_FiducialCandidates.erase(m_FiducialCandidates.begin() + i);
      return true;
    }
  }
  return false;
}

bool QmitkAutomaticFiducialmarkerRegistrationWidget::FindFiducialNo2And3()
{
  if (m_FiducialMarkerCentroids.find(1) == m_FiducialMarkerCentroids.end())
  {
    MITK_WARN << "Cannot find fiducial No 2 and 3. Before must be found fiducial No 1.";
    return false;
  }

  mitk::Point3D fiducialNo1(m_FiducialMarkerCentroids.at(1));
  mitk::Point3D fiducialVectorA;
  mitk::Point3D fiducialVectorB;
  mitk::Point3D fiducialPointA;
  mitk::Point3D fiducialPointB;
  bool foundFiducialA = false;
  bool foundFiducialB = false;
  mitk::Vector3D vectorFiducial1ToFiducialA;
  mitk::Vector3D vectorFiducial1ToFiducialB;

  for (unsigned int i = 0; i < m_FiducialCandidates.size(); ++i)
  {
    mitk::Point3D fiducialCentroid(m_FiducialCandidates.at(i).first);
    double distance = fiducialNo1.EuclideanDistanceTo(fiducialCentroid);
    if (distance <= this->GetCharacteristicDistanceAWithUpperMargin())
    {
      fiducialVectorA = m_FiducialCandidates.at(i).first;
      fiducialPointA = fiducialCentroid;
      m_FiducialCandidates.erase(m_FiducialCandidates.begin() + i);
      foundFiducialA = true;
      break;
    }
  }

  for (unsigned int i = 0; i < m_FiducialCandidates.size(); ++i)
  {
    mitk::Point3D fiducialCentroid(m_FiducialCandidates.at(i).first);
    double distance = fiducialNo1.EuclideanDistanceTo(fiducialCentroid);
    if (distance <= this->GetCharacteristicDistanceAWithUpperMargin())
    {
      fiducialVectorB = m_FiducialCandidates.at(i).first;
      fiducialPointB = fiducialCentroid;
      m_FiducialCandidates.erase(m_FiducialCandidates.begin() + i);
      foundFiducialB = true;
      break;
    }
  }

  if (!foundFiducialA || !foundFiducialB)
  {
    MITK_WARN << "Cannot identify fiducial candidates 2 and 3";
    return false;
  }
  else if (m_FiducialCandidates.size() == 0)
  {
    MITK_WARN << "Too less fiducials detected. Cannot identify fiducial candidates 2 and 3";
    return false;
  }

  vectorFiducial1ToFiducialA = fiducialVectorA - m_FiducialMarkerCentroids.at(1);
  vectorFiducial1ToFiducialB = fiducialVectorB - m_FiducialMarkerCentroids.at(1);

  vnl_vector<double> crossProductVnl =
    vnl_cross_3d(vectorFiducial1ToFiducialA.GetVnlVector(), vectorFiducial1ToFiducialB.GetVnlVector());
  mitk::Vector3D crossProduct;
  crossProduct.SetVnlVector(crossProductVnl);

  mitk::Vector3D vectorFiducial1ToRandomLeftFiducial =
    m_FiducialCandidates.at(0).first - m_FiducialMarkerCentroids.at(1);

  double scalarProduct = (crossProduct * vectorFiducial1ToRandomLeftFiducial) /
                         (crossProduct.GetNorm() * vectorFiducial1ToRandomLeftFiducial.GetNorm());

  double alpha = acos(scalarProduct) * 57.29578; // Transform into degree
  MITK_INFO << "Scalar Product = " << alpha;

  if (alpha <= 90)
  {
    m_FiducialMarkerCentroids[3] = fiducialVectorA;
    m_FiducialMarkerCentroids[2] = fiducialVectorB;
  }
  else
  {
    m_FiducialMarkerCentroids[2] = fiducialVectorA;
    m_FiducialMarkerCentroids[3] = fiducialVectorB;
  }

  MITK_INFO << "Found Fiducial 2, PointSet: " << m_FiducialMarkerCentroids.at(2);
  MITK_INFO << "Found Fiducial 3, PointSet: " << m_FiducialMarkerCentroids.at(3);

  return true;
}

bool QmitkAutomaticFiducialmarkerRegistrationWidget::FindFiducialNo4(
  std::vector<std::vector<double>> &distanceVectorsFiducials)
{
  double characteristicDistanceBWithLowerMargin = this->GetCharacteristicDistanceBWithLowerMargin();
  double characteristicDistanceBWithUpperMargin = this->GetCharacteristicDistanceBWithUpperMargin();

  for (unsigned int i = 0; i < distanceVectorsFiducials.size(); ++i)
  {
    std::vector<double> &distances = distanceVectorsFiducials.at(i);
    if ((distances.size() < 7 && m_Controls->fiducialMarkerConfigurationComboBox->currentIndex() <= 1) ||
        (distances.size() < 3 && m_Controls->fiducialMarkerConfigurationComboBox->currentIndex() == 2))
    {
      MITK_WARN << "Cannot find fiducial 4, there aren't found enough fiducial candidates.";
      return false;
    }

    if (distances.at(0) > characteristicDistanceBWithLowerMargin &&
        distances.at(0) <= characteristicDistanceBWithUpperMargin &&
        distances.at(1) > characteristicDistanceBWithLowerMargin &&
        distances.at(1) <= characteristicDistanceBWithUpperMargin)
    {
      MITK_INFO << "Found Fiducial 4 (PointSet number " << i << ")";
      m_FiducialMarkerCentroids.insert(std::pair<int, mitk::Point3D>(4, m_FiducialCandidates.at(i).first));
      distanceVectorsFiducials.erase(distanceVectorsFiducials.begin() + i);
      m_FiducialCandidates.erase(m_FiducialCandidates.begin() + i);
      return true;
    }
  }
  return false;
}

bool QmitkAutomaticFiducialmarkerRegistrationWidget::FindFiducialNo5()
{
  if (m_FiducialMarkerCentroids.find(2) == m_FiducialMarkerCentroids.end())
  {
    MITK_WARN << "To find fiducial No 5, fiducial No 2 has to be found before.";
    return false;
  }

  double characteristicDistanceBWithUpperMargin = this->GetCharacteristicDistanceBWithUpperMargin();

  mitk::Point3D fiducialNo2(m_FiducialMarkerCentroids.at(2));

  for (unsigned int counter = 0; counter < m_FiducialCandidates.size(); ++counter)
  {
    mitk::Point3D fiducialCentroid(m_FiducialCandidates.at(counter).first);
    double distance = fiducialNo2.EuclideanDistanceTo(fiducialCentroid);
    if (distance <= characteristicDistanceBWithUpperMargin)
    {
      m_FiducialMarkerCentroids[5] = m_FiducialCandidates.at(counter).first;
      m_FiducialCandidates.erase(m_FiducialCandidates.begin() + counter);
      MITK_INFO << "Found Fiducial No 5, PointSet: " << m_FiducialMarkerCentroids[5];
      return true;
    }
  }

  MITK_WARN << "Cannot find fiducial No 5.";
  return false;
}

bool QmitkAutomaticFiducialmarkerRegistrationWidget::FindFiducialNo6()
{
  if (m_FiducialMarkerCentroids.find(5) == m_FiducialMarkerCentroids.end())
  {
    MITK_WARN << "To find fiducial No 6, fiducial No 5 has to be found before.";
    return false;
  }

  double characteristicDistanceAWithUpperMargin = this->GetCharacteristicDistanceAWithUpperMargin();

  mitk::Point3D fiducialNo5(m_FiducialMarkerCentroids.at(5));

  for (unsigned int counter = 0; counter < m_FiducialCandidates.size(); ++counter)
  {
    mitk::Point3D fiducialCentroid(m_FiducialCandidates.at(counter).first);
    double distance = fiducialNo5.EuclideanDistanceTo(fiducialCentroid);
    if (distance <= characteristicDistanceAWithUpperMargin)
    {
      m_FiducialMarkerCentroids[6] = m_FiducialCandidates.at(counter).first;
      m_FiducialCandidates.erase(m_FiducialCandidates.begin() + counter);
      MITK_INFO << "Found Fiducial No 6, PointSet: " << m_FiducialMarkerCentroids[6];
      return true;
    }
  }

  MITK_WARN << "Cannot find fiducial No 6.";
  return false;
}

bool QmitkAutomaticFiducialmarkerRegistrationWidget::FindFiducialNo7()
{
  if (m_FiducialMarkerCentroids.find(8) == m_FiducialMarkerCentroids.end())
  {
    MITK_WARN << "To find fiducial No 7, fiducial No 8 has to be found before.";
    return false;
  }

  double characteristicDistanceAWithUpperMargin = this->GetCharacteristicDistanceAWithUpperMargin();

  mitk::Point3D fiducialNo8(m_FiducialMarkerCentroids.at(8));

  for (unsigned int counter = 0; counter < m_FiducialCandidates.size(); ++counter)
  {
    mitk::Point3D fiducialCentroid(m_FiducialCandidates.at(counter).first);
    double distance = fiducialNo8.EuclideanDistanceTo(fiducialCentroid);
    if (distance <= characteristicDistanceAWithUpperMargin)
    {
      m_FiducialMarkerCentroids[7] = m_FiducialCandidates.at(counter).first;
      m_FiducialCandidates.erase(m_FiducialCandidates.begin() + counter);
      MITK_INFO << "Found Fiducial No 7, PointSet: " << m_FiducialMarkerCentroids[7];
      return true;
    }
  }

  MITK_WARN << "Cannot find fiducial No 7.";
  return false;
}

bool QmitkAutomaticFiducialmarkerRegistrationWidget::FindFiducialNo8()
{
  if (m_FiducialMarkerCentroids.find(3) == m_FiducialMarkerCentroids.end())
  {
    MITK_WARN << "To find fiducial No 8, fiducial No 3 has to be found before.";
    return false;
  }

  double characteristicDistanceBWithUpperMargin = this->GetCharacteristicDistanceBWithUpperMargin();

  mitk::Point3D fiducialNo3(m_FiducialMarkerCentroids.at(3));

  for (unsigned int counter = 0; counter < m_FiducialCandidates.size(); ++counter)
  {
    mitk::Point3D fiducialCentroid(m_FiducialCandidates.at(counter).first);
    double distance = fiducialNo3.EuclideanDistanceTo(fiducialCentroid);
    if (distance <= characteristicDistanceBWithUpperMargin)
    {
      m_FiducialMarkerCentroids[8] = m_FiducialCandidates.at(counter).first;
      m_FiducialCandidates.erase(m_FiducialCandidates.begin() + counter);
      MITK_INFO << "Found Fiducial No 8, PointSet: " << m_FiducialMarkerCentroids[8];
      return true;
    }
  }

  MITK_WARN << "Cannot find fiducial No 8.";
  return false;
}

void QmitkAutomaticFiducialmarkerRegistrationWidget::DefineDataStorageImageFilter()
{
  m_IsAPointSetPredicate = mitk::TNodePredicateDataType<mitk::PointSet>::New();
  mitk::TNodePredicateDataType<mitk::Image>::Pointer isImage = mitk::TNodePredicateDataType<mitk::Image>::New();

  auto isSegmentation = mitk::NodePredicateDataType::New("Segment");
  m_IsASurfacePredicate = mitk::NodePredicateDataType::New("Surface");

  mitk::NodePredicateOr::Pointer validImages = mitk::NodePredicateOr::New();
  validImages->AddPredicate(mitk::NodePredicateAnd::New(isImage, mitk::NodePredicateNot::New(isSegmentation)));

  mitk::NodePredicateNot::Pointer isNotAHelperObject =
    mitk::NodePredicateNot::New(mitk::NodePredicateProperty::New("helper object", mitk::BoolProperty::New(true)));

  m_IsOfTypeImagePredicate = mitk::NodePredicateAnd::New(validImages, isNotAHelperObject);

  mitk::NodePredicateProperty::Pointer isBinaryPredicate =
    mitk::NodePredicateProperty::New("binary", mitk::BoolProperty::New(true));
  mitk::NodePredicateNot::Pointer isNotBinaryPredicate = mitk::NodePredicateNot::New(isBinaryPredicate);

  mitk::NodePredicateAnd::Pointer isABinaryImagePredicate =
    mitk::NodePredicateAnd::New(m_IsOfTypeImagePredicate, isBinaryPredicate);
  mitk::NodePredicateAnd::Pointer isNotABinaryImagePredicate =
    mitk::NodePredicateAnd::New(m_IsOfTypeImagePredicate, isNotBinaryPredicate);

  m_IsASegmentationImagePredicate =
    mitk::NodePredicateOr::New(isABinaryImagePredicate, mitk::TNodePredicateDataType<mitk::LabelSetImage>::New());
  m_IsAPatientImagePredicate = mitk::NodePredicateAnd::New(
    isNotABinaryImagePredicate, mitk::NodePredicateNot::New(mitk::TNodePredicateDataType<mitk::LabelSetImage>::New()));
}
