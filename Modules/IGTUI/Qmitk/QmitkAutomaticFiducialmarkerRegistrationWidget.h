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

#ifndef _QmitkAutomaticFiducialmarkerRegistrationWidget_H_INCLUDED
#define _QmitkAutomaticFiducialmarkerRegistrationWidget_H_INCLUDED


#include "QWidget"
#include "ui_QmitkAutomaticFiducialmarkerRegistrationWidget.h"
#include "MitkIGTUIExports.h"

#include <mitkNodePredicateDataType.h>
#include <mitkNodePredicateAnd.h>
#include <mitkNodePredicateOr.h>

#include <mitkPointSet.h>
#include <mitkImage.h>
#include <mitkPlaneFit.h>

#include <itkImage.h>
#include <itkThresholdImageFilter.h>
#include <itkBinaryThresholdImageFilter.h>
#include <itkGradientMagnitudeImageFilter.h>
#include <itkLaplacianRecursiveGaussianImageFilter.h>
#include "itkVotingBinaryIterativeHoleFillingImageFilter.h"
#include <itkBinaryImageToShapeLabelMapFilter.h>

namespace itk {
  template<class T> class SmartPointer;

}

namespace mitk {
  class NodeDisplacementFilter;
  class NavigationDataSource;
  class PlaneFit;
  class FloatingImageToUltrasoundRegistrationFilter;
}

namespace Ui {
  class QmitkAutomaticFiducialmarkerRegistrationWidget;
}

// Declare typedefs:
typedef itk::Image<int, 3>  ImageType;
typedef itk::ThresholdImageFilter<ImageType> ThresholdImageFilterType;
typedef itk::BinaryThresholdImageFilter <ImageType, ImageType> BinaryThresholdImageFilterType;
typedef itk::LaplacianRecursiveGaussianImageFilter<ImageType, ImageType> LaplacianRecursiveGaussianImageFilterType;
typedef itk::VotingBinaryIterativeHoleFillingImageFilter<ImageType> VotingBinaryIterativeHoleFillingImageFilterType;
typedef itk::BinaryImageToShapeLabelMapFilter<ImageType> BinaryImageToShapeLabelMapFilterType;
typedef BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType ShapeLabelObjectType;
typedef std::pair<mitk::Point3D, BinaryImageToShapeLabelMapFilterType::OutputImageType::LabelObjectType *> FiducialCandidatePairType;




/*!
\brief QmitkAutomaticFiducialmarkerRegistrationWidget

Widget for performing an automatic fiducial marker based registration.

*/
class MITKIGTUI_EXPORT QmitkAutomaticFiducialmarkerRegistrationWidget : public QWidget
{
  Q_OBJECT // this is needed for all Qt objects that should have a MOC object (everything that derives from QObject)
public:

  /*!
  \brief default constructor
  */
  QmitkAutomaticFiducialmarkerRegistrationWidget( QWidget* parent = 0 );

  /*!
  \brief default destructor
  */
  ~QmitkAutomaticFiducialmarkerRegistrationWidget() override;

  void Initialize(itk::SmartPointer<mitk::DataStorage> dataStorage);

  mitk::AffineTransform3D::Pointer GetTransformMarkerCSToImageCS();

  mitk::AffineTransform3D::Pointer GetTransformMarkerCSToSensorCS();

  mitk::AffineTransform3D::Pointer GetInverseTransform(mitk::AffineTransform3D::Pointer transform);

  mitk::DataNode::Pointer GetImageNode();

  mitk::DataNode::Pointer GetSurfaceNode();

  mitk::DataNode::Pointer GetPointSetNode();

  std::vector<mitk::DataNode::Pointer> GetPointSetsDataNodes();



signals:

public slots:




protected:
  void UnsetFloatingImageGeometry();
  void SetImageGeometryInformation(mitk::Image* image);
  double GetVoxelVolume();
  double GetFiducialVolume(double radius);

  bool FilterImage();
  void InitializeImageFilters();

  double GetCharacteristicDistanceAWithUpperMargin();
  double GetCharacteristicDistanceBWithLowerMargin();
  double GetCharacteristicDistanceBWithUpperMargin();
  double GetMinimalFiducialConfigurationDistance();
  void CreateMarkerModelCoordinateSystemPointSet();

  void EliminateTooSmallLabeledObjects(ImageType::Pointer binaryImage);
  bool EliminateFiducialCandidatesByEuclideanDistances();
  void EliminateNearFiducialCandidatesByMaxDistanceToCentroids();
  void ClassifyFiducialCandidates();
  void GetCentroidsOfLabeledObjects();
  void NumerateFiducialMarks();
  void CalculateDistancesBetweenFiducials(std::vector<std::vector<double>> &distanceVectorsFiducials);
  double GetMaxDistanceToCentroidOfFiducialCandidate(unsigned int &index);
  bool FindFiducialNo1(std::vector<std::vector<double>> &distanceVectorsFiducials);
  bool FindFiducialNo2And3();
  bool FindFiducialNo4(std::vector<std::vector<double>> &distanceVectorsFiducials);
  bool FindFiducialNo5();
  bool FindFiducialNo6();
  bool FindFiducialNo7();
  bool FindFiducialNo8();
  void DefineDataStorageImageFilter();

  void CreateQtPartControl( QWidget *parent );
  Ui::QmitkAutomaticFiducialmarkerRegistrationWidget* m_Controls;


protected slots:
  void OnImageToRegisterComboBoxSelectionChanged(const mitk::DataNode *node);
  void OnRegisterMarkerCSToImageCS();
  void OnLocalizeFiducials();
  void OnAddPointSetClicked();
  void OnRemovePointSetClicked();

private:

  mitk::NodePredicateAnd::Pointer m_IsOfTypeImagePredicate;
  mitk::NodePredicateOr::Pointer m_IsASegmentationImagePredicate;
  mitk::NodePredicateAnd::Pointer m_IsAPatientImagePredicate;
  mitk::NodePredicateDataType::Pointer m_IsASurfacePredicate;
  mitk::TNodePredicateDataType<mitk::PointSet>::Pointer m_IsAPointSetPredicate;

  mitk::Image::Pointer m_ImageToRegister;
  mitk::PointSet::Pointer m_MarkerModelCoordinateSystemPointSet;
  mitk::PointSet::Pointer m_MarkerImageToRegisterCoordinateSystemPointSet;
  std::vector<mitk::DataNode::Pointer> m_PointSetsDataNodes;


  ThresholdImageFilterType::Pointer m_ThresholdFilter;
  BinaryThresholdImageFilterType::Pointer m_BinaryThresholdFilter;
  LaplacianRecursiveGaussianImageFilterType::Pointer m_LaplacianFilter1;
  LaplacianRecursiveGaussianImageFilterType::Pointer m_LaplacianFilter2;
  VotingBinaryIterativeHoleFillingImageFilterType::Pointer m_HoleFillingFilter;
  BinaryImageToShapeLabelMapFilterType::Pointer m_BinaryImageToShapeLabelMapFilter;

  //itk::SmartPointer<mitk::FloatingImageToUltrasoundRegistrationFilter> m_FloatingImageToUltrasoundRegistrationFilter;

  std::vector<FiducialCandidatePairType> m_FiducialCandidates;
  std::map<double, mitk::Vector3D> m_EigenVectorsFiducialCandidates;
  std::vector<double> m_EigenValuesFiducialCandidates;
  mitk::Vector3D m_MeanCentroidFiducialCandidates;
  std::map<int, mitk::Point3D> m_FiducialMarkerCentroids;


  mitk::AffineTransform3D::Pointer m_TransformMarkerCSToImageCS;
  mitk::AffineTransform3D::Pointer m_TransformMarkerCSToSensorCS;


  /*!
  \brief The 3D  dimension of the CT or MRI image given in index size.
  */
  mitk::Vector3D m_ImageDimension;
  mitk::Vector3D m_ImageSpacing;

  itk::SmartPointer<mitk::DataStorage> m_DataStorage;
};
#endif // _QmitkAutomaticFiducialmarkerRegistrationWidget_H_INCLUDED
