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

//misc
#define _USE_MATH_DEFINES
#include <math.h>

// Blueberry
#include <berryISelectionService.h>
#include <berryIWorkbenchWindow.h>

// Qmitk
#include "QmitkFiberfoxView.h"

// MITK
#include <mitkImage.h>
#include <mitkDiffusionImage.h>
#include <mitkImageToItk.h>
#include <itkDwiPhantomGenerationFilter.h>
#include <mitkImageCast.h>
#include <mitkProperties.h>
#include <mitkPlanarFigureInteractor.h>
#include <mitkDataStorage.h>
#include <itkFibersFromPlanarFiguresFilter.h>
#include <itkTractsToDWIImageFilter.h>
#include <mitkTensorImage.h>
#include <mitkILinkedRenderWindowPart.h>
#include <mitkGlobalInteraction.h>
#include <mitkImageToItk.h>
#include <mitkImageCast.h>
#include <mitkImageGenerator.h>
#include <mitkNodePredicateDataType.h>
#include <mitkTensorModel.h>
#include <mitkBallModel.h>
#include <mitkStickModel.h>
#include <mitkAstroStickModel.h>
#include <mitkDotModel.h>
#include <mitkRicianNoiseModel.h>
#include <mitkGibbsRingingArtifact.h>
#include <mitkSignalDecay.h>
#include <itkScalableAffineTransform.h>
#include <mitkLevelWindowProperty.h>
#include <mitkNodePredicateOr.h>
#include <mitkNodePredicateAnd.h>
#include <mitkNodePredicateNot.h>

#include <QMessageBox>

#define _USE_MATH_DEFINES
#include <math.h>

const std::string QmitkFiberfoxView::VIEW_ID = "org.mitk.views.fiberfoxview";

QmitkFiberfoxView::QmitkFiberfoxView()
    : QmitkAbstractView()
    , m_Controls( 0 )
    , m_SelectedImage( NULL )
{

}

// Destructor
QmitkFiberfoxView::~QmitkFiberfoxView()
{

}

void QmitkFiberfoxView::CreateQtPartControl( QWidget *parent )
{
    // build up qt view, unless already done
    if ( !m_Controls )
    {
        // create GUI widgets from the Qt Designer's .ui file
        m_Controls = new Ui::QmitkFiberfoxViewControls;
        m_Controls->setupUi( parent );

        m_Controls->m_StickWidget1->setVisible(true);
        m_Controls->m_StickWidget2->setVisible(false);
        m_Controls->m_ZeppelinWidget1->setVisible(false);
        m_Controls->m_ZeppelinWidget2->setVisible(false);
        m_Controls->m_TensorWidget1->setVisible(false);
        m_Controls->m_TensorWidget2->setVisible(false);

        m_Controls->m_BallWidget1->setVisible(true);
        m_Controls->m_BallWidget2->setVisible(false);
        m_Controls->m_AstrosticksWidget1->setVisible(false);
        m_Controls->m_AstrosticksWidget2->setVisible(false);
        m_Controls->m_DotWidget1->setVisible(false);
        m_Controls->m_DotWidget2->setVisible(false);

        m_Controls->m_Comp4FractionFrame->setVisible(false);
        m_Controls->m_DiffusionPropsMessage->setVisible(false);
        m_Controls->m_GeometryMessage->setVisible(false);
        m_Controls->m_AdvancedSignalOptionsFrame->setVisible(false);
        m_Controls->m_AdvancedFiberOptionsFrame->setVisible(false);
        m_Controls->m_VarianceBox->setVisible(false);
        m_Controls->m_GibbsRingingFrame->setVisible(false);
        m_Controls->m_NoiseFrame->setVisible(true);
        m_Controls->m_GhostFrame->setVisible(false);
        m_Controls->m_DistortionsFrame->setVisible(false);

        m_Controls->m_FrequencyMapBox->SetDataStorage(this->GetDataStorage());
        mitk::TNodePredicateDataType<mitk::Image>::Pointer isMitkImage = mitk::TNodePredicateDataType<mitk::Image>::New();
        mitk::NodePredicateDataType::Pointer isDwi = mitk::NodePredicateDataType::New("DiffusionImage");
        mitk::NodePredicateDataType::Pointer isDti = mitk::NodePredicateDataType::New("TensorImage");
        mitk::NodePredicateDataType::Pointer isQbi = mitk::NodePredicateDataType::New("QBallImage");
        mitk::NodePredicateOr::Pointer isDiffusionImage = mitk::NodePredicateOr::New(isDwi, isDti);
        isDiffusionImage = mitk::NodePredicateOr::New(isDiffusionImage, isQbi);
        mitk::NodePredicateNot::Pointer noDiffusionImage = mitk::NodePredicateNot::New(isDiffusionImage);
        mitk::NodePredicateAnd::Pointer finalPredicate = mitk::NodePredicateAnd::New(isMitkImage, noDiffusionImage);
        m_Controls->m_FrequencyMapBox->SetPredicate(finalPredicate);

        connect((QObject*) m_Controls->m_GenerateImageButton, SIGNAL(clicked()), (QObject*) this, SLOT(GenerateImage()));
        connect((QObject*) m_Controls->m_GenerateFibersButton, SIGNAL(clicked()), (QObject*) this, SLOT(GenerateFibers()));
        connect((QObject*) m_Controls->m_CircleButton, SIGNAL(clicked()), (QObject*) this, SLOT(OnDrawROI()));
        connect((QObject*) m_Controls->m_FlipButton, SIGNAL(clicked()), (QObject*) this, SLOT(OnFlipButton()));
        connect((QObject*) m_Controls->m_JoinBundlesButton, SIGNAL(clicked()), (QObject*) this, SLOT(JoinBundles()));
        connect((QObject*) m_Controls->m_VarianceBox, SIGNAL(valueChanged(double)), (QObject*) this, SLOT(OnVarianceChanged(double)));
        connect((QObject*) m_Controls->m_DistributionBox, SIGNAL(currentIndexChanged(int)), (QObject*) this, SLOT(OnDistributionChanged(int)));
        connect((QObject*) m_Controls->m_FiberDensityBox, SIGNAL(valueChanged(int)), (QObject*) this, SLOT(OnFiberDensityChanged(int)));
        connect((QObject*) m_Controls->m_FiberSamplingBox, SIGNAL(valueChanged(int)), (QObject*) this, SLOT(OnFiberSamplingChanged(int)));
        connect((QObject*) m_Controls->m_TensionBox, SIGNAL(valueChanged(double)), (QObject*) this, SLOT(OnTensionChanged(double)));
        connect((QObject*) m_Controls->m_ContinuityBox, SIGNAL(valueChanged(double)), (QObject*) this, SLOT(OnContinuityChanged(double)));
        connect((QObject*) m_Controls->m_BiasBox, SIGNAL(valueChanged(double)), (QObject*) this, SLOT(OnBiasChanged(double)));
        connect((QObject*) m_Controls->m_AddGibbsRinging, SIGNAL(stateChanged(int)), (QObject*) this, SLOT(OnAddGibbsRinging(int)));
        connect((QObject*) m_Controls->m_AddNoise, SIGNAL(stateChanged(int)), (QObject*) this, SLOT(OnAddNoise(int)));
        connect((QObject*) m_Controls->m_AddGhosts, SIGNAL(stateChanged(int)), (QObject*) this, SLOT(OnAddGhosts(int)));
        connect((QObject*) m_Controls->m_AddDistortions, SIGNAL(stateChanged(int)), (QObject*) this, SLOT(OnAddDistortions(int)));

        connect((QObject*) m_Controls->m_ConstantRadiusBox, SIGNAL(stateChanged(int)), (QObject*) this, SLOT(OnConstantRadius(int)));
        connect((QObject*) m_Controls->m_CopyBundlesButton, SIGNAL(clicked()), (QObject*) this, SLOT(CopyBundles()));
        connect((QObject*) m_Controls->m_TransformBundlesButton, SIGNAL(clicked()), (QObject*) this, SLOT(ApplyTransform()));
        connect((QObject*) m_Controls->m_AlignOnGrid, SIGNAL(clicked()), (QObject*) this, SLOT(AlignOnGrid()));

        connect((QObject*) m_Controls->m_Compartment1Box, SIGNAL(currentIndexChanged(int)), (QObject*) this, SLOT(Comp1ModelFrameVisibility(int)));
        connect((QObject*) m_Controls->m_Compartment2Box, SIGNAL(currentIndexChanged(int)), (QObject*) this, SLOT(Comp2ModelFrameVisibility(int)));
        connect((QObject*) m_Controls->m_Compartment3Box, SIGNAL(currentIndexChanged(int)), (QObject*) this, SLOT(Comp3ModelFrameVisibility(int)));
        connect((QObject*) m_Controls->m_Compartment4Box, SIGNAL(currentIndexChanged(int)), (QObject*) this, SLOT(Comp4ModelFrameVisibility(int)));

        connect((QObject*) m_Controls->m_AdvancedOptionsBox, SIGNAL( stateChanged(int)), (QObject*) this, SLOT(ShowAdvancedOptions(int)));
        connect((QObject*) m_Controls->m_AdvancedOptionsBox_2, SIGNAL( stateChanged(int)), (QObject*) this, SLOT(ShowAdvancedOptions(int)));
    }
}

void QmitkFiberfoxView::ShowAdvancedOptions(int state)
{
    if (state)
    {
        m_Controls->m_AdvancedFiberOptionsFrame->setVisible(true);
        m_Controls->m_AdvancedSignalOptionsFrame->setVisible(true);
        m_Controls->m_AdvancedOptionsBox->setChecked(true);
        m_Controls->m_AdvancedOptionsBox_2->setChecked(true);
    }
    else
    {
        m_Controls->m_AdvancedFiberOptionsFrame->setVisible(false);
        m_Controls->m_AdvancedSignalOptionsFrame->setVisible(false);
        m_Controls->m_AdvancedOptionsBox->setChecked(false);
        m_Controls->m_AdvancedOptionsBox_2->setChecked(false);
    }
}

void QmitkFiberfoxView::Comp1ModelFrameVisibility(int index)
{
    m_Controls->m_StickWidget1->setVisible(false);
    m_Controls->m_ZeppelinWidget1->setVisible(false);
    m_Controls->m_TensorWidget1->setVisible(false);

    switch (index)
    {
    case 0:
        m_Controls->m_StickWidget1->setVisible(true);
        break;
    case 1:
        m_Controls->m_ZeppelinWidget1->setVisible(true);
        break;
    case 2:
        m_Controls->m_TensorWidget1->setVisible(true);
        break;
    }
}

void QmitkFiberfoxView::Comp2ModelFrameVisibility(int index)
{
    m_Controls->m_StickWidget2->setVisible(false);
    m_Controls->m_ZeppelinWidget2->setVisible(false);
    m_Controls->m_TensorWidget2->setVisible(false);

    switch (index)
    {
    case 0:
        break;
    case 1:
        m_Controls->m_StickWidget2->setVisible(true);
        break;
    case 2:
        m_Controls->m_ZeppelinWidget2->setVisible(true);
        break;
    case 3:
        m_Controls->m_TensorWidget2->setVisible(true);
        break;
    }
}

void QmitkFiberfoxView::Comp3ModelFrameVisibility(int index)
{
    m_Controls->m_BallWidget1->setVisible(false);
    m_Controls->m_AstrosticksWidget1->setVisible(false);
    m_Controls->m_DotWidget1->setVisible(false);

    switch (index)
    {
    case 0:
        m_Controls->m_BallWidget1->setVisible(true);
        break;
    case 1:
        m_Controls->m_AstrosticksWidget1->setVisible(true);
        break;
    case 2:
        m_Controls->m_DotWidget1->setVisible(true);
        break;
    }
}

void QmitkFiberfoxView::Comp4ModelFrameVisibility(int index)
{
    m_Controls->m_BallWidget2->setVisible(false);
    m_Controls->m_AstrosticksWidget2->setVisible(false);
    m_Controls->m_DotWidget2->setVisible(false);
    m_Controls->m_Comp4FractionFrame->setVisible(false);

    switch (index)
    {
    case 0:
        break;
    case 1:
        m_Controls->m_BallWidget2->setVisible(true);
        m_Controls->m_Comp4FractionFrame->setVisible(true);
        break;
    case 2:
        m_Controls->m_AstrosticksWidget2->setVisible(true);
        m_Controls->m_Comp4FractionFrame->setVisible(true);
        break;
    case 3:
        m_Controls->m_DotWidget2->setVisible(true);
        m_Controls->m_Comp4FractionFrame->setVisible(true);
        break;
    }
}

void QmitkFiberfoxView::OnConstantRadius(int value)
{
    if (value>0 && m_Controls->m_RealTimeFibers->isChecked())
        GenerateFibers();
}

void QmitkFiberfoxView::OnAddDistortions(int value)
{
    if (value>0)
        m_Controls->m_DistortionsFrame->setVisible(true);
    else
        m_Controls->m_DistortionsFrame->setVisible(false);
}

void QmitkFiberfoxView::OnAddGhosts(int value)
{
    if (value>0)
        m_Controls->m_GhostFrame->setVisible(true);
    else
        m_Controls->m_GhostFrame->setVisible(false);
}

void QmitkFiberfoxView::OnAddNoise(int value)
{
    if (value>0)
        m_Controls->m_NoiseFrame->setVisible(true);
    else
        m_Controls->m_NoiseFrame->setVisible(false);
}

void QmitkFiberfoxView::OnAddGibbsRinging(int value)
{
    if (value>0)
        m_Controls->m_GibbsRingingFrame->setVisible(true);
    else
        m_Controls->m_GibbsRingingFrame->setVisible(false);
}

void QmitkFiberfoxView::OnDistributionChanged(int value)
{
    if (value==1)
        m_Controls->m_VarianceBox->setVisible(true);
    else
        m_Controls->m_VarianceBox->setVisible(false);

    if (m_Controls->m_RealTimeFibers->isChecked())
        GenerateFibers();
}

void QmitkFiberfoxView::OnVarianceChanged(double value)
{
    if (m_Controls->m_RealTimeFibers->isChecked())
        GenerateFibers();
}

void QmitkFiberfoxView::OnFiberDensityChanged(int value)
{
    if (m_Controls->m_RealTimeFibers->isChecked())
        GenerateFibers();
}

void QmitkFiberfoxView::OnFiberSamplingChanged(int value)
{
    if (m_Controls->m_RealTimeFibers->isChecked())
        GenerateFibers();
}

void QmitkFiberfoxView::OnTensionChanged(double value)
{
    if (m_Controls->m_RealTimeFibers->isChecked())
        GenerateFibers();
}

void QmitkFiberfoxView::OnContinuityChanged(double value)
{
    if (m_Controls->m_RealTimeFibers->isChecked())
        GenerateFibers();
}

void QmitkFiberfoxView::OnBiasChanged(double value)
{
    if (m_Controls->m_RealTimeFibers->isChecked())
        GenerateFibers();
}

void QmitkFiberfoxView::AlignOnGrid()
{
    for (int i=0; i<m_SelectedFiducials.size(); i++)
    {
        mitk::PlanarEllipse::Pointer pe = dynamic_cast<mitk::PlanarEllipse*>(m_SelectedFiducials.at(i)->GetData());
        mitk::Point3D wc0 = pe->GetWorldControlPoint(0);

        mitk::DataStorage::SetOfObjects::ConstPointer parentFibs = GetDataStorage()->GetSources(m_SelectedFiducials.at(i));
        for( mitk::DataStorage::SetOfObjects::const_iterator it = parentFibs->begin(); it != parentFibs->end(); ++it )
        {
            mitk::DataNode::Pointer pFibNode = *it;
            if ( pFibNode.IsNotNull() && dynamic_cast<mitk::FiberBundleX*>(pFibNode->GetData()) )
            {
                mitk::DataStorage::SetOfObjects::ConstPointer parentImgs = GetDataStorage()->GetSources(pFibNode);
                for( mitk::DataStorage::SetOfObjects::const_iterator it2 = parentImgs->begin(); it2 != parentImgs->end(); ++it2 )
                {
                    mitk::DataNode::Pointer pImgNode = *it2;
                    if ( pImgNode.IsNotNull() && dynamic_cast<mitk::Image*>(pImgNode->GetData()) )
                    {
                        mitk::Image::Pointer img = dynamic_cast<mitk::Image*>(pImgNode->GetData());
                        mitk::Geometry3D::Pointer geom = img->GetGeometry();
                        itk::Index<3> idx;
                        geom->WorldToIndex(wc0, idx);

                        mitk::Point3D cIdx; cIdx[0]=idx[0]; cIdx[1]=idx[1]; cIdx[2]=idx[2];
                        mitk::Point3D world;
                        geom->IndexToWorld(cIdx,world);

                        mitk::Vector3D trans = world - wc0;
                        pe->GetGeometry()->Translate(trans);

                        break;
                    }
                }
                break;
            }
        }
    }

    for( int i=0; i<m_SelectedBundles2.size(); i++ )
    {
        mitk::DataNode::Pointer fibNode = m_SelectedBundles2.at(i);

        mitk::DataStorage::SetOfObjects::ConstPointer sources = GetDataStorage()->GetSources(fibNode);
        for( mitk::DataStorage::SetOfObjects::const_iterator it = sources->begin(); it != sources->end(); ++it )
        {
            mitk::DataNode::Pointer imgNode = *it;
            if ( imgNode.IsNotNull() && dynamic_cast<mitk::Image*>(imgNode->GetData()) )
            {
                mitk::DataStorage::SetOfObjects::ConstPointer derivations = GetDataStorage()->GetDerivations(fibNode);
                for( mitk::DataStorage::SetOfObjects::const_iterator it2 = derivations->begin(); it2 != derivations->end(); ++it2 )
                {
                    mitk::DataNode::Pointer fiducialNode = *it2;
                    if ( fiducialNode.IsNotNull() && dynamic_cast<mitk::PlanarEllipse*>(fiducialNode->GetData()) )
                    {
                        mitk::PlanarEllipse::Pointer pe = dynamic_cast<mitk::PlanarEllipse*>(fiducialNode->GetData());
                        mitk::Point3D wc0 = pe->GetWorldControlPoint(0);

                        mitk::Image::Pointer img = dynamic_cast<mitk::Image*>(imgNode->GetData());
                        mitk::Geometry3D::Pointer geom = img->GetGeometry();
                        itk::Index<3> idx;
                        geom->WorldToIndex(wc0, idx);
                        mitk::Point3D cIdx; cIdx[0]=idx[0]; cIdx[1]=idx[1]; cIdx[2]=idx[2];
                        mitk::Point3D world;
                        geom->IndexToWorld(cIdx,world);

                        mitk::Vector3D trans = world - wc0;
                        pe->GetGeometry()->Translate(trans);
                    }
                }

                break;
            }
        }
    }

    for( int i=0; i<m_SelectedImages.size(); i++ )
    {
        mitk::Image::Pointer img = dynamic_cast<mitk::Image*>(m_SelectedImages.at(i)->GetData());

        mitk::DataStorage::SetOfObjects::ConstPointer derivations = GetDataStorage()->GetDerivations(m_SelectedImages.at(i));
        for( mitk::DataStorage::SetOfObjects::const_iterator it = derivations->begin(); it != derivations->end(); ++it )
        {
            mitk::DataNode::Pointer fibNode = *it;
            if ( fibNode.IsNotNull() && dynamic_cast<mitk::FiberBundleX*>(fibNode->GetData()) )
            {
                mitk::DataStorage::SetOfObjects::ConstPointer derivations2 = GetDataStorage()->GetDerivations(fibNode);
                for( mitk::DataStorage::SetOfObjects::const_iterator it2 = derivations2->begin(); it2 != derivations2->end(); ++it2 )
                {
                    mitk::DataNode::Pointer fiducialNode = *it2;
                    if ( fiducialNode.IsNotNull() && dynamic_cast<mitk::PlanarEllipse*>(fiducialNode->GetData()) )
                    {
                        mitk::PlanarEllipse::Pointer pe = dynamic_cast<mitk::PlanarEllipse*>(fiducialNode->GetData());
                        mitk::Point3D wc0 = pe->GetWorldControlPoint(0);

                        mitk::Geometry3D::Pointer geom = img->GetGeometry();
                        itk::Index<3> idx;
                        geom->WorldToIndex(wc0, idx);
                        mitk::Point3D cIdx; cIdx[0]=idx[0]; cIdx[1]=idx[1]; cIdx[2]=idx[2];
                        mitk::Point3D world;
                        geom->IndexToWorld(cIdx,world);

                        mitk::Vector3D trans = world - wc0;
                        pe->GetGeometry()->Translate(trans);
                    }
                }
            }
        }
    }

    mitk::RenderingManager::GetInstance()->RequestUpdateAll();
    if (m_Controls->m_RealTimeFibers->isChecked())
        GenerateFibers();
}

void QmitkFiberfoxView::OnFlipButton()
{
    if (m_SelectedFiducial.IsNull())
        return;

    std::map<mitk::DataNode*, QmitkPlanarFigureData>::iterator it = m_DataNodeToPlanarFigureData.find(m_SelectedFiducial.GetPointer());
    if( it != m_DataNodeToPlanarFigureData.end() )
    {
        QmitkPlanarFigureData& data = it->second;
        data.m_Flipped += 1;
        data.m_Flipped %= 2;
    }

    if (m_Controls->m_RealTimeFibers->isChecked())
        GenerateFibers();
}

QmitkFiberfoxView::GradientListType QmitkFiberfoxView::GenerateHalfShell(int NPoints)
{
    NPoints *= 2;
    GradientListType pointshell;

    int numB0 = NPoints/20;
    if (numB0==0)
        numB0=1;
    GradientType g;
    g.Fill(0.0);
    for (int i=0; i<numB0; i++)
        pointshell.push_back(g);

    if (NPoints==0)
        return pointshell;

    vnl_vector<double> theta; theta.set_size(NPoints);
    vnl_vector<double> phi; phi.set_size(NPoints);
    double C = sqrt(4*M_PI);
    phi(0) = 0.0;
    phi(NPoints-1) = 0.0;
    for(int i=0; i<NPoints; i++)
    {
        theta(i) = acos(-1.0+2.0*i/(NPoints-1.0)) - M_PI / 2.0;
        if( i>0 && i<NPoints-1)
        {
            phi(i) = (phi(i-1) + C /
                      sqrt(NPoints*(1-(-1.0+2.0*i/(NPoints-1.0))*(-1.0+2.0*i/(NPoints-1.0)))));
            // % (2*DIST_POINTSHELL_PI);
        }
    }

    for(int i=0; i<NPoints; i++)
    {
        g[2] = sin(theta(i));
        if (g[2]<0)
            continue;
        g[0] = cos(theta(i)) * cos(phi(i));
        g[1] = cos(theta(i)) * sin(phi(i));
        pointshell.push_back(g);
    }

    return pointshell;
}

template<int ndirs>
std::vector<itk::Vector<double,3> > QmitkFiberfoxView::MakeGradientList()
{
    std::vector<itk::Vector<double,3> > retval;
    vnl_matrix_fixed<double, 3, ndirs>* U =
            itk::PointShell<ndirs, vnl_matrix_fixed<double, 3, ndirs> >::DistributePointShell();


    // Add 0 vector for B0
    int numB0 = ndirs/10;
    if (numB0==0)
        numB0=1;
    itk::Vector<double,3> v;
    v.Fill(0.0);
    for (int i=0; i<numB0; i++)
    {
        retval.push_back(v);
    }

    for(int i=0; i<ndirs;i++)
    {
        itk::Vector<double,3> v;
        v[0] = U->get(0,i); v[1] = U->get(1,i); v[2] = U->get(2,i);
        retval.push_back(v);
    }

    return retval;
}

void QmitkFiberfoxView::OnAddBundle()
{
    if (m_SelectedImage.IsNull())
        return;

    mitk::DataStorage::SetOfObjects::ConstPointer children = GetDataStorage()->GetDerivations(m_SelectedImage);

    mitk::FiberBundleX::Pointer bundle = mitk::FiberBundleX::New();
    mitk::DataNode::Pointer node = mitk::DataNode::New();
    node->SetData( bundle );
    QString name = QString("Bundle_%1").arg(children->size());
    node->SetName(name.toStdString());
    m_SelectedBundles.push_back(node);
    UpdateGui();

    GetDataStorage()->Add(node, m_SelectedImage);
}

void QmitkFiberfoxView::OnDrawROI()
{
    if (m_SelectedBundles.empty())
        OnAddBundle();
    if (m_SelectedBundles.empty())
        return;

    mitk::DataStorage::SetOfObjects::ConstPointer children = GetDataStorage()->GetDerivations(m_SelectedBundles.at(0));

    mitk::PlanarEllipse::Pointer figure = mitk::PlanarEllipse::New();

    mitk::DataNode::Pointer node = mitk::DataNode::New();
    node->SetData( figure );


    QList<mitk::DataNode::Pointer> nodes = this->GetDataManagerSelection();
    for( int i=0; i<nodes.size(); i++)
        nodes.at(i)->SetSelected(false);

    m_SelectedFiducial = node;

    QString name = QString("Fiducial_%1").arg(children->size());
    node->SetName(name.toStdString());
    node->SetSelected(true);
    GetDataStorage()->Add(node, m_SelectedBundles.at(0));

    this->DisableCrosshairNavigation();

    mitk::PlanarFigureInteractor::Pointer figureInteractor = dynamic_cast<mitk::PlanarFigureInteractor*>(node->GetInteractor());
    if(figureInteractor.IsNull())
        figureInteractor = mitk::PlanarFigureInteractor::New("PlanarFigureInteractor", node);
    mitk::GlobalInteraction::GetInstance()->AddInteractor(figureInteractor);

    UpdateGui();
}

bool CompareLayer(mitk::DataNode::Pointer i,mitk::DataNode::Pointer j)
{
    int li = -1;
    i->GetPropertyValue("layer", li);
    int lj = -1;
    j->GetPropertyValue("layer", lj);
    return li<lj;
}

void QmitkFiberfoxView::GenerateFibers()
{
    if (m_SelectedBundles.empty())
    {
        if (m_SelectedFiducial.IsNull())
            return;

        mitk::DataStorage::SetOfObjects::ConstPointer parents = GetDataStorage()->GetSources(m_SelectedFiducial);
        for( mitk::DataStorage::SetOfObjects::const_iterator it = parents->begin(); it != parents->end(); ++it )
            if(dynamic_cast<mitk::FiberBundleX*>((*it)->GetData()))
                m_SelectedBundles.push_back(*it);

        if (m_SelectedBundles.empty())
            return;
    }

    vector< vector< mitk::PlanarEllipse::Pointer > > fiducials;
    vector< vector< unsigned int > > fliplist;
    for (int i=0; i<m_SelectedBundles.size(); i++)
    {
        mitk::DataStorage::SetOfObjects::ConstPointer children = GetDataStorage()->GetDerivations(m_SelectedBundles.at(i));
        std::vector< mitk::DataNode::Pointer > childVector;
        for( mitk::DataStorage::SetOfObjects::const_iterator it = children->begin(); it != children->end(); ++it )
            childVector.push_back(*it);
        sort(childVector.begin(), childVector.end(), CompareLayer);

        vector< mitk::PlanarEllipse::Pointer > fib;
        vector< unsigned int > flip;
        float radius = 1;
        int count = 0;
        for( std::vector< mitk::DataNode::Pointer >::const_iterator it = childVector.begin(); it != childVector.end(); ++it )
        {
            mitk::DataNode::Pointer node = *it;

            if ( node.IsNotNull() && dynamic_cast<mitk::PlanarEllipse*>(node->GetData()) )
            {
                mitk::PlanarEllipse* ellipse = dynamic_cast<mitk::PlanarEllipse*>(node->GetData());
                if (m_Controls->m_ConstantRadiusBox->isChecked())
                {
                    ellipse->SetTreatAsCircle(true);
                    mitk::Point2D c = ellipse->GetControlPoint(0);
                    mitk::Point2D p = ellipse->GetControlPoint(1);
                    mitk::Vector2D v = p-c;
                    if (count==0)
                    {
                        radius = v.GetVnlVector().magnitude();
                        ellipse->SetControlPoint(1, p);
                    }
                    else
                    {
                        v.Normalize();
                        v *= radius;
                        ellipse->SetControlPoint(1, c+v);
                    }
                }
                fib.push_back(ellipse);

                std::map<mitk::DataNode*, QmitkPlanarFigureData>::iterator it = m_DataNodeToPlanarFigureData.find(node.GetPointer());
                if( it != m_DataNodeToPlanarFigureData.end() )
                {
                    QmitkPlanarFigureData& data = it->second;
                    flip.push_back(data.m_Flipped);
                }
                else
                    flip.push_back(0);
            }
            count++;
        }
        if (fib.size()>1)
        {
            fiducials.push_back(fib);
            fliplist.push_back(flip);
        }
        else if (fib.size()>0)
            m_SelectedBundles.at(i)->SetData( mitk::FiberBundleX::New() );

        mitk::RenderingManager::GetInstance()->RequestUpdateAll();
    }

    itk::FibersFromPlanarFiguresFilter::Pointer filter = itk::FibersFromPlanarFiguresFilter::New();
    filter->SetFiducials(fiducials);
    filter->SetFlipList(fliplist);

    switch(m_Controls->m_DistributionBox->currentIndex()){
    case 0:
        filter->SetFiberDistribution(itk::FibersFromPlanarFiguresFilter::DISTRIBUTE_UNIFORM);
        break;
    case 1:
        filter->SetFiberDistribution(itk::FibersFromPlanarFiguresFilter::DISTRIBUTE_GAUSSIAN);
        filter->SetVariance(m_Controls->m_VarianceBox->value());
        break;
    }

    filter->SetDensity(m_Controls->m_FiberDensityBox->value());
    filter->SetTension(m_Controls->m_TensionBox->value());
    filter->SetContinuity(m_Controls->m_ContinuityBox->value());
    filter->SetBias(m_Controls->m_BiasBox->value());
    filter->SetFiberSampling(m_Controls->m_FiberSamplingBox->value());
    filter->Update();
    vector< mitk::FiberBundleX::Pointer > fiberBundles = filter->GetFiberBundles();

    for (int i=0; i<fiberBundles.size(); i++)
    {
        m_SelectedBundles.at(i)->SetData( fiberBundles.at(i) );
        if (fiberBundles.at(i)->GetNumFibers()>50000)
            m_SelectedBundles.at(i)->SetVisibility(false);
    }

    mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}

void QmitkFiberfoxView::GenerateImage()
{
    itk::ImageRegion<3> imageRegion;
    imageRegion.SetSize(0, m_Controls->m_SizeX->value());
    imageRegion.SetSize(1, m_Controls->m_SizeY->value());
    imageRegion.SetSize(2, m_Controls->m_SizeZ->value());
    mitk::Vector3D spacing;
    spacing[0] = m_Controls->m_SpacingX->value();
    spacing[1] = m_Controls->m_SpacingY->value();
    spacing[2] = m_Controls->m_SpacingZ->value();

    mitk::Point3D   origin;
    origin[0] = spacing[0]/2;
    origin[1] = spacing[1]/2;
    origin[2] = spacing[2]/2;
    itk::Matrix<double, 3, 3>           directionMatrix; directionMatrix.SetIdentity();

    if (m_SelectedBundles.empty())
    {
        mitk::Image::Pointer image = mitk::ImageGenerator::GenerateGradientImage<unsigned int>(
                    m_Controls->m_SizeX->value(),
                    m_Controls->m_SizeY->value(),
                    m_Controls->m_SizeZ->value(),
                    m_Controls->m_SpacingX->value(),
                    m_Controls->m_SpacingY->value(),
                    m_Controls->m_SpacingZ->value());

        mitk::Geometry3D* geom = image->GetGeometry();
        geom->SetOrigin(origin);

        mitk::DataNode::Pointer node = mitk::DataNode::New();
        node->SetData( image );
        node->SetName("Dummy");
        unsigned int window = m_Controls->m_SizeX->value()*m_Controls->m_SizeY->value()*m_Controls->m_SizeZ->value();
        unsigned int level = window/2;
        mitk::LevelWindow lw; lw.SetLevelWindow(level, window);
        node->SetProperty( "levelwindow", mitk::LevelWindowProperty::New( lw ) );
        GetDataStorage()->Add(node);
        m_SelectedImage = node;

        mitk::BaseData::Pointer basedata = node->GetData();
        if (basedata.IsNotNull())
        {
            mitk::RenderingManager::GetInstance()->InitializeViews(
                        basedata->GetTimeSlicedGeometry(), mitk::RenderingManager::REQUEST_UPDATE_ALL, true );
            mitk::RenderingManager::GetInstance()->RequestUpdateAll();
        }
        UpdateGui();

        return;
    }

    if (m_SelectedImage.IsNotNull())
    {
        mitk::Image* img = dynamic_cast<mitk::Image*>(m_SelectedImage->GetData());
        itk::Image< float, 3 >::Pointer itkImg = itk::Image< float, 3 >::New();
        CastToItkImage< itk::Image< float, 3 > >(img, itkImg);

        imageRegion = itkImg->GetLargestPossibleRegion();
        spacing = itkImg->GetSpacing();
        origin = itkImg->GetOrigin();
        directionMatrix = itkImg->GetDirection();
    }

    DiffusionSignalModel<double>::GradientListType gradientList;
    double bVal = 1000;
    if (m_SelectedDWI.IsNull())
    {
        gradientList = GenerateHalfShell(m_Controls->m_NumGradientsBox->value());;
        bVal = m_Controls->m_BvalueBox->value();
    }
    else
    {
        mitk::DiffusionImage<short>::Pointer dwi = dynamic_cast<mitk::DiffusionImage<short>*>(m_SelectedDWI->GetData());
        imageRegion = dwi->GetVectorImage()->GetLargestPossibleRegion();
        spacing = dwi->GetVectorImage()->GetSpacing();
        origin = dwi->GetVectorImage()->GetOrigin();
        directionMatrix = dwi->GetVectorImage()->GetDirection();
        bVal = dwi->GetB_Value();
        mitk::DiffusionImage<short>::GradientDirectionContainerType::Pointer dirs = dwi->GetDirections();
        for (int i=0; i<dirs->Size(); i++)
        {
            DiffusionSignalModel<double>::GradientType g;
            g[0] = dirs->at(i)[0];
            g[1] = dirs->at(i)[1];
            g[2] = dirs->at(i)[2];
            gradientList.push_back(g);
        }
    }


    for (int i=0; i<m_SelectedBundles.size(); i++)
    {
        itk::TractsToDWIImageFilter::Pointer tractsToDwiFilter = itk::TractsToDWIImageFilter::New();

        // storage for generated phantom image
        mitk::DataNode::Pointer resultNode = mitk::DataNode::New();
        QString signalModelString("");
        itk::TractsToDWIImageFilter::DiffusionModelList fiberModelList, nonFiberModelList;

        // signal models
        double comp3Weight = 1;
        double comp4Weight = 0;
        if (m_Controls->m_Compartment4Box->currentIndex()>0)
        {
            comp4Weight = m_Controls->m_Comp4FractionBox->value();
            comp3Weight -= comp4Weight;
        }

        mitk::StickModel<double> stickModel1;
        mitk::StickModel<double> stickModel2;
        mitk::TensorModel<double> zeppelinModel1;
        mitk::TensorModel<double> zeppelinModel2;
        mitk::TensorModel<double> tensorModel1;
        mitk::TensorModel<double> tensorModel2;
        mitk::BallModel<double> ballModel1;
        mitk::BallModel<double> ballModel2;
        mitk::AstroStickModel<double> astrosticksModel1;
        mitk::AstroStickModel<double> astrosticksModel2;
        mitk::DotModel<double> dotModel1;
        mitk::DotModel<double> dotModel2;

        // compartment 1
        switch (m_Controls->m_Compartment1Box->currentIndex())
        {
        case 0:
            MITK_INFO << "Using stick model";
            stickModel1.SetGradientList(gradientList);
            stickModel1.SetDiffusivity(m_Controls->m_StickWidget1->GetD());
            stickModel1.SetT2(m_Controls->m_StickWidget1->GetT2());
            fiberModelList.push_back(&stickModel1);
            signalModelString += "Stick";
            resultNode->AddProperty("Fiberfox.Compartment1.Description", StringProperty::New("Intra-axonal compartment") );
            resultNode->AddProperty("Fiberfox.Compartment1.Model", StringProperty::New("Stick") );
            resultNode->AddProperty("Fiberfox.Compartment1.D", DoubleProperty::New(m_Controls->m_StickWidget1->GetD()) );
            resultNode->AddProperty("Fiberfox.Compartment1.T2", DoubleProperty::New(stickModel1.GetT2()) );
            break;
        case 1:
            MITK_INFO << "Using zeppelin model";
            zeppelinModel1.SetGradientList(gradientList);
            zeppelinModel1.SetBvalue(bVal);
            zeppelinModel1.SetDiffusivity1(m_Controls->m_ZeppelinWidget1->GetD1());
            zeppelinModel1.SetDiffusivity2(m_Controls->m_ZeppelinWidget1->GetD2());
            zeppelinModel1.SetDiffusivity3(m_Controls->m_ZeppelinWidget1->GetD2());
            zeppelinModel1.SetT2(m_Controls->m_ZeppelinWidget1->GetT2());
            fiberModelList.push_back(&zeppelinModel1);
            signalModelString += "Zeppelin";
            resultNode->AddProperty("Fiberfox.Compartment1.Description", StringProperty::New("Intra-axonal compartment") );
            resultNode->AddProperty("Fiberfox.Compartment1.Model", StringProperty::New("Zeppelin") );
            resultNode->AddProperty("Fiberfox.Compartment1.D1", DoubleProperty::New(m_Controls->m_ZeppelinWidget1->GetD1()) );
            resultNode->AddProperty("Fiberfox.Compartment1.D2", DoubleProperty::New(m_Controls->m_ZeppelinWidget1->GetD2()) );
            resultNode->AddProperty("Fiberfox.Compartment1.T2", DoubleProperty::New(zeppelinModel1.GetT2()) );
            break;
        case 2:
            MITK_INFO << "Using tensor model";
            tensorModel1.SetGradientList(gradientList);
            tensorModel1.SetBvalue(bVal);
            tensorModel1.SetDiffusivity1(m_Controls->m_TensorWidget1->GetD1());
            tensorModel1.SetDiffusivity2(m_Controls->m_TensorWidget1->GetD2());
            tensorModel1.SetDiffusivity3(m_Controls->m_TensorWidget1->GetD3());
            tensorModel1.SetT2(m_Controls->m_TensorWidget1->GetT2());
            fiberModelList.push_back(&tensorModel1);
            signalModelString += "Tensor";
            resultNode->AddProperty("Fiberfox.Compartment1.Description", StringProperty::New("Intra-axonal compartment") );
            resultNode->AddProperty("Fiberfox.Compartment1.Model", StringProperty::New("Tensor") );
            resultNode->AddProperty("Fiberfox.Compartment1.D1", DoubleProperty::New(m_Controls->m_TensorWidget1->GetD1()) );
            resultNode->AddProperty("Fiberfox.Compartment1.D2", DoubleProperty::New(m_Controls->m_TensorWidget1->GetD2()) );
            resultNode->AddProperty("Fiberfox.Compartment1.D3", DoubleProperty::New(m_Controls->m_TensorWidget1->GetD3()) );
            resultNode->AddProperty("Fiberfox.Compartment1.T2", DoubleProperty::New(zeppelinModel1.GetT2()) );
            break;
        }

        // compartment 2
        switch (m_Controls->m_Compartment2Box->currentIndex())
        {
        case 0:
            break;
        case 1:
            stickModel2.SetGradientList(gradientList);
            stickModel2.SetDiffusivity(m_Controls->m_StickWidget2->GetD());
            stickModel2.SetT2(m_Controls->m_StickWidget2->GetT2());
            fiberModelList.push_back(&stickModel2);
            signalModelString += "Stick";
            resultNode->AddProperty("Fiberfox.Compartment2.Description", StringProperty::New("Inter-axonal compartment") );
            resultNode->AddProperty("Fiberfox.Compartment2.Model", StringProperty::New("Stick") );
            resultNode->AddProperty("Fiberfox.Compartment2.D", DoubleProperty::New(m_Controls->m_StickWidget2->GetD()) );
            resultNode->AddProperty("Fiberfox.Compartment2.T2", DoubleProperty::New(stickModel2.GetT2()) );
            break;
        case 2:
            zeppelinModel2.SetGradientList(gradientList);
            zeppelinModel2.SetBvalue(bVal);
            zeppelinModel2.SetDiffusivity1(m_Controls->m_ZeppelinWidget2->GetD1());
            zeppelinModel2.SetDiffusivity2(m_Controls->m_ZeppelinWidget2->GetD2());
            zeppelinModel2.SetDiffusivity3(m_Controls->m_ZeppelinWidget2->GetD2());
            zeppelinModel2.SetT2(m_Controls->m_ZeppelinWidget2->GetT2());
            fiberModelList.push_back(&zeppelinModel2);
            signalModelString += "Zeppelin";
            resultNode->AddProperty("Fiberfox.Compartment2.Description", StringProperty::New("Inter-axonal compartment") );
            resultNode->AddProperty("Fiberfox.Compartment2.Model", StringProperty::New("Zeppelin") );
            resultNode->AddProperty("Fiberfox.Compartment2.D1", DoubleProperty::New(m_Controls->m_ZeppelinWidget2->GetD1()) );
            resultNode->AddProperty("Fiberfox.Compartment2.D2", DoubleProperty::New(m_Controls->m_ZeppelinWidget2->GetD2()) );
            resultNode->AddProperty("Fiberfox.Compartment2.T2", DoubleProperty::New(zeppelinModel2.GetT2()) );
            break;
        case 3:
            tensorModel2.SetGradientList(gradientList);
            tensorModel2.SetBvalue(bVal);
            tensorModel2.SetDiffusivity1(m_Controls->m_TensorWidget2->GetD1());
            tensorModel2.SetDiffusivity2(m_Controls->m_TensorWidget2->GetD2());
            tensorModel2.SetDiffusivity3(m_Controls->m_TensorWidget2->GetD3());
            tensorModel2.SetT2(m_Controls->m_TensorWidget2->GetT2());
            fiberModelList.push_back(&tensorModel2);
            signalModelString += "Tensor";
            resultNode->AddProperty("Fiberfox.Compartment2.Description", StringProperty::New("Inter-axonal compartment") );
            resultNode->AddProperty("Fiberfox.Compartment2.Model", StringProperty::New("Tensor") );
            resultNode->AddProperty("Fiberfox.Compartment2.D1", DoubleProperty::New(m_Controls->m_TensorWidget2->GetD1()) );
            resultNode->AddProperty("Fiberfox.Compartment2.D2", DoubleProperty::New(m_Controls->m_TensorWidget2->GetD2()) );
            resultNode->AddProperty("Fiberfox.Compartment2.D3", DoubleProperty::New(m_Controls->m_TensorWidget2->GetD3()) );
            resultNode->AddProperty("Fiberfox.Compartment2.T2", DoubleProperty::New(zeppelinModel2.GetT2()) );
            break;
        }

        // compartment 3
        switch (m_Controls->m_Compartment3Box->currentIndex())
        {
        case 0:
            ballModel1.SetGradientList(gradientList);
            ballModel1.SetBvalue(bVal);
            ballModel1.SetDiffusivity(m_Controls->m_BallWidget1->GetD());
            ballModel1.SetT2(m_Controls->m_BallWidget1->GetT2());
            ballModel1.SetWeight(comp3Weight);
            nonFiberModelList.push_back(&ballModel1);
            signalModelString += "Ball";
            resultNode->AddProperty("Fiberfox.Compartment3.Description", StringProperty::New("Extra-axonal compartment 1") );
            resultNode->AddProperty("Fiberfox.Compartment3.Model", StringProperty::New("Ball") );
            resultNode->AddProperty("Fiberfox.Compartment3.D", DoubleProperty::New(m_Controls->m_BallWidget1->GetD()) );
            resultNode->AddProperty("Fiberfox.Compartment3.T2", DoubleProperty::New(ballModel1.GetT2()) );
            break;
        case 1:
            astrosticksModel1.SetGradientList(gradientList);
            astrosticksModel1.SetBvalue(bVal);
            astrosticksModel1.SetDiffusivity(m_Controls->m_AstrosticksWidget1->GetD());
            astrosticksModel1.SetT2(m_Controls->m_AstrosticksWidget1->GetT2());
            astrosticksModel1.SetRandomizeSticks(m_Controls->m_AstrosticksWidget1->GetRandomizeSticks());
            astrosticksModel1.SetWeight(comp3Weight);
            nonFiberModelList.push_back(&astrosticksModel1);
            signalModelString += "Astrosticks";
            resultNode->AddProperty("Fiberfox.Compartment3.Description", StringProperty::New("Extra-axonal compartment 1") );
            resultNode->AddProperty("Fiberfox.Compartment3.Model", StringProperty::New("Astrosticks") );
            resultNode->AddProperty("Fiberfox.Compartment3.D", DoubleProperty::New(m_Controls->m_AstrosticksWidget1->GetD()) );
            resultNode->AddProperty("Fiberfox.Compartment3.T2", DoubleProperty::New(astrosticksModel1.GetT2()) );
            resultNode->AddProperty("Fiberfox.Compartment3.RandomSticks", BoolProperty::New(m_Controls->m_AstrosticksWidget1->GetRandomizeSticks()) );
            break;
        case 2:
            dotModel1.SetGradientList(gradientList);
            dotModel1.SetT2(m_Controls->m_DotWidget1->GetT2());
            dotModel1.SetWeight(comp3Weight);
            nonFiberModelList.push_back(&dotModel1);
            signalModelString += "Dot";
            resultNode->AddProperty("Fiberfox.Compartment3.Description", StringProperty::New("Extra-axonal compartment 1") );
            resultNode->AddProperty("Fiberfox.Compartment3.Model", StringProperty::New("Dot") );
            resultNode->AddProperty("Fiberfox.Compartment3.T2", DoubleProperty::New(dotModel1.GetT2()) );
            break;
        }

        // compartment 4
        switch (m_Controls->m_Compartment4Box->currentIndex())
        {
        case 0:
            break;
        case 1:
            ballModel2.SetGradientList(gradientList);
            ballModel2.SetBvalue(bVal);
            ballModel2.SetDiffusivity(m_Controls->m_BallWidget2->GetD());
            ballModel2.SetT2(m_Controls->m_BallWidget2->GetT2());
            ballModel2.SetWeight(comp4Weight);
            nonFiberModelList.push_back(&ballModel2);
            signalModelString += "Ball";
            resultNode->AddProperty("Fiberfox.Compartment4.Description", StringProperty::New("Extra-axonal compartment 2") );
            resultNode->AddProperty("Fiberfox.Compartment4.Model", StringProperty::New("Ball") );
            resultNode->AddProperty("Fiberfox.Compartment4.D", DoubleProperty::New(m_Controls->m_BallWidget2->GetD()) );
            resultNode->AddProperty("Fiberfox.Compartment4.T2", DoubleProperty::New(ballModel2.GetT2()) );
            break;
        case 2:
            astrosticksModel2.SetGradientList(gradientList);
            astrosticksModel2.SetBvalue(bVal);
            astrosticksModel2.SetDiffusivity(m_Controls->m_AstrosticksWidget2->GetD());
            astrosticksModel2.SetT2(m_Controls->m_AstrosticksWidget2->GetT2());
            astrosticksModel2.SetRandomizeSticks(m_Controls->m_AstrosticksWidget2->GetRandomizeSticks());
            astrosticksModel2.SetWeight(comp4Weight);
            nonFiberModelList.push_back(&astrosticksModel2);
            signalModelString += "Astrosticks";
            resultNode->AddProperty("Fiberfox.Compartment4.Description", StringProperty::New("Extra-axonal compartment 2") );
            resultNode->AddProperty("Fiberfox.Compartment4.Model", StringProperty::New("Astrosticks") );
            resultNode->AddProperty("Fiberfox.Compartment4.D", DoubleProperty::New(m_Controls->m_AstrosticksWidget2->GetD()) );
            resultNode->AddProperty("Fiberfox.Compartment4.T2", DoubleProperty::New(astrosticksModel2.GetT2()) );
            resultNode->AddProperty("Fiberfox.Compartment4.RandomSticks", BoolProperty::New(m_Controls->m_AstrosticksWidget2->GetRandomizeSticks()) );
            break;
        case 3:
            dotModel2.SetGradientList(gradientList);
            dotModel2.SetT2(m_Controls->m_DotWidget2->GetT2());
            dotModel2.SetWeight(comp4Weight);
            nonFiberModelList.push_back(&dotModel2);
            signalModelString += "Dot";
            resultNode->AddProperty("Fiberfox.Compartment4.Description", StringProperty::New("Extra-axonal compartment 2") );
            resultNode->AddProperty("Fiberfox.Compartment4.Model", StringProperty::New("Dot") );
            resultNode->AddProperty("Fiberfox.Compartment4.T2", DoubleProperty::New(dotModel2.GetT2()) );
            break;
        }

        itk::TractsToDWIImageFilter::KspaceArtifactList artifactList;

        // artifact models
        QString artifactModelString("");
        double noiseVariance = 0;
        if (m_Controls->m_AddNoise->isChecked())
        {
            noiseVariance = m_Controls->m_NoiseLevel->value();
            artifactModelString += "_NOISE";
            artifactModelString += QString::number(noiseVariance);
            resultNode->AddProperty("Fiberfox.Noise-Variance", DoubleProperty::New(noiseVariance));
        }
        mitk::RicianNoiseModel<double> noiseModel;
        noiseModel.SetNoiseVariance(noiseVariance);

        mitk::GibbsRingingArtifact<double> gibbsModel;
        if (m_Controls->m_AddGibbsRinging->isChecked())
        {
            artifactModelString += "_RINGING";
            resultNode->AddProperty("Fiberfox.k-Space-Undersampling", IntProperty::New(m_Controls->m_KspaceUndersamplingBox->currentText().toInt()));
            gibbsModel.SetKspaceCropping((double)m_Controls->m_KspaceUndersamplingBox->currentText().toInt());
            artifactList.push_back(&gibbsModel);
        }

        if ( this->m_Controls->m_TEbox->value() < imageRegion.GetSize(1)*m_Controls->m_LineReadoutTimeBox->value() )
        {
            this->m_Controls->m_TEbox->setValue( imageRegion.GetSize(1)*m_Controls->m_LineReadoutTimeBox->value() );
            QMessageBox::information( NULL, "Warning", "Echo time is too short! Time not sufficient to read slice. Automaticall adjusted to "+QString::number(this->m_Controls->m_TEbox->value())+" ms");
        }

        double lineReadoutTime = m_Controls->m_LineReadoutTimeBox->value();

        // adjusting line readout time to the adapted image size needed for the DFT
        int y = imageRegion.GetSize(1);
        if ( y%2 == 1 )
            y += 1;
        if ( y>imageRegion.GetSize(1) )
            lineReadoutTime *= (double)imageRegion.GetSize(1)/y;

        // add signal contrast model
        mitk::SignalDecay<double> contrastModel;
        contrastModel.SetTinhom(this->m_Controls->m_T2starBox->value());
        contrastModel.SetTE(this->m_Controls->m_TEbox->value());
        contrastModel.SetTline(lineReadoutTime);
        artifactList.push_back(&contrastModel);

        // add N/2 ghosting
        double kOffset = 0;
        if (m_Controls->m_AddGhosts->isChecked())
        {
            artifactModelString += "_GHOST";
            kOffset = m_Controls->m_kOffsetBox->value();
            resultNode->AddProperty("Fiberfox.Line-Offset", DoubleProperty::New(kOffset));
        }

        // add distortions
        if (m_Controls->m_AddDistortions->isChecked() && m_Controls->m_FrequencyMapBox->GetSelectedNode().IsNotNull())
        {
            mitk::DataNode::Pointer fMapNode = m_Controls->m_FrequencyMapBox->GetSelectedNode();
            mitk::Image* img = dynamic_cast<mitk::Image*>(fMapNode->GetData());
            ItkDoubleImgType::Pointer itkImg = ItkDoubleImgType::New();
            CastToItkImage< ItkDoubleImgType >(img, itkImg);

            if (imageRegion.GetSize(0)==itkImg->GetLargestPossibleRegion().GetSize(0) &&
                imageRegion.GetSize(1)==itkImg->GetLargestPossibleRegion().GetSize(1) &&
                imageRegion.GetSize(2)==itkImg->GetLargestPossibleRegion().GetSize(2))
                tractsToDwiFilter->SetFrequencyMap(itkImg);
        }

        mitk::FiberBundleX::Pointer fiberBundle = dynamic_cast<mitk::FiberBundleX*>(m_SelectedBundles.at(i)->GetData());
        if (fiberBundle->GetNumFibers()<=0)
            continue;

        tractsToDwiFilter->SetImageRegion(imageRegion);
        tractsToDwiFilter->SetSpacing(spacing);
        tractsToDwiFilter->SetOrigin(origin);
        tractsToDwiFilter->SetDirectionMatrix(directionMatrix);
        tractsToDwiFilter->SetFiberBundle(fiberBundle);
        tractsToDwiFilter->SetFiberModels(fiberModelList);
        tractsToDwiFilter->SetNonFiberModels(nonFiberModelList);
        tractsToDwiFilter->SetNoiseModel(&noiseModel);
        tractsToDwiFilter->SetKspaceArtifacts(artifactList);
        tractsToDwiFilter->SetkOffset(kOffset);
        tractsToDwiFilter->SettLine(m_Controls->m_LineReadoutTimeBox->value());
        tractsToDwiFilter->SetNumberOfRepetitions(m_Controls->m_RepetitionsBox->value());
        tractsToDwiFilter->SetEnforcePureFiberVoxels(m_Controls->m_EnforcePureFiberVoxelsBox->isChecked());
        tractsToDwiFilter->SetInterpolationShrink(m_Controls->m_InterpolationShrink->value());
        tractsToDwiFilter->SetFiberRadius(m_Controls->m_FiberRadius->value());
        tractsToDwiFilter->SetSignalScale(m_Controls->m_SignalScaleBox->value());

        if (m_TissueMask.IsNotNull())
        {
            ItkUcharImgType::Pointer mask = ItkUcharImgType::New();
            mitk::CastToItkImage<ItkUcharImgType>(m_TissueMask, mask);
            tractsToDwiFilter->SetTissueMask(mask);
        }
        tractsToDwiFilter->Update();

        mitk::DiffusionImage<short>::Pointer image = mitk::DiffusionImage<short>::New();
        image->SetVectorImage( tractsToDwiFilter->GetOutput() );
        image->SetB_Value(bVal);
        image->SetDirections(gradientList);
        image->InitializeFromVectorImage();
        resultNode->SetData( image );
        resultNode->SetName(m_SelectedBundles.at(i)->GetName()
                            +"_D"+QString::number(imageRegion.GetSize(0)).toStdString()
                            +"-"+QString::number(imageRegion.GetSize(1)).toStdString()
                            +"-"+QString::number(imageRegion.GetSize(2)).toStdString()
                            +"_S"+QString::number(spacing[0]).toStdString()
                            +"-"+QString::number(spacing[1]).toStdString()
                            +"-"+QString::number(spacing[2]).toStdString()
                            +"_b"+QString::number(bVal).toStdString()
                            +"_"+signalModelString.toStdString()
                            +artifactModelString.toStdString());
        GetDataStorage()->Add(resultNode, m_SelectedBundles.at(i));

        resultNode->AddProperty("Fiberfox.InterpolationShrink", IntProperty::New(m_Controls->m_InterpolationShrink->value()));
        resultNode->AddProperty("Fiberfox.SignalScale", IntProperty::New(m_Controls->m_SignalScaleBox->value()));
        resultNode->AddProperty("Fiberfox.FiberRadius", IntProperty::New(m_Controls->m_FiberRadius->value()));
        resultNode->AddProperty("Fiberfox.Tinhom", IntProperty::New(m_Controls->m_T2starBox->value()));
        resultNode->AddProperty("Fiberfox.Repetitions", IntProperty::New(m_Controls->m_RepetitionsBox->value()));
        resultNode->AddProperty("Fiberfox.b-value", DoubleProperty::New(bVal));
        resultNode->AddProperty("Fiberfox.Model", StringProperty::New(signalModelString.toStdString()));
        resultNode->AddProperty("Fiberfox.PureFiberVoxels", BoolProperty::New(m_Controls->m_EnforcePureFiberVoxelsBox->isChecked()));
        resultNode->AddProperty("binary", BoolProperty::New(false));
        resultNode->SetProperty( "levelwindow", mitk::LevelWindowProperty::New(tractsToDwiFilter->GetLevelWindow()) );

        if (m_Controls->m_KspaceImageBox->isChecked())
        {
            itk::Image<double, 3>::Pointer kspace = tractsToDwiFilter->GetKspaceImage();
            mitk::Image::Pointer image = mitk::Image::New();
            image->InitializeByItk(kspace.GetPointer());
            image->SetVolume(kspace->GetBufferPointer());

            mitk::DataNode::Pointer node = mitk::DataNode::New();
            node->SetData( image );
            node->SetName(m_SelectedBundles.at(i)->GetName()+"_k-space");
            GetDataStorage()->Add(node, m_SelectedBundles.at(i));
        }

        mitk::BaseData::Pointer basedata = resultNode->GetData();
        if (basedata.IsNotNull())
        {
            mitk::RenderingManager::GetInstance()->InitializeViews(
                        basedata->GetTimeSlicedGeometry(), mitk::RenderingManager::REQUEST_UPDATE_ALL, true );
            mitk::RenderingManager::GetInstance()->RequestUpdateAll();
        }
    }
}

void QmitkFiberfoxView::ApplyTransform()
{
    vector< mitk::DataNode::Pointer > selectedBundles;
    for( int i=0; i<m_SelectedImages.size(); i++ )
    {
        mitk::DataStorage::SetOfObjects::ConstPointer derivations = GetDataStorage()->GetDerivations(m_SelectedImages.at(i));
        for( mitk::DataStorage::SetOfObjects::const_iterator it = derivations->begin(); it != derivations->end(); ++it )
        {
            mitk::DataNode::Pointer fibNode = *it;
            if ( fibNode.IsNotNull() && dynamic_cast<mitk::FiberBundleX*>(fibNode->GetData()) )
                selectedBundles.push_back(fibNode);
        }
    }
    if (selectedBundles.empty())
        selectedBundles = m_SelectedBundles2;

    if (!selectedBundles.empty())
    {
        std::vector<mitk::DataNode::Pointer>::const_iterator it = selectedBundles.begin();
        for (it; it!=selectedBundles.end(); ++it)
        {
            mitk::FiberBundleX::Pointer fib = dynamic_cast<mitk::FiberBundleX*>((*it)->GetData());
            fib->RotateAroundAxis(m_Controls->m_XrotBox->value(), m_Controls->m_YrotBox->value(), m_Controls->m_ZrotBox->value());
            fib->TranslateFibers(m_Controls->m_XtransBox->value(), m_Controls->m_YtransBox->value(), m_Controls->m_ZtransBox->value());
            fib->ScaleFibers(m_Controls->m_XscaleBox->value(), m_Controls->m_YscaleBox->value(), m_Controls->m_ZscaleBox->value());

            // handle child fiducials
            if (m_Controls->m_IncludeFiducials->isChecked())
            {
                mitk::DataStorage::SetOfObjects::ConstPointer derivations = GetDataStorage()->GetDerivations(*it);
                for( mitk::DataStorage::SetOfObjects::const_iterator it2 = derivations->begin(); it2 != derivations->end(); ++it2 )
                {
                    mitk::DataNode::Pointer fiducialNode = *it2;
                    if ( fiducialNode.IsNotNull() && dynamic_cast<mitk::PlanarEllipse*>(fiducialNode->GetData()) )
                    {
                        mitk::PlanarEllipse* pe = dynamic_cast<mitk::PlanarEllipse*>(fiducialNode->GetData());
                        mitk::Geometry3D* geom = pe->GetGeometry();

                        // translate
                        mitk::Vector3D world;
                        world[0] = m_Controls->m_XtransBox->value();
                        world[1] = m_Controls->m_YtransBox->value();
                        world[2] = m_Controls->m_ZtransBox->value();
                        geom->Translate(world);

                        // calculate rotation matrix
                        double x = m_Controls->m_XrotBox->value()*M_PI/180;
                        double y = m_Controls->m_YrotBox->value()*M_PI/180;
                        double z = m_Controls->m_ZrotBox->value()*M_PI/180;

                        itk::Matrix< float, 3, 3 > rotX; rotX.SetIdentity();
                        rotX[1][1] = cos(x);
                        rotX[2][2] = rotX[1][1];
                        rotX[1][2] = -sin(x);
                        rotX[2][1] = -rotX[1][2];

                        itk::Matrix< float, 3, 3 > rotY; rotY.SetIdentity();
                        rotY[0][0] = cos(y);
                        rotY[2][2] = rotY[0][0];
                        rotY[0][2] = sin(y);
                        rotY[2][0] = -rotY[0][2];

                        itk::Matrix< float, 3, 3 > rotZ; rotZ.SetIdentity();
                        rotZ[0][0] = cos(z);
                        rotZ[1][1] = rotZ[0][0];
                        rotZ[0][1] = -sin(z);
                        rotZ[1][0] = -rotZ[0][1];

                        itk::Matrix< float, 3, 3 > rot = rotZ*rotY*rotX;

                        // transform control point coordinate into geometry translation
                        geom->SetOrigin(pe->GetWorldControlPoint(0));
                        mitk::Point2D cp; cp.Fill(0.0);
                        pe->SetControlPoint(0, cp);

                        // rotate fiducial
                        geom->GetIndexToWorldTransform()->SetMatrix(rot*geom->GetIndexToWorldTransform()->GetMatrix());

                        // implicit translation
                        mitk::Vector3D trans;
                        trans[0] = geom->GetOrigin()[0]-fib->GetGeometry()->GetCenter()[0];
                        trans[1] = geom->GetOrigin()[1]-fib->GetGeometry()->GetCenter()[1];
                        trans[2] = geom->GetOrigin()[2]-fib->GetGeometry()->GetCenter()[2];
                        mitk::Vector3D newWc = rot*trans;
                        newWc = newWc-trans;
                        geom->Translate(newWc);
                    }
                }
            }
        }
    }
    else
    {
        for (int i=0; i<m_SelectedFiducials.size(); i++)
        {
            mitk::PlanarEllipse* pe = dynamic_cast<mitk::PlanarEllipse*>(m_SelectedFiducials.at(i)->GetData());
            mitk::Geometry3D* geom = pe->GetGeometry();

            // translate
            mitk::Vector3D world;
            world[0] = m_Controls->m_XtransBox->value();
            world[1] = m_Controls->m_YtransBox->value();
            world[2] = m_Controls->m_ZtransBox->value();
            geom->Translate(world);

            // calculate rotation matrix
            double x = m_Controls->m_XrotBox->value()*M_PI/180;
            double y = m_Controls->m_YrotBox->value()*M_PI/180;
            double z = m_Controls->m_ZrotBox->value()*M_PI/180;
            itk::Matrix< float, 3, 3 > rotX; rotX.SetIdentity();
            rotX[1][1] = cos(x);
            rotX[2][2] = rotX[1][1];
            rotX[1][2] = -sin(x);
            rotX[2][1] = -rotX[1][2];
            itk::Matrix< float, 3, 3 > rotY; rotY.SetIdentity();
            rotY[0][0] = cos(y);
            rotY[2][2] = rotY[0][0];
            rotY[0][2] = sin(y);
            rotY[2][0] = -rotY[0][2];
            itk::Matrix< float, 3, 3 > rotZ; rotZ.SetIdentity();
            rotZ[0][0] = cos(z);
            rotZ[1][1] = rotZ[0][0];
            rotZ[0][1] = -sin(z);
            rotZ[1][0] = -rotZ[0][1];
            itk::Matrix< float, 3, 3 > rot = rotZ*rotY*rotX;

            // transform control point coordinate into geometry translation
            geom->SetOrigin(pe->GetWorldControlPoint(0));
            mitk::Point2D cp; cp.Fill(0.0);
            pe->SetControlPoint(0, cp);

            // rotate fiducial
            geom->GetIndexToWorldTransform()->SetMatrix(rot*geom->GetIndexToWorldTransform()->GetMatrix());
        }
        if (m_Controls->m_RealTimeFibers->isChecked())
            GenerateFibers();
    }
    mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}

void QmitkFiberfoxView::CopyBundles()
{
    if ( m_SelectedBundles.size()<1 ){
        QMessageBox::information( NULL, "Warning", "Select at least one fiber bundle!");
        MITK_WARN("QmitkFiberProcessingView") << "Select at least one fiber bundle!";
        return;
    }

    std::vector<mitk::DataNode::Pointer>::const_iterator it = m_SelectedBundles.begin();
    for (it; it!=m_SelectedBundles.end(); ++it)
    {
        // find parent image
        mitk::DataNode::Pointer parentNode;
        mitk::DataStorage::SetOfObjects::ConstPointer parentImgs = GetDataStorage()->GetSources(*it);
        for( mitk::DataStorage::SetOfObjects::const_iterator it2 = parentImgs->begin(); it2 != parentImgs->end(); ++it2 )
        {
            mitk::DataNode::Pointer pImgNode = *it2;
            if ( pImgNode.IsNotNull() && dynamic_cast<mitk::Image*>(pImgNode->GetData()) )
            {
                parentNode = pImgNode;
                break;
            }
        }

        mitk::FiberBundleX::Pointer fib = dynamic_cast<mitk::FiberBundleX*>((*it)->GetData());
        mitk::FiberBundleX::Pointer newBundle = fib->GetDeepCopy();
        QString name((*it)->GetName().c_str());
        name += "_copy";

        mitk::DataNode::Pointer fbNode = mitk::DataNode::New();
        fbNode->SetData(newBundle);
        fbNode->SetName(name.toStdString());
        fbNode->SetVisibility(true);
        if (parentNode.IsNotNull())
            GetDataStorage()->Add(fbNode, parentNode);
        else
            GetDataStorage()->Add(fbNode);

        // copy child fiducials
        if (m_Controls->m_IncludeFiducials->isChecked())
        {
            mitk::DataStorage::SetOfObjects::ConstPointer derivations = GetDataStorage()->GetDerivations(*it);
            for( mitk::DataStorage::SetOfObjects::const_iterator it2 = derivations->begin(); it2 != derivations->end(); ++it2 )
            {
                mitk::DataNode::Pointer fiducialNode = *it2;
                if ( fiducialNode.IsNotNull() && dynamic_cast<mitk::PlanarEllipse*>(fiducialNode->GetData()) )
                {
                    mitk::PlanarEllipse::Pointer pe = mitk::PlanarEllipse::New();
                    pe->DeepCopy(dynamic_cast<mitk::PlanarEllipse*>(fiducialNode->GetData()));
                    mitk::DataNode::Pointer newNode = mitk::DataNode::New();
                    newNode->SetData(pe);
                    newNode->SetName(fiducialNode->GetName());
                    GetDataStorage()->Add(newNode, fbNode);
                }
            }
        }
    }
    mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}

void QmitkFiberfoxView::JoinBundles()
{
    if ( m_SelectedBundles.size()<2 ){
        QMessageBox::information( NULL, "Warning", "Select at least two fiber bundles!");
        MITK_WARN("QmitkFiberProcessingView") << "Select at least two fiber bundles!";
        return;
    }

    std::vector<mitk::DataNode::Pointer>::const_iterator it = m_SelectedBundles.begin();
    mitk::FiberBundleX::Pointer newBundle = dynamic_cast<mitk::FiberBundleX*>((*it)->GetData());
    QString name("");
    name += QString((*it)->GetName().c_str());
    ++it;
    for (it; it!=m_SelectedBundles.end(); ++it)
    {
        newBundle = newBundle->AddBundle(dynamic_cast<mitk::FiberBundleX*>((*it)->GetData()));
        name += "+"+QString((*it)->GetName().c_str());
    }

    mitk::DataNode::Pointer fbNode = mitk::DataNode::New();
    fbNode->SetData(newBundle);
    fbNode->SetName(name.toStdString());
    fbNode->SetVisibility(true);
    GetDataStorage()->Add(fbNode);
    mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}

void QmitkFiberfoxView::UpdateGui()
{
    m_Controls->m_FiberBundleLabel->setText("<font color='red'>mandatory</font>");
    m_Controls->m_GeometryFrame->setEnabled(true);
    m_Controls->m_GeometryMessage->setVisible(false);
    m_Controls->m_DiffusionPropsMessage->setVisible(false);
    m_Controls->m_FiberGenMessage->setVisible(true);

    m_Controls->m_TransformBundlesButton->setEnabled(false);
    m_Controls->m_CopyBundlesButton->setEnabled(false);
    m_Controls->m_GenerateFibersButton->setEnabled(false);
    m_Controls->m_FlipButton->setEnabled(false);
    m_Controls->m_CircleButton->setEnabled(false);
    m_Controls->m_BvalueBox->setEnabled(true);
    m_Controls->m_NumGradientsBox->setEnabled(true);
    m_Controls->m_JoinBundlesButton->setEnabled(false);
    m_Controls->m_AlignOnGrid->setEnabled(false);

    if (m_SelectedFiducial.IsNotNull())
    {
        m_Controls->m_TransformBundlesButton->setEnabled(true);
        m_Controls->m_FlipButton->setEnabled(true);
        m_Controls->m_AlignOnGrid->setEnabled(true);
    }

    if (m_SelectedImage.IsNotNull() || !m_SelectedBundles.empty())
    {
        m_Controls->m_TransformBundlesButton->setEnabled(true);
        m_Controls->m_CircleButton->setEnabled(true);
        m_Controls->m_FiberGenMessage->setVisible(false);
        m_Controls->m_AlignOnGrid->setEnabled(true);
    }

    if (m_TissueMask.IsNotNull() || m_SelectedImage.IsNotNull())
    {
        m_Controls->m_GeometryMessage->setVisible(true);
        m_Controls->m_GeometryFrame->setEnabled(false);
    }

    if (m_SelectedDWI.IsNotNull())
    {
        m_Controls->m_DiffusionPropsMessage->setVisible(true);
        m_Controls->m_BvalueBox->setEnabled(false);
        m_Controls->m_NumGradientsBox->setEnabled(false);
        m_Controls->m_GeometryMessage->setVisible(true);
        m_Controls->m_GeometryFrame->setEnabled(false);
    }

    if (!m_SelectedBundles.empty())
    {
        m_Controls->m_CopyBundlesButton->setEnabled(true);
        m_Controls->m_GenerateFibersButton->setEnabled(true);
        m_Controls->m_FiberBundleLabel->setText(m_SelectedBundles.at(0)->GetName().c_str());

        if (m_SelectedBundles.size()>1)
            m_Controls->m_JoinBundlesButton->setEnabled(true);
    }
}

void QmitkFiberfoxView::OnSelectionChanged( berry::IWorkbenchPart::Pointer, const QList<mitk::DataNode::Pointer>& nodes )
{
    m_SelectedBundles2.clear();
    m_SelectedImages.clear();
    m_SelectedFiducials.clear();
    m_SelectedFiducial = NULL;
    m_TissueMask = NULL;
    m_SelectedBundles.clear();
    m_SelectedImage = NULL;
    m_SelectedDWI = NULL;
    m_Controls->m_TissueMaskLabel->setText("<font color='grey'>optional</font>");

    // iterate all selected objects, adjust warning visibility
    for( int i=0; i<nodes.size(); i++)
    {
        mitk::DataNode::Pointer node = nodes.at(i);

        if ( node.IsNotNull() && dynamic_cast<mitk::DiffusionImage<short>*>(node->GetData()) )
        {
            m_SelectedDWI = node;
            m_SelectedImage = node;
            m_SelectedImages.push_back(node);
        }
        else if( node.IsNotNull() && dynamic_cast<mitk::Image*>(node->GetData()) )
        {
            m_SelectedImages.push_back(node);
            m_SelectedImage = node;
            bool isBinary = false;
            node->GetPropertyValue<bool>("binary", isBinary);
            if (isBinary)
            {
                m_TissueMask = dynamic_cast<mitk::Image*>(node->GetData());
                m_Controls->m_TissueMaskLabel->setText(node->GetName().c_str());
            }
        }
        else if ( node.IsNotNull() && dynamic_cast<mitk::FiberBundleX*>(node->GetData()) )
        {
            m_SelectedBundles2.push_back(node);
            if (m_Controls->m_RealTimeFibers->isChecked())
            {
                m_SelectedBundles.push_back(node);
                mitk::FiberBundleX::Pointer newFib = dynamic_cast<mitk::FiberBundleX*>(node->GetData());
                if (newFib->GetNumFibers()!=m_Controls->m_FiberDensityBox->value())
                    GenerateFibers();
            }
            else
                m_SelectedBundles.push_back(node);
        }
        else if ( node.IsNotNull() && dynamic_cast<mitk::PlanarEllipse*>(node->GetData()) )
        {
            m_SelectedFiducials.push_back(node);
            m_SelectedFiducial = node;
            m_SelectedBundles.clear();
            mitk::DataStorage::SetOfObjects::ConstPointer parents = GetDataStorage()->GetSources(node);
            for( mitk::DataStorage::SetOfObjects::const_iterator it = parents->begin(); it != parents->end(); ++it )
            {
                mitk::DataNode::Pointer pNode = *it;
                if ( pNode.IsNotNull() && dynamic_cast<mitk::FiberBundleX*>(pNode->GetData()) )
                    m_SelectedBundles.push_back(pNode);
            }
        }
    }
    UpdateGui();
}


void QmitkFiberfoxView::EnableCrosshairNavigation()
{
    MITK_DEBUG << "EnableCrosshairNavigation";

    // enable the crosshair navigation
    if (mitk::ILinkedRenderWindowPart* linkedRenderWindow =
            dynamic_cast<mitk::ILinkedRenderWindowPart*>(this->GetRenderWindowPart()))
    {
        MITK_DEBUG << "enabling linked navigation";
        linkedRenderWindow->EnableLinkedNavigation(true);
        //        linkedRenderWindow->EnableSlicingPlanes(true);
    }

    if (m_Controls->m_RealTimeFibers->isChecked())
        GenerateFibers();
}

void QmitkFiberfoxView::DisableCrosshairNavigation()
{
    MITK_DEBUG << "DisableCrosshairNavigation";

    // disable the crosshair navigation during the drawing
    if (mitk::ILinkedRenderWindowPart* linkedRenderWindow =
            dynamic_cast<mitk::ILinkedRenderWindowPart*>(this->GetRenderWindowPart()))
    {
        MITK_DEBUG << "disabling linked navigation";
        linkedRenderWindow->EnableLinkedNavigation(false);
        //        linkedRenderWindow->EnableSlicingPlanes(false);
    }
}

void QmitkFiberfoxView::NodeRemoved(const mitk::DataNode* node)
{
    mitk::DataNode* nonConstNode = const_cast<mitk::DataNode*>(node);
    std::map<mitk::DataNode*, QmitkPlanarFigureData>::iterator it = m_DataNodeToPlanarFigureData.find(nonConstNode);

    if( it != m_DataNodeToPlanarFigureData.end() )
    {
        QmitkPlanarFigureData& data = it->second;

        // remove observers
        data.m_Figure->RemoveObserver( data.m_EndPlacementObserverTag );
        data.m_Figure->RemoveObserver( data.m_SelectObserverTag );
        data.m_Figure->RemoveObserver( data.m_StartInteractionObserverTag );
        data.m_Figure->RemoveObserver( data.m_EndInteractionObserverTag );

        m_DataNodeToPlanarFigureData.erase( it );
    }
}

void QmitkFiberfoxView::NodeAdded( const mitk::DataNode* node )
{
    // add observer for selection in renderwindow
    mitk::PlanarFigure* figure = dynamic_cast<mitk::PlanarFigure*>(node->GetData());
    bool isPositionMarker (false);
    node->GetBoolProperty("isContourMarker", isPositionMarker);
    if( figure && !isPositionMarker )
    {
        MITK_DEBUG << "figure added. will add interactor if needed.";
        mitk::PlanarFigureInteractor::Pointer figureInteractor
                = dynamic_cast<mitk::PlanarFigureInteractor*>(node->GetInteractor());

        mitk::DataNode* nonConstNode = const_cast<mitk::DataNode*>( node );
        if(figureInteractor.IsNull())
        {
            figureInteractor = mitk::PlanarFigureInteractor::New("PlanarFigureInteractor", nonConstNode);
        }
        else
        {
            // just to be sure that the interactor is not added twice
            mitk::GlobalInteraction::GetInstance()->RemoveInteractor(figureInteractor);
        }

        MITK_DEBUG << "adding interactor to globalinteraction";
        mitk::GlobalInteraction::GetInstance()->AddInteractor(figureInteractor);

        MITK_DEBUG << "will now add observers for planarfigure";
        QmitkPlanarFigureData data;
        data.m_Figure = figure;

        //        // add observer for event when figure has been placed
        typedef itk::SimpleMemberCommand< QmitkFiberfoxView > SimpleCommandType;
        //        SimpleCommandType::Pointer initializationCommand = SimpleCommandType::New();
        //        initializationCommand->SetCallbackFunction( this, &QmitkFiberfoxView::PlanarFigureInitialized );
        //        data.m_EndPlacementObserverTag = figure->AddObserver( mitk::EndPlacementPlanarFigureEvent(), initializationCommand );

        // add observer for event when figure is picked (selected)
        typedef itk::MemberCommand< QmitkFiberfoxView > MemberCommandType;
        MemberCommandType::Pointer selectCommand = MemberCommandType::New();
        selectCommand->SetCallbackFunction( this, &QmitkFiberfoxView::PlanarFigureSelected );
        data.m_SelectObserverTag = figure->AddObserver( mitk::SelectPlanarFigureEvent(), selectCommand );

        // add observer for event when interaction with figure starts
        SimpleCommandType::Pointer startInteractionCommand = SimpleCommandType::New();
        startInteractionCommand->SetCallbackFunction( this, &QmitkFiberfoxView::DisableCrosshairNavigation);
        data.m_StartInteractionObserverTag = figure->AddObserver( mitk::StartInteractionPlanarFigureEvent(), startInteractionCommand );

        // add observer for event when interaction with figure starts
        SimpleCommandType::Pointer endInteractionCommand = SimpleCommandType::New();
        endInteractionCommand->SetCallbackFunction( this, &QmitkFiberfoxView::EnableCrosshairNavigation);
        data.m_EndInteractionObserverTag = figure->AddObserver( mitk::EndInteractionPlanarFigureEvent(), endInteractionCommand );

        m_DataNodeToPlanarFigureData[nonConstNode] = data;
    }
}

void QmitkFiberfoxView::PlanarFigureSelected( itk::Object* object, const itk::EventObject& )
{
    mitk::TNodePredicateDataType<mitk::PlanarFigure>::Pointer isPf = mitk::TNodePredicateDataType<mitk::PlanarFigure>::New();

    mitk::DataStorage::SetOfObjects::ConstPointer allPfs = this->GetDataStorage()->GetSubset( isPf );
    for ( mitk::DataStorage::SetOfObjects::const_iterator it = allPfs->begin(); it!=allPfs->end(); ++it)
    {
        mitk::DataNode* node = *it;

        if( node->GetData() == object )
        {
            node->SetSelected(true);
            m_SelectedFiducial = node;
        }
        else
            node->SetSelected(false);
    }
    UpdateGui();
    this->RequestRenderWindowUpdate();
}

void QmitkFiberfoxView::SetFocus()
{
    m_Controls->m_CircleButton->setFocus();
}
