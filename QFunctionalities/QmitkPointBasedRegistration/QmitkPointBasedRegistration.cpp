/*=========================================================================
 
Program:   Medical Imaging & Interaction Toolkit
Module:    $RCSfile$
Language:  C++
Date:      $Date$
Version:   $Revision$
 
Copyright (c) German Cancer Research Center, Division of Medical and
Biological Informatics. All rights reserved.
See MITKCopyright.txt or http://www.mitk.org/copyright.html for details.
 
This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.
 
=========================================================================*/

#include "QmitkPointBasedRegistration.h"
#include "QmitkPointBasedRegistrationControls.h"
#include "PointBasedRegistration.xpm"
#include "QmitkTreeNodeSelector.h"
#include "QmitkPointListWidget.h"

#include <vtkIterativeClosestPointTransform.h>
#include <vtkMatrix4x4.h>
#include <vtkMath.h>

#include <vtkLandmarkTransform.h>
#include <QmitkStdMultiWidget.h>
#include "QmitkCommonFunctionality.h"
#include "qradiobutton.h"
#include "qapplication.h"
#include <qcursor.h>
#include <qinputdialog.h>
#include <qlcdnumber.h>
#include <qlabel.h>
#include "qmessagebox.h"

#include <LandmarkWarping.h>
#include <mitkPointOperation.h>
#include <mitkPositionEvent.h>
#include "mitkOperationEvent.h"
#include "mitkUndoController.h"
#include <mitkPointSetWriter.h>
#include <mitkPointSetReader.h>


QmitkPointBasedRegistration::QmitkPointBasedRegistration(QObject *parent, const char *name, QmitkStdMultiWidget *mitkStdMultiWidget, mitk::DataTreeIteratorBase* it)
: QmitkFunctionality(parent, name, it), m_MultiWidget(mitkStdMultiWidget), m_Controls(NULL), m_FixedLandmarks(NULL), m_MovingLandmarks(NULL), m_MovingNode(NULL), 
  m_FixedNode(NULL), m_SetInvisible(true), m_ShowRedGreen(false), m_Opacity(0.5), m_OriginalOpacity(1.0), m_OldMovingLayer(0), m_NewMovingLayer(0), 
  m_OldMovingLayerSet(false), m_NewMovingLayerSet(false), m_Transformation(0), m_HideFixedImage(false), m_HideMovingImage(false),
  m_OldFixedLabel(""), m_OldMovingLabel(""), m_Deactivated (false)
{
  SetAvailability(true);
}

QmitkPointBasedRegistration::~QmitkPointBasedRegistration()
{
  if (m_FixedPointSetNode.IsNotNull())
    m_FixedPointSetNode->SetProperty("label", mitk::StringProperty::New(m_OldFixedLabel));
  m_Controls->m_FixedPointListWidget->SetPointSetNode(NULL);
  if (m_MovingPointSetNode.IsNotNull())
    m_MovingPointSetNode->SetProperty("label", mitk::StringProperty::New(m_OldMovingLabel));
  m_Controls->m_MovingPointListWidget->SetPointSetNode(NULL);
}

QWidget * QmitkPointBasedRegistration::CreateMainWidget(QWidget* /*parent*/)
{
  return NULL;
}

QWidget * QmitkPointBasedRegistration::CreateControlWidget(QWidget *parent)
{
  if (m_Controls == NULL)
  {
    m_Controls = new QmitkPointBasedRegistrationControls(parent);
  }
  return m_Controls;
}

static bool Filter( mitk::DataTreeNode * node )
{
  if(node==NULL)
    return false;
  mitk::BaseData* data = node->GetData();
  if(data==NULL)
    return false;
  if(!data->IsInitialized())
    return false;
  if(dynamic_cast<mitk::PointSet*>(data)!=NULL)
    return false;
  if(dynamic_cast<mitk::Geometry2DData*>(data)!=NULL)
    return false;
  if(data->GetTimeSlicedGeometry()==NULL)
    return false;
  if(data->GetTimeSlicedGeometry()->GetTimeSteps()==0)
    return false;
  return true;
}

static bool FilterForMoving( mitk::DataTreeNode * node )
{
  if(node==NULL)
    return false;
  mitk::BaseData* data = node->GetData();
  if(data==NULL) 
    return false;
  if(!data->IsInitialized())
    return false;
  if(dynamic_cast<mitk::PointSet*>(data)!=NULL) 
    return false;
  if(dynamic_cast<mitk::Geometry2DData*>(data)!=NULL) 
    return false;
  if(data->GetTimeSlicedGeometry()==NULL)
    return false;
  if(data->GetTimeSlicedGeometry()->GetTimeSteps()==0)
    return false;
  bool selectedFixed = false;
  node->GetBoolProperty("selectedFixedImage", selectedFixed);
  if(selectedFixed)
    return false;
  return true;
}

void QmitkPointBasedRegistration::CreateConnections()
{
  if ( m_Controls )
  {
    connect( (QObject*)(m_Controls->m_FixedPointListWidget->m_SetPoints), SIGNAL(toggled(bool)), (QObject*)(m_Controls->m_MovingPointListWidget), SLOT(DeactivateInteractor(bool)));
    connect( (QObject*)(m_Controls->m_MovingPointListWidget->m_SetPoints), SIGNAL(toggled(bool)), (QObject*)(m_Controls->m_FixedPointListWidget), SLOT(DeactivateInteractor(bool)));
    connect( (QObject*)(m_Controls->m_FixedPointListWidget), SIGNAL(PointListChanged()), this, SLOT(updateFixedLandmarksList()));
    connect( (QObject*)(m_Controls->m_MovingPointListWidget), SIGNAL(PointListChanged()), this, SLOT(updateMovingLandmarksList()));
    connect( (QObject*)(m_Controls->m_FixedSelector), SIGNAL(Activated(mitk::DataTreeIteratorClone)),(QObject*) this, SLOT(FixedSelected(mitk::DataTreeIteratorClone)) );
    connect( (QObject*)(m_Controls->m_MovingSelector), SIGNAL(Activated(mitk::DataTreeIteratorClone)),(QObject*) this, SLOT(MovingSelected(mitk::DataTreeIteratorClone)) );
    connect( (QObject*)(m_Controls->m_FixedPointListWidget->m_SetPoints), SIGNAL(toggled(bool)), this, SLOT(HideMovingImage(bool)));
    connect( (QObject*)(m_Controls->m_MovingPointListWidget->m_SetPoints), SIGNAL(toggled(bool)), this, SLOT(HideFixedImage(bool)));
    connect(m_Controls,SIGNAL(calculateIt(PointBasedRegistrationControlParameters*)),this,SLOT(calculateLandmarkbased(PointBasedRegistrationControlParameters*)));
    connect(m_Controls,SIGNAL(calculateWithICP(PointBasedRegistrationControlParameters*)),this,SLOT(calculateLandmarkbasedWithICP(PointBasedRegistrationControlParameters*)));
    connect(m_Controls,SIGNAL(calculateLandmarkWarping(PointBasedRegistrationControlParameters*)),this,SLOT(calculateLandmarkWarping(PointBasedRegistrationControlParameters*)));
    connect(m_Controls,SIGNAL(SaveModel()),this,SLOT(SaveModel()));
    connect(m_Controls,SIGNAL(undoTransformation()),this,SLOT(UndoTransformation()));
    connect(m_Controls,SIGNAL(redoTransformation()),this,SLOT(RedoTransformation()));
    connect(m_Controls,SIGNAL(showRedGreenValues(bool)),this,SLOT(showRedGreen(bool)));
    connect(m_Controls,SIGNAL(setImagesInvisible(bool)),this,SLOT(setInvisible(bool)));
    connect(m_Controls,SIGNAL(opacityValueChanged(float)),this,SLOT(OpacityUpdate(float)));
    connect(m_Controls,SIGNAL(transformationChanged(int)), this,SLOT(transformationChanged(int)));
    
    connect(m_Controls, SIGNAL( reinitFixed(const mitk::Geometry3D*) ), m_MultiWidget, SLOT( InitializeStandardViews(const mitk::Geometry3D*) ));
    connect(m_Controls, SIGNAL( reinitMoving(const mitk::Geometry3D*) ), m_MultiWidget, SLOT( InitializeStandardViews(const mitk::Geometry3D*) ));
    connect(m_Controls, SIGNAL( globalReinit(mitk::DataTreeIteratorBase*) ), m_MultiWidget, SLOT( InitializeStandardViews(mitk::DataTreeIteratorBase*) ));

    m_Controls->m_FixedSelector->m_FilterFunction = Filter;
    m_Controls->m_MovingSelector->m_FilterFunction = FilterForMoving;
  }
}

QAction * QmitkPointBasedRegistration::CreateAction(QActionGroup *parent)
{
  QAction* action;
  action = new QAction( tr( "PointBasedRegistration" ), QPixmap((const char**)PointBasedRegistration_xpm), tr( "&PointBasedRegistration" ), 0, parent, "PointBasedRegistration" );
  return action;
}

void QmitkPointBasedRegistration::TreeChanged()
{
  if (!m_Deactivated)
  {
    mitk::DataTreeIteratorClone iter = this->GetDataTreeIterator();
    iter->GoToBegin();
    if (m_SetInvisible == true)
    {
      while ( !iter->IsAtEnd() )
      {
        if ( (iter->Get().IsNotNull()) && (iter->Get()->GetProperty("visible")) && dynamic_cast<mitk::Geometry2DData*>(iter->Get()->GetData())==NULL)
        {
          iter->Get()->SetVisibility(false);
        }
        ++iter;
      }
    }
    m_Controls->m_FixedSelector->SetDataTreeNodeIterator(this->GetDataTreeIterator());
    m_Controls->m_MovingSelector->SetDataTreeNodeIterator(this->GetDataTreeIterator());
    this->FixedSelected(iter);
    this->MovingSelected(iter);
  }
}

void QmitkPointBasedRegistration::Activated()
{
  m_Deactivated = false;
  m_MultiWidget->SetWidgetPlanesVisibility(true);
  mitk::DataTreeIteratorClone iter = this->GetDataTreeIterator();
  iter->GoToBegin();
  if (m_SetInvisible == true)
  {
    bool visible = false;
    m_InvisibleNodesList.clear();
    // set all nodes to invisible
    while ( !iter->IsAtEnd() )
    {
      if ( (iter->Get().IsNotNull()) && (iter->Get()->GetProperty("visible")) && dynamic_cast<mitk::Geometry2DData*>(iter->Get()->GetData())==NULL)
      {
        iter->Get()->GetBoolProperty("visible", visible);
        // safe all previous invisible nodes in a set
        if (visible == false)
        {
          m_InvisibleNodesList.insert(iter->Get());
        }
        else
        {
          iter->Get()->SetVisibility(false);
        }
      }
      ++iter;
    }
  }
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
  QmitkFunctionality::Activated();
  m_Controls->m_FixedSelector->SetDataTreeNodeIterator(this->GetDataTreeIterator());
  m_Controls->m_MovingSelector->SetDataTreeNodeIterator(this->GetDataTreeIterator());
  this->clearTransformationLists();
  this->FixedSelected(iter);
  this->MovingSelected(iter);
  this->OpacityUpdate(m_Opacity);
}

void QmitkPointBasedRegistration::Deactivated()
{
  m_Deactivated = true;
  if (m_FixedPointSetNode.IsNotNull())
    m_FixedPointSetNode->SetProperty("label", mitk::StringProperty::New(m_OldFixedLabel));
  m_Controls->m_FixedPointListWidget->SetPointSetNode(NULL);
  if (m_MovingPointSetNode.IsNotNull())
    m_MovingPointSetNode->SetProperty("label", mitk::StringProperty::New(m_OldMovingLabel));
  m_Controls->m_MovingPointListWidget->SetPointSetNode(NULL);
  // reset previous invisible nodes to invisible and previous visible nodes to visible
  mitk::DataTreeIteratorClone iter = this->GetDataTreeIterator();
  iter->GoToBegin();
  if (m_SetInvisible == true)
  {
    std::set<mitk::DataTreeNode*>::iterator setIter;
    while ( !iter->IsAtEnd() )
    {
      if ( (iter->Get().IsNotNull()) && (iter->Get()->GetProperty("visible")) && dynamic_cast<mitk::Geometry2DData*>(iter->Get()->GetData())==NULL)
      {
        setIter = m_InvisibleNodesList.find(iter->Get());
        if (setIter != m_InvisibleNodesList.end())
        {
          iter->Get()->SetVisibility(false);
        }
        else
        {
          iter->Get()->SetVisibility(true);
        }
      }
      ++iter;
    }
    m_InvisibleNodesList.clear();
  }
  this->setImageColor(false);
  if (m_MovingNode != NULL)
  {
    m_MovingNode->SetOpacity(m_OriginalOpacity);
    if (m_OldMovingLayerSet)
    {
      m_MovingNode->SetIntProperty("layer", m_OldMovingLayer);
      m_OldMovingLayerSet = false;
    }
  }
  this->clearTransformationLists();
  if (m_FixedLandmarks.IsNotNull() && m_FixedLandmarks->GetSize() == 0)
  {
    mitk::DataTreeIteratorClone it = mitk::DataTree::GetIteratorToNode( GetDataTreeIterator()->GetTree(), m_FixedPointSetNode );
    it->Remove();
  }
  if (m_MovingLandmarks.IsNotNull() && m_MovingLandmarks->GetSize() == 0)
  {
    mitk::DataTreeIteratorClone it2 = mitk::DataTree::GetIteratorToNode( GetDataTreeIterator()->GetTree(), m_MovingPointSetNode );
    it2->Remove();
  }
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
  m_FixedNode = NULL;
  m_MovingNode = NULL;
  m_FixedLandmarks = NULL;
  m_MovingLandmarks = NULL;
  m_FixedPointSetNode = NULL;
  m_MovingPointSetNode = NULL;
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
  QmitkFunctionality::Deactivated();
}

void QmitkPointBasedRegistration::FixedSelected(mitk::DataTreeIteratorClone /*imageIt*/)
{
  if (m_Controls->m_FixedSelector->GetSelectedNode() != NULL)
  {
    mitk::DataTreeIteratorClone it;
    it = *m_Controls->m_FixedSelector->GetSelectedIterator();
    if (m_FixedNode != it->Get())
    {
      // remove changes on previous selected node
      if (m_FixedNode != NULL)
      {
        this->setImageColor(false);
        m_FixedNode->SetVisibility(false);
        m_FixedNode->SetProperty("selectedFixedImage", mitk::BoolProperty::New(false));
        if (m_FixedPointSetNode.IsNotNull())
        {
          m_FixedPointSetNode->SetVisibility(false);
          m_FixedPointSetNode->SetProperty("label", mitk::StringProperty::New(m_OldFixedLabel));
          m_FixedNode->SetProperty("selectedFixedImage", mitk::BoolProperty::New(false));
        }
      }
      // get selected node
      m_FixedNode = it->Get();
      mitk::ColorProperty::Pointer colorProperty;
      colorProperty = dynamic_cast<mitk::ColorProperty*>(m_FixedNode->GetProperty("color"));
      if ( colorProperty.IsNotNull() )
      {
        m_FixedColor = colorProperty->GetColor();
      }
      this->setImageColor(m_ShowRedGreen);
      m_FixedNode->SetProperty("selectedFixedImage", mitk::BoolProperty::New(true));

      if (m_MovingNode != NULL)
      {
        // safe MovingLayer only if we do not overwrite original MovingLayer
        if (!m_OldMovingLayerSet)
        {
          m_MovingNode->GetIntProperty("layer", m_OldMovingLayer);
          m_OldMovingLayerSet = true;
        }
        m_FixedNode->GetIntProperty("layer", m_NewMovingLayer);
        m_NewMovingLayer += 1;
        m_MovingNode->SetIntProperty("layer", m_NewMovingLayer);
      }
      bool hasPointSetNode = false;
      for (int i = 0; i < it->CountChildren(); ++i)
      {
        it->GoToChild(i);
        mitk::StringProperty::Pointer nameProp = dynamic_cast<mitk::StringProperty*>(it->Get()->GetProperty("name"));
        if(nameProp.IsNotNull() && nameProp->GetValueAsString()=="PointBasedRegistrationNode")
        {
          m_FixedPointSetNode=it->Get();
          m_FixedLandmarks = dynamic_cast<mitk::PointSet*> (m_FixedPointSetNode->GetData());
          it->GoToParent();
          it->RemoveChild(i);
          hasPointSetNode = true;
          break;
        }
        it->GoToParent();
      }
      if (!hasPointSetNode)
      {
        m_FixedLandmarks = mitk::PointSet::New();
        m_FixedPointSetNode = mitk::DataTreeNode::New();
        m_FixedPointSetNode->SetData(m_FixedLandmarks);
      }
      m_FixedPointSetNode->GetStringProperty("label", m_OldFixedLabel);
      m_FixedPointSetNode->SetProperty("label", mitk::StringProperty::New("F "));
      m_Controls->m_FixedPointListWidget->SetPointSetNode(m_FixedPointSetNode);
      it->Add(m_FixedPointSetNode);
      mitk::RenderingManager::GetInstance()->RequestUpdateAll();
    }
    if (!m_HideFixedImage && m_FixedNode != NULL && m_FixedPointSetNode.IsNotNull())
    {
      if (m_FixedNode != NULL)
        m_FixedNode->SetVisibility(true);
      if (m_FixedPointSetNode.IsNotNull())
        m_FixedPointSetNode->SetVisibility(true);
    }
    else
    {
      if (m_FixedNode != NULL)
        m_FixedNode->SetVisibility(false);
    }
  }
  else
  {
    m_FixedNode = NULL;
    if (m_FixedPointSetNode.IsNotNull())
      m_FixedPointSetNode->SetProperty("label", mitk::StringProperty::New(m_OldFixedLabel));
    m_FixedPointSetNode = NULL;
    m_FixedLandmarks = NULL;
   m_Controls->m_FixedPointListWidget->SetPointSetNode(m_FixedPointSetNode);
  }
}

void QmitkPointBasedRegistration::MovingSelected(mitk::DataTreeIteratorClone /*imageIt*/)
{
  if (m_Controls->m_MovingSelector->GetSelectedNode() != NULL)
  {
    mitk::DataTreeIteratorClone it;
    it = *m_Controls->m_MovingSelector->GetSelectedIterator();
    if (m_MovingNode != it->Get())
    {
      if (m_MovingNode != NULL)
      {
        m_MovingNode->SetOpacity(m_OriginalOpacity);
        this->setImageColor(false);
        if (m_MovingNode != m_FixedNode)
        {
          m_MovingNode->SetVisibility(false);
          m_MovingPointSetNode->SetProperty("label", mitk::StringProperty::New(m_OldMovingLabel));
          m_MovingPointSetNode->SetVisibility(false);
        }
        if (m_OldMovingLayerSet)
        {
          m_MovingNode->SetIntProperty("layer", m_OldMovingLayer);
          m_OldMovingLayerSet = false;
        }
      }
      if (m_MovingPointSetNode.IsNotNull())
        m_MovingPointSetNode->SetProperty("label", mitk::StringProperty::New(m_OldMovingLabel));
      m_MovingNode = it->Get();
      mitk::ColorProperty::Pointer colorProperty;
      colorProperty = dynamic_cast<mitk::ColorProperty*>(m_MovingNode->GetProperty("color"));
      if ( colorProperty.IsNotNull() )
      {
        m_MovingColor = colorProperty->GetColor();
      }
      this->setImageColor(m_ShowRedGreen);
      m_MovingNode->GetFloatProperty("opacity", m_OriginalOpacity);
      this->OpacityUpdate(m_Opacity);
      // safe MovingLayer only if we do not overwrite original MovingLayer
      m_MovingNode->GetIntProperty("layer", m_OldMovingLayer);
      m_OldMovingLayerSet = true;
      // change MovingLayer to be one higher than FixedLayer -> MovingImage will be the upper image
      if (m_FixedNode != NULL)
      {
        m_FixedNode->GetIntProperty("layer", m_NewMovingLayer);
        m_NewMovingLayer += 1;
        m_MovingNode->SetIntProperty("layer", m_NewMovingLayer);
      }
      bool hasPointSetNode = false;
      for (int i = 0; i < it->CountChildren(); ++i)
      {
        it->GoToChild(i);
        mitk::StringProperty::Pointer nameProp = dynamic_cast<mitk::StringProperty*>(it->Get()->GetProperty("name"));
        if(nameProp.IsNotNull() && nameProp->GetValueAsString()=="PointBasedRegistrationNode")
        {
          m_MovingPointSetNode=it->Get();
          m_MovingLandmarks = dynamic_cast<mitk::PointSet*> (m_MovingPointSetNode->GetData());
          it->GoToParent();
          it->RemoveChild(i);
          hasPointSetNode = true;
          break;
        }
        it->GoToParent();
      }
      if (!hasPointSetNode)
      {
        m_MovingLandmarks = mitk::PointSet::New();
        m_MovingPointSetNode = mitk::DataTreeNode::New();
        m_MovingPointSetNode->SetData(m_MovingLandmarks);
      }
      it->Add(m_MovingPointSetNode);
      m_MovingPointSetNode->SetVisibility(true);
      m_MovingPointSetNode->GetStringProperty("label", m_OldMovingLabel);
      m_MovingPointSetNode->SetProperty("label", mitk::StringProperty::New("M "));
      m_MovingPointSetNode->SetProperty("layer",mitk::IntProperty::New(m_NewMovingLayer+1));
      if (m_FixedPointSetNode.IsNotNull())
        m_FixedPointSetNode->SetProperty("layer",mitk::IntProperty::New(m_NewMovingLayer+1));
      m_Controls->m_MovingPointListWidget->SetPointSetNode(m_MovingPointSetNode);
      mitk::RenderingManager::GetInstance()->RequestUpdateAll();
      this->clearTransformationLists();
      this->OpacityUpdate(m_Opacity);
    }
    if (!m_HideMovingImage)
    {
      if (m_MovingNode != NULL)
        m_MovingNode->SetVisibility(true);
      if (m_MovingPointSetNode.IsNotNull())
        m_MovingPointSetNode->SetVisibility(true);
    }
    else
    {
      if (m_MovingNode != NULL)
        m_MovingNode->SetVisibility(false);
    }
  }
  else
  {
    m_MovingNode = NULL;
    if (m_MovingPointSetNode.IsNotNull())
      m_MovingPointSetNode->SetProperty("label", mitk::StringProperty::New(m_OldMovingLabel));
    m_MovingPointSetNode = NULL;
    m_MovingLandmarks = NULL;
    m_Controls->m_MovingPointListWidget->SetPointSetNode(m_MovingPointSetNode);
  }
}

void QmitkPointBasedRegistration::updateMovingLandmarksList()
{
  if(m_MovingPointSetNode.IsNotNull())
  {
    m_MovingLandmarks = (mitk::PointSet*)(m_MovingPointSetNode->GetData());
  }
  this->clearTransformationLists();
  this->checkLandmarkError();
  this->CheckCalculate();
}

void QmitkPointBasedRegistration::updateFixedLandmarksList()
{
  if(m_FixedPointSetNode.IsNotNull())
  {
    m_FixedLandmarks = (mitk::PointSet*)(m_FixedPointSetNode->GetData());
  }
  this->clearTransformationLists();
  this->checkLandmarkError();
  this->CheckCalculate();
}

void QmitkPointBasedRegistration::HideFixedImage(bool hide)
{ 
  m_HideFixedImage = hide;
  if(m_FixedNode != NULL)
  {
    m_FixedNode->SetVisibility(!hide);
  }
  if (hide)
  {
    m_Controls->reinitMovingClicked();
  }
  if (!m_HideMovingImage && !m_HideFixedImage)
  {
    m_Controls->globalReinitClicked();
  }
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}

void QmitkPointBasedRegistration::HideMovingImage(bool hide)
{ 
  m_HideMovingImage = hide;
  if(m_MovingNode != NULL)
  {
    m_MovingNode->SetVisibility(!hide);
  }
  if (hide)
  {
    m_Controls->reinitFixedClicked();
  }
  if (!m_HideMovingImage && !m_HideFixedImage)
  {
    m_Controls->globalReinitClicked();
  }
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}

bool QmitkPointBasedRegistration::CheckCalculate()
{
  if((m_MovingPointSetNode.IsNull())||(m_FixedPointSetNode.IsNull()||m_FixedLandmarks.IsNull()||m_MovingLandmarks.IsNull())) 
     return false;
  if(m_MovingNode==m_FixedNode)
    return false;
  return m_Controls->checkCalculateEnabled(m_FixedLandmarks->GetSize(), m_MovingLandmarks->GetSize());
}

void QmitkPointBasedRegistration::SaveModel()
{
  if(m_MovingNode != NULL)
  {
    mitk::BaseData::Pointer data=m_MovingNode->GetData();
    if (data.IsNotNull())
    {
      CommonFunctionality::SaveBaseData( data.GetPointer(), "RegistrationResult");
    }
  }
}

void QmitkPointBasedRegistration::UndoTransformation()
{
  if(!m_UndoPointsGeometryList.empty())
  {
    mitk::AffineGeometryFrame3D::Pointer movingLandmarksGeometry = m_MovingLandmarks->GetGeometry(0)->Clone();
    m_RedoPointsGeometryList.push_back(static_cast<mitk::Geometry3D *>(movingLandmarksGeometry.GetPointer()));
    m_MovingLandmarks->SetGeometry(m_UndoPointsGeometryList.back());
    m_UndoPointsGeometryList.pop_back();
    //\FIXME when geometry is substituted the matrix referenced by the actor created by the mapper
    //is still pointing to the old one. Workaround: delete mapper
    m_MovingPointSetNode->SetMapper(1, NULL);

    mitk::BaseData::Pointer movingData = m_MovingNode->GetData();
    mitk::AffineGeometryFrame3D::Pointer movingGeometry = movingData->GetGeometry(0)->Clone();
    m_RedoGeometryList.push_back(static_cast<mitk::Geometry3D *>(movingGeometry.GetPointer()));
    movingData->SetGeometry(m_UndoGeometryList.back());
    m_UndoGeometryList.pop_back();
    //\FIXME when geometry is substituted the matrix referenced by the actor created by the mapper
    //is still pointing to the old one. Workaround: delete mapper
    m_MovingNode->SetMapper(1, NULL);
    mitk::RenderingManager::GetInstance()->RequestUpdate(m_MultiWidget->mitkWidget4->GetRenderWindow());

    movingData->GetTimeSlicedGeometry()->UpdateInformation();
    m_MovingLandmarks->GetTimeSlicedGeometry()->UpdateInformation();
    m_Controls->setRedoTransformationEnabled(true);
    mitk::RenderingManager::GetInstance()->RequestUpdateAll();
    this->checkLandmarkError();
    
  }
  if(!m_UndoPointsGeometryList.empty())
  {
    m_Controls->setUndoTransformationEnabled(true);
  }
  else
  {
    m_Controls->setUndoTransformationEnabled(false);
  }
}

void QmitkPointBasedRegistration::RedoTransformation()
{
  if(!m_RedoPointsGeometryList.empty())
  {
    mitk::AffineGeometryFrame3D::Pointer movingLandmarksGeometry = m_MovingLandmarks->GetGeometry(0)->Clone();
    m_UndoPointsGeometryList.push_back(static_cast<mitk::Geometry3D *>(movingLandmarksGeometry.GetPointer()));
    m_MovingLandmarks->SetGeometry(m_RedoPointsGeometryList.back());
    m_RedoPointsGeometryList.pop_back();
    //\FIXME when geometry is substituted the matrix referenced by the actor created by the mapper
    //is still pointing to the old one. Workaround: delete mapper
    m_MovingPointSetNode->SetMapper(1, NULL);

    mitk::BaseData::Pointer movingData = m_MovingNode->GetData();
    mitk::AffineGeometryFrame3D::Pointer movingGeometry = movingData->GetGeometry(0)->Clone();
    m_UndoGeometryList.push_back(static_cast<mitk::Geometry3D *>(movingGeometry.GetPointer()));
    movingData->SetGeometry(m_RedoGeometryList.back());
    m_RedoGeometryList.pop_back();
    //\FIXME when geometry is substituted the matrix referenced by the actor created by the mapper
    //is still pointing to the old one. Workaround: delete mapper
    m_MovingNode->SetMapper(1, NULL);
    mitk::RenderingManager::GetInstance()->RequestUpdate(m_MultiWidget->mitkWidget4->GetRenderWindow());

    movingData->GetTimeSlicedGeometry()->UpdateInformation();
    m_MovingLandmarks->GetTimeSlicedGeometry()->UpdateInformation();
    m_Controls->setUndoTransformationEnabled(true);
    mitk::RenderingManager::GetInstance()->RequestUpdateAll();
    this->checkLandmarkError();
  }
  if(!m_RedoPointsGeometryList.empty())
  {
    m_Controls->setRedoTransformationEnabled(true);
  }
  else
  {
    m_Controls->setRedoTransformationEnabled(false);
  }
}

void QmitkPointBasedRegistration::showRedGreen(bool redGreen)
{
  m_ShowRedGreen = redGreen;
  this->setImageColor(m_ShowRedGreen);
}

void QmitkPointBasedRegistration::setImageColor(bool redGreen)
{
  if (!redGreen && m_FixedNode != NULL)
  {
    m_FixedNode->SetColor(m_FixedColor);
  }
  if (!redGreen && m_MovingNode != NULL)
  {
    m_MovingNode->SetColor(m_MovingColor);
  }
  if (redGreen && m_FixedNode != NULL)
  {
    m_FixedNode->SetColor(1.0f, 0.0f, 0.0f);
  }
  if (redGreen && m_MovingNode != NULL)
  {
    m_MovingNode->SetColor(0.0f, 1.0f, 0.0f);
  }
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}

void QmitkPointBasedRegistration::OpacityUpdate(float opacity)
{
  m_Opacity = opacity;
  if (m_MovingNode != NULL)
  {
    m_MovingNode->SetOpacity(m_Opacity);
  }
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}

void QmitkPointBasedRegistration::clearTransformationLists()
{
  m_Controls->setUndoTransformationEnabled(false);
  m_Controls->setRedoTransformationEnabled(false);
  m_Controls->m_MeanErrorLCD->hide();
  m_Controls->m_MeanError->hide();
  m_UndoGeometryList.clear();
  m_UndoPointsGeometryList.clear();
  m_RedoGeometryList.clear();
  m_RedoPointsGeometryList.clear();
}

void QmitkPointBasedRegistration::checkLandmarkError()
{
  double totalDist = 0, dist = 0, dist2 = 0;
  mitk::Point3D point1, point2, point3;
  double p1[3], p2[3];
  if(m_Transformation < 3)
  {
    if (m_MovingLandmarks.IsNotNull() && m_FixedLandmarks.IsNotNull()&& m_MovingLandmarks->GetSize() != 0 && m_FixedLandmarks->GetSize() != 0)
    {
      for(int pointId = 0; pointId < m_MovingLandmarks->GetSize(); ++pointId)
      {
        point1 = m_MovingLandmarks->GetPoint(pointId);
        point2 = m_FixedLandmarks->GetPoint(0);
        p1[0] = point1[0]; p1[1] = point1[1]; p1[2] = point1[2];
        p2[0] = point2[0]; p2[1] = point2[1]; p2[2] = point2[2];
        dist = vtkMath::Distance2BetweenPoints(p1, p2);
        for(int pointId2 = 1; pointId2 < m_FixedLandmarks->GetSize(); ++pointId2)
        {
          point2 = m_FixedLandmarks->GetPoint(pointId2);
          p1[0] = point1[0]; p1[1] = point1[1]; p1[2] = p1[2];
          p2[0] = point2[0]; p2[1] = point2[1]; p2[2] = p2[2];
          dist2 = vtkMath::Distance2BetweenPoints(p1, p2);
          if (dist2 < dist)
          {
            dist = dist2;
          }
        }
        totalDist += dist;
      }
      m_Controls->m_MeanErrorLCD->display(sqrt(totalDist/m_FixedLandmarks->GetSize()));
      m_Controls->m_MeanErrorLCD->show();
      m_Controls->m_MeanError->show();
    }
    else
    {
      m_Controls->m_MeanErrorLCD->hide();
      m_Controls->m_MeanError->hide();
    }
  }
  else
  {
    if (m_MovingLandmarks.IsNotNull() && m_FixedLandmarks.IsNotNull() && m_MovingLandmarks->GetSize() != 0 && m_FixedLandmarks->GetSize() != 0 && m_MovingLandmarks->GetSize() == m_FixedLandmarks->GetSize())
    {
      for(int pointId = 0; pointId < m_MovingLandmarks->GetSize(); ++pointId)
      {
        point1 = m_MovingLandmarks->GetPoint(pointId);
        point2 = m_FixedLandmarks->GetPoint(pointId);
        p1[0] = point1[0]; p1[1] = point1[1]; p1[2] = point1[2];
        p2[0] = point2[0]; p2[1] = point2[1]; p2[2] = point2[2];
        totalDist += vtkMath::Distance2BetweenPoints(p1, p2);
      }
      m_Controls->m_MeanErrorLCD->display(sqrt(totalDist/m_FixedLandmarks->GetSize()));
      m_Controls->m_MeanErrorLCD->show();
      m_Controls->m_MeanError->show();
    }
    else
    {
      m_Controls->m_MeanErrorLCD->hide();
      m_Controls->m_MeanError->hide();
    }
  }
}

void QmitkPointBasedRegistration::transformationChanged(int transform)
{
  m_Transformation = transform;
  this->checkLandmarkError();
}

// ICP with vtkLandmarkTransformation
void QmitkPointBasedRegistration::calculateLandmarkbasedWithICP(PointBasedRegistrationControlParameters* params)
{
  if(CheckCalculate())
  {
    mitk::Geometry3D::Pointer pointsGeometry = m_MovingLandmarks->GetGeometry(0);
    mitk::AffineGeometryFrame3D::Pointer movingLandmarksGeometry = m_MovingLandmarks->GetGeometry(0)->Clone();
    m_UndoPointsGeometryList.push_back(static_cast<mitk::Geometry3D *>(movingLandmarksGeometry.GetPointer()));
      
    mitk::BaseData::Pointer originalData = m_MovingNode->GetData();
    mitk::AffineGeometryFrame3D::Pointer originalDataGeometry = originalData->GetGeometry(0)->Clone();
    m_UndoGeometryList.push_back(static_cast<mitk::Geometry3D *>(originalDataGeometry.GetPointer()));

    vtkIdType pointId;
    vtkPoints* vPointsSource=vtkPoints::New();
    vtkCellArray* vCellsSource=vtkCellArray::New();
    for(pointId=0; pointId<m_MovingLandmarks->GetSize();++pointId)
    {
      mitk::Point3D pointSource=m_MovingLandmarks->GetPoint(pointId);
      vPointsSource->InsertNextPoint(pointSource[0],pointSource[1],pointSource[2]);
      vCellsSource->InsertNextCell(1, &pointId);
    } 
  
    vtkPoints* vPointsTarget=vtkPoints::New();
    vtkCellArray* vCellsTarget = vtkCellArray::New();
    for(pointId=0; pointId<m_FixedLandmarks->GetSize();++pointId)
    {
      mitk::Point3D pointTarget=m_FixedLandmarks->GetPoint(pointId);
      vPointsTarget->InsertNextPoint(pointTarget[0],pointTarget[1],pointTarget[2]);
      vCellsTarget->InsertNextCell(1, &pointId);
    }
  
    vtkPolyData* vPointSetSource=vtkPolyData::New();
    vtkPolyData* vPointSetTarget=vtkPolyData::New();

    vPointSetTarget->SetPoints(vPointsTarget);
    vPointSetTarget->SetVerts(vCellsTarget);
    vPointSetSource->SetPoints(vPointsSource);
    vPointSetSource->SetVerts(vCellsSource);

    vtkIterativeClosestPointTransform * icp=vtkIterativeClosestPointTransform::New();
    icp->SetCheckMeanDistance(1);
    icp->SetSource(vPointSetSource);
    icp->SetTarget(vPointSetTarget);
    icp->SetMaximumNumberOfIterations(50);
    icp->StartByMatchingCentroidsOn();
    vtkLandmarkTransform * transform=icp->GetLandmarkTransform();
    if((params->getTransformationType())==0)
    {
      transform->SetModeToRigidBody();
    }
    if((params->getTransformationType())==1)
    {
      transform->SetModeToSimilarity();
    }
    if((params->getTransformationType())==2)
    {
      transform->SetModeToAffine();
    }
    
    vtkMatrix4x4 * matrix=icp->GetMatrix();
  
    double determinant = fabs(matrix->Determinant());
    if((determinant < mitk::eps) || (determinant > 100) || (determinant < 0.01)
      || (determinant==itk::NumericTraits<double>::infinity()) 
      || (determinant==itk::NumericTraits<double>::quiet_NaN())
      || (determinant==itk::NumericTraits<double>::signaling_NaN())
      || (determinant==-itk::NumericTraits<double>::infinity()) 
      || (determinant==-itk::NumericTraits<double>::quiet_NaN())
      || (determinant==-itk::NumericTraits<double>::signaling_NaN())
      || (!(determinant <= 0) && !(determinant > 0)))
    {
      QMessageBox::critical( m_Controls, "PointBasedRegistration", "Suspicious determinant of matrix calculated by ICP.\n"
                                                             "Please select more points or other points!" );
      return;
    }

    pointsGeometry->Compose(matrix);
    m_MovingLandmarks->GetTimeSlicedGeometry()->UpdateInformation();
 
    mitk::BaseData::Pointer movingData = m_MovingNode->GetData();
    mitk::Geometry3D::Pointer movingGeometry = movingData->GetGeometry(0);
    movingGeometry->Compose(matrix);
    movingData->GetTimeSlicedGeometry()->UpdateInformation();
    m_Controls->setUndoTransformationEnabled(true);
    m_Controls->setRedoTransformationEnabled(false);
    m_RedoGeometryList.clear();
    m_RedoPointsGeometryList.clear();
    mitk::RenderingManager::GetInstance()->RequestUpdateAll();
    this->checkLandmarkError();
  }
}

// only vtkLandmarkTransformation
void QmitkPointBasedRegistration::calculateLandmarkbased(PointBasedRegistrationControlParameters* params)
{
  if(CheckCalculate())
  {
    mitk::Geometry3D::Pointer pointsGeometry = m_MovingLandmarks->GetGeometry(0);
    mitk::AffineGeometryFrame3D::Pointer movingLandmarksGeometry = m_MovingLandmarks->GetGeometry(0)->Clone();
    m_UndoPointsGeometryList.push_back(static_cast<mitk::Geometry3D *>(movingLandmarksGeometry.GetPointer()));
      
    mitk::BaseData::Pointer originalData = m_MovingNode->GetData();
    mitk::AffineGeometryFrame3D::Pointer originalDataGeometry = originalData->GetGeometry(0)->Clone();
    m_UndoGeometryList.push_back(static_cast<mitk::Geometry3D *>(originalDataGeometry.GetPointer()));

    vtkIdType pointId;
    vtkPoints* vPointsSource=vtkPoints::New();
    for(pointId = 0; pointId < m_MovingLandmarks->GetSize(); ++pointId)
    {
      mitk::Point3D sourcePoint = m_MovingLandmarks->GetPoint(pointId);
      vPointsSource->InsertNextPoint(sourcePoint[0],sourcePoint[1],sourcePoint[2]);
    }
    vtkPoints* vPointsTarget=vtkPoints::New();
    for(pointId=0; pointId<m_FixedLandmarks->GetSize();++pointId)
    {
      mitk::Point3D targetPoint=m_FixedLandmarks->GetPoint(pointId);
      vPointsTarget->InsertNextPoint(targetPoint[0],targetPoint[1],targetPoint[2]);
    }

    vtkLandmarkTransform * transform= vtkLandmarkTransform::New();
    transform->SetSourceLandmarks(vPointsSource);
    transform->SetTargetLandmarks(vPointsTarget);
    if((params->getTransformationType())==0)
    {
      transform->SetModeToRigidBody();
    }
    if((params->getTransformationType())==1)
    {
      transform->SetModeToRigidBody();
    }
    if((params->getTransformationType())==2)
    {
      transform->SetModeToRigidBody();
    }
    if((params->getTransformationType())==3)
    {
      transform->SetModeToRigidBody();
    }
    if((params->getTransformationType())==4)
    {
      transform->SetModeToSimilarity();
    }
    if((params->getTransformationType())==5)
    {
      transform->SetModeToAffine();
    }
    vtkMatrix4x4 * matrix=transform->GetMatrix();
    double determinant = fabs(matrix->Determinant());
    if((determinant < mitk::eps) || (determinant > 100) || (determinant < 0.01)
      || (determinant==itk::NumericTraits<double>::infinity()) 
      || (determinant==itk::NumericTraits<double>::quiet_NaN())
      || (determinant==itk::NumericTraits<double>::signaling_NaN())
      || (determinant==-itk::NumericTraits<double>::infinity()) 
      || (determinant==-itk::NumericTraits<double>::quiet_NaN())
      || (determinant==-itk::NumericTraits<double>::signaling_NaN())
      || (!(determinant <= 0) && !(determinant > 0)))
    {
      QMessageBox::critical( m_Controls, "PointBasedRegistration", "Suspicious determinant of matrix calculated.\n"
                                                             "Please select more points or other points!" );
      return;
    }
    pointsGeometry->Compose(matrix);
    m_MovingLandmarks->GetTimeSlicedGeometry()->UpdateInformation();

    mitk::BaseData::Pointer movingData = m_MovingNode->GetData();
    mitk::Geometry3D::Pointer movingGeometry = movingData->GetGeometry(0);
    movingGeometry->Compose(matrix);
    movingData->GetTimeSlicedGeometry()->UpdateInformation();
    m_Controls->setUndoTransformationEnabled(true);
    m_Controls->setRedoTransformationEnabled(false);
    m_RedoGeometryList.clear();
    m_RedoPointsGeometryList.clear();
    mitk::RenderingManager::GetInstance()->RequestUpdateAll();
    this->checkLandmarkError();
  }
}

void QmitkPointBasedRegistration::calculateLandmarkWarping(PointBasedRegistrationControlParameters* /*params*/)
{
  LandmarkWarping* registration = new LandmarkWarping();

  LandmarkWarping::FixedImageType::Pointer fixedImage = LandmarkWarping::FixedImageType::New();
  mitk::Image::Pointer fimage = dynamic_cast<mitk::Image*>(m_FixedNode->GetData());
  if (fimage.IsNotNull() && (fimage->GetDimension() == 2 || fimage->GetDimension() == 3))
  {
    if(fimage.IsNotNull())
    {
      mitk::CastToItkImage(fimage, fixedImage);
    }
    LandmarkWarping::MovingImageType::Pointer movingImage = LandmarkWarping::MovingImageType::New();
    mitk::Image::Pointer mimage = dynamic_cast<mitk::Image*>(m_MovingNode->GetData());
    if(mimage.IsNotNull())
    {
      mitk::CastToItkImage(mimage, movingImage);
    }
    registration->SetFixedImage(fixedImage);
    registration->SetMovingImage(movingImage);
    unsigned int pointId;
    mitk::Point3D sourcePoint, targetPoint;
    LandmarkWarping::LandmarkContainerType::Pointer fixedLandmarks = LandmarkWarping::LandmarkContainerType::New();
    LandmarkWarping::LandmarkPointType point;
    for(pointId = 0; pointId < (unsigned int)m_FixedLandmarks->GetSize(); ++pointId)
    {
      fimage->GetGeometry(0)->WorldToItkPhysicalPoint(m_FixedLandmarks->GetPoint(pointId), point);
      fixedLandmarks->InsertElement( pointId, point);
    }
    LandmarkWarping::LandmarkContainerType::Pointer movingLandmarks = LandmarkWarping::LandmarkContainerType::New();
    for(pointId = 0; pointId < (unsigned int)m_MovingLandmarks->GetSize(); ++pointId)
    {
      mitk::BaseData::Pointer fixedData = m_FixedNode->GetData();
      mitk::Geometry3D::Pointer fixedGeometry = fixedData->GetGeometry(0);
      fixedGeometry->WorldToItkPhysicalPoint(m_MovingLandmarks->GetPoint(pointId), point);
      movingLandmarks->InsertElement( pointId, point);
    }
    registration->SetLandmarks(fixedLandmarks.GetPointer(), movingLandmarks.GetPointer());
    LandmarkWarping::MovingImageType::Pointer output = registration->Register();
    if (output.IsNotNull()) 
    {
      mitk::Image::Pointer image = mitk::Image::New();
      mitk::CastToMitkImage(output, image);
      m_MovingNode->SetData(image);
      mitk::LevelWindowProperty::Pointer levWinProp = mitk::LevelWindowProperty::New();
      mitk::LevelWindow levelWindow;
      levelWindow.SetAuto( image );
      levWinProp->SetLevelWindow(levelWindow);
      m_MovingNode->GetPropertyList()->SetProperty("levelwindow",levWinProp);
      movingLandmarks = registration->GetTransformedTargetLandmarks();
      mitk::PointSet::PointDataIterator it;
      it = m_MovingLandmarks->GetPointSet()->GetPointData()->Begin();
      //increase the eventId to capsulate the coming operations
      mitk::OperationEvent::IncCurrObjectEventId();
      mitk::OperationEvent::ExecuteIncrement();
      for(pointId=0; pointId<movingLandmarks->Size();++pointId, ++it)
      {
        int position = it->Index();
        mitk::PointSet::PointType pt = m_MovingLandmarks->GetPoint(position);
        mitk::Point3D undoPoint = ( pt );
        point = movingLandmarks->GetElement(pointId);
        fimage->GetGeometry(0)->ItkPhysicalPointToWorld(point, pt);
        mitk::PointOperation* doOp = new mitk::PointOperation(mitk::OpMOVE, pt, position);
        //undo operation
        mitk::PointOperation* undoOp = new mitk::PointOperation(mitk::OpMOVE, undoPoint, position);
        mitk::OperationEvent* operationEvent = new mitk::OperationEvent(m_MovingLandmarks, doOp, undoOp, "Move point");
        mitk::UndoController::GetCurrentUndoModel()->SetOperationEvent(operationEvent);

        //execute the Operation
        m_MovingLandmarks->ExecuteOperation(doOp);
      }
      mitk::RenderingManager::GetInstance()->RequestUpdateAll();
      this->clearTransformationLists();
      this->checkLandmarkError();
    }
  }
}

void QmitkPointBasedRegistration::setInvisible(bool invisible)
{
  m_SetInvisible = invisible;
  if (m_SetInvisible)
  {
    this->TreeChanged();
  }
  else
  {
    mitk::DataTreeIteratorClone iter = this->GetDataTreeIterator();
    iter->GoToBegin();
    while ( !iter->IsAtEnd() )
    {
      if ( (iter->Get().IsNotNull()) && (iter->Get()->GetProperty("visible")) )
      {
        iter->Get()->SetVisibility(true);
      }
      ++iter;
    }
  }
  mitk::RenderingManager::GetInstance()->RequestUpdateAll();
}
