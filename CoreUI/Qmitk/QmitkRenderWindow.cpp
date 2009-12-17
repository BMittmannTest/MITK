/*=========================================================================

Program:   Medical Imaging & Interaction Toolkit
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


#include "QmitkRenderWindow.h"

#include <QCursor>
#include <QLayout>
#include <QMouseEvent>
#include <QWheelEvent>
#include <QKeyEvent>
#include <QResizeEvent>

#include "QmitkEventAdapter.h"

#include "mitkDisplayPositionEvent.h"
#include "mitkVtkLayerController.h"
#include "mitkRenderingManager.h"
#include "vtkRenderer.h"

#include "QmitkRenderWindowMenu.h"

QmitkRenderWindow::QmitkRenderWindow(QWidget *parent, QString name, mitk::VtkPropRenderer* renderer)
: QVTKWidget(parent), m_Renderer(renderer), m_ResendQtEvents(true), m_MenuWidgetActivated(false), m_MenuWidget(NULL)
{
  if(m_Renderer.IsNull())
  {
      m_Renderer = mitk::VtkPropRenderer::New( qPrintable(name), GetRenderWindow());
  }

  m_Renderer->InitRenderer(this->GetRenderWindow());

  mitk::BaseRenderer::AddInstance(GetRenderWindow(),m_Renderer);

  mitk::RenderingManager::GetInstance()->AddRenderWindow(GetRenderWindow());

  m_RenderProp = vtkMitkRenderProp::New();
  m_RenderProp->SetPropRenderer(m_Renderer);
  m_Renderer->GetVtkRenderer()->AddViewProp(m_RenderProp);

  if((this->GetRenderWindow()->GetSize()[0] > 10)
      && (this->GetRenderWindow()->GetSize()[1] > 10))
    m_Renderer->InitSize(this->GetRenderWindow()->GetSize()[0], this->GetRenderWindow()->GetSize()[1]);


  setFocusPolicy(Qt::StrongFocus);
  setMouseTracking(true);

  m_InResize = false;  

  //crate Renderwindow MenuBar for split, close Window or set new setting.
  m_MenuWidget = new QmitkRenderWindowMenu(this);
  
  //create Signal/Slot Connection
  connect( m_MenuWidget, SIGNAL( SignalChangeLayoutDesign(int) ), this, SLOT(OnChangeLayoutDesign(int)) );
  connect( m_MenuWidget, SIGNAL( ShowCrosshair(bool) ), this, SIGNAL( ShowCrosshair(bool)) );
  connect( m_MenuWidget, SIGNAL( ResetView() ), this, SIGNAL( ResetView()) );
  connect( m_MenuWidget, SIGNAL( ChangeCrosshairRotationMode(int) ), this, SIGNAL( ChangeCrosshairRotationMode(int)) );
  connect( m_MenuWidget, SIGNAL( SetCrosshairRotationLinked(bool) ), this, SIGNAL( SetCrosshairRotationLinked(bool)) );
}

QmitkRenderWindow::~QmitkRenderWindow()
{
  mitk::RenderingManager::GetInstance()->RemoveRenderWindow(GetRenderWindow());
  mitk::BaseRenderer::RemoveInstance(GetRenderWindow());
  m_Renderer->GetVtkRenderer()->RemoveViewProp(m_RenderProp);
  m_RenderProp->Delete();
}

mitk::VtkPropRenderer* QmitkRenderWindow::GetRenderer()
{
  return m_Renderer;
}

void QmitkRenderWindow::SetResendQtEvents(bool resend)
{
  m_ResendQtEvents = resend;
}

void QmitkRenderWindow::SetLayoutIndex( unsigned int layoutIndex )
{
  if( m_MenuWidget )
    m_MenuWidget->SetLayoutIndex(layoutIndex);
}

unsigned int QmitkRenderWindow::GetLayoutIndex()
{ 
  if( m_MenuWidget )
    return m_MenuWidget->GetLayoutIndex();
  else
    return NULL;
}

void QmitkRenderWindow::LayoutDesignListChanged( int layoutDesignIndex )
{ 
  if( m_MenuWidget ) 
    m_MenuWidget->UpdateLayoutDesignList( layoutDesignIndex );  
}

void QmitkRenderWindow::mousePressEvent(QMouseEvent *me)
{
  if (m_Renderer.IsNotNull())
  {
    mitk::MouseEvent event(QmitkEventAdapter::AdaptMouseEvent(m_Renderer, me));
    m_Renderer->MousePressEvent(&event);
  }
  QVTKWidget::mousePressEvent(me);

  if (m_ResendQtEvents) me->ignore();
}

void QmitkRenderWindow::mouseReleaseEvent(QMouseEvent *me)
{
  if (m_Renderer.IsNotNull())
  {
    mitk::MouseEvent event(QmitkEventAdapter::AdaptMouseEvent(m_Renderer, me));
    m_Renderer->MouseReleaseEvent(&event);
  }
  QVTKWidget::mouseReleaseEvent(me);

  if (m_ResendQtEvents) me->ignore();
}

void QmitkRenderWindow::mouseMoveEvent(QMouseEvent *me)
{
  if (m_Renderer.IsNotNull()) {
    mitk::MouseEvent event(QmitkEventAdapter::AdaptMouseEvent(m_Renderer, me));
    m_Renderer->MouseMoveEvent(&event);
  }
  QVTKWidget::mouseMoveEvent(me);

  if (m_ResendQtEvents) me->ignore();

  //Show/Hide Menu Widget
  if( m_MenuWidgetActivated )
  {
    //Show Menu Widget when mouse is inside of the define region of the top right corner
    if( m_MenuWidget->GetLayoutIndex() <= QmitkRenderWindowMenu::CORONAL 
      && me->pos().x() >= 0
      && me->pos().y() <= m_MenuWidget->height() + 20 )
    {
      m_MenuWidget->MoveWidgetToCorrectPos();
      m_MenuWidget->show();
      m_MenuWidget->update();
    }
    else if( m_MenuWidget->GetLayoutIndex() == QmitkRenderWindowMenu::THREE_D  
      && me->pos().x() >= this->width() - m_MenuWidget->width() - 20 
      && me->pos().y() <= m_MenuWidget->height() + 20 )
    {
      m_MenuWidget->MoveWidgetToCorrectPos();
      m_MenuWidget->show();
      m_MenuWidget->update();
    }
    //Hide Menu Widget when mouse is outside of the define region of the the right corner
    else if( !m_MenuWidget->GetSettingsMenuVisibilty() )
    {
      m_MenuWidget->hide();
    }    
  }
}

void QmitkRenderWindow::wheelEvent(QWheelEvent *we)
{
  //QVTKWidget::wheelEvent(we);

  if ( !GetSliceNavigationController()->GetSliceLocked() )
  {
    mitk::Stepper* stepper = GetSliceNavigationController()->GetSlice();

    if (stepper->GetSteps() <= 1)
    {
      stepper = GetSliceNavigationController()->GetTime();
    }

    if (we->orientation() * we->delta()  > 0)
    {
      stepper->Next();
    }
    else
    {
      stepper->Previous();
    }

    //also send to Renderer to send if to MITK interaction mechanism
    if (m_Renderer.IsNotNull()) 
    {
      mitk::WheelEvent event(QmitkEventAdapter::AdaptWheelEvent(m_Renderer, we));
      m_Renderer->WheelEvent(&event);
    }
  }

  QVTKWidget::wheelEvent(we);

  if (m_ResendQtEvents) 
    we->ignore();
}

void QmitkRenderWindow::keyPressEvent(QKeyEvent *ke)
{
  if (m_Renderer.IsNotNull())
  {
    QPoint cp = mapFromGlobal(QCursor::pos());
    mitk::KeyEvent mke(QmitkEventAdapter::AdaptKeyEvent(m_Renderer, ke, cp));
    m_Renderer->KeyPressEvent(&mke);
    ke->accept();
  }

  QVTKWidget::keyPressEvent(ke);

  if (m_ResendQtEvents) ke->ignore();
}

void QmitkRenderWindow::enterEvent( QEvent *e )
{
  //show Menu Widget
  //if( m_MenuWidget->isHidden() && m_MenuWidgetActivated )
  //{
  //  m_MenuWidget->MoveWidgetToCorrectPos();
  //  m_MenuWidget->show();
  //  m_MenuWidget->update();
  //}

  QVTKWidget::enterEvent(e);
}

void QmitkRenderWindow::leaveEvent( QEvent *e )
{
  //hide Menu Widget
  if( m_MenuWidget->isVisible() && !m_MenuWidget->GetSettingsMenuVisibilty() && m_MenuWidgetActivated )
    m_MenuWidget->hide();

  QVTKWidget::leaveEvent(e);
}

void QmitkRenderWindow::InitRenderer()
{
}

void QmitkRenderWindow::resizeEvent(QResizeEvent* event)
{
  QVTKWidget::resizeEvent(event);

  if(m_InResize) //@FIXME CRITICAL probably related to VtkSizeBug
    return;
  m_InResize = true;

  if(this->isVisible())
  {
    if(m_Renderer.IsNotNull())
    {
      m_Renderer->Resize(event->size().width(), event->size().height());
    }

    //this->update();
    //updateGL();
  }

  m_InResize = false;
}

mitk::SliceNavigationController * QmitkRenderWindow::GetSliceNavigationController()
{
  return mitk::BaseRenderer::GetInstance(this->GetRenderWindow())->GetSliceNavigationController();
}

mitk::CameraRotationController * QmitkRenderWindow::GetCameraRotationController()
{
  return mitk::BaseRenderer::GetInstance(this->GetRenderWindow())->GetCameraRotationController();
}

mitk::BaseController * QmitkRenderWindow::GetController()
{
  mitk::BaseRenderer * renderer = mitk::BaseRenderer::GetInstance(GetRenderWindow());
  switch ( renderer->GetMapperID() )
  {
    case mitk::BaseRenderer::Standard2D:
      return GetSliceNavigationController();

    case mitk::BaseRenderer::Standard3D:
      return GetCameraRotationController();

    default:
      return GetSliceNavigationController();
  }
}


void QmitkRenderWindow::OnChangeLayoutDesign( int layoutDesignIndex )
{
  
  emit SignalLayoutDesignChanged( layoutDesignIndex );
}

void QmitkRenderWindow::HideMenuWidget()
{
  if( m_MenuWidget )
    m_MenuWidget->hide();
}

void QmitkRenderWindow::OnUpdateCrosshairState( bool state )
{
  if( m_MenuWidget )
    m_MenuWidget->UpdateCrosshairState( state );
}

void QmitkRenderWindow::FullScreenMode(bool state)
{
  if( m_MenuWidget )
    m_MenuWidget->ChangeFullScreenMode( state );
}
  
QMenu* QmitkRenderWindow::GetCrossHairMenu()
{
  if (m_MenuWidget)
  {
    return m_MenuWidget->GetCrossHairMenu();
  }
}
  
void QmitkRenderWindow::SetCrossHairMenu(QMenu* menu)
{
  if (m_MenuWidget)
  {
    m_MenuWidget->SetCrossHairMenu( menu );
  }
}

