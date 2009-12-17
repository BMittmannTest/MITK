/*=========================================================================

Program:   Medical Imaging & Interaction Toolkit
Language:  C++
Date:      $Date: 2009-05-12 20:04:59 +0200 (Di, 12 Mai 2009) $
Version:   $Revision: 17180 $

Copyright (c) German Cancer Research Center, Division of Medical and
Biological Informatics. All rights reserved.
See MITKCopyright.txt or http://www.mitk.org/copyright.html for details.

This software is distributed WITHOUT ANY WARRANTY; without even
the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
PURPOSE.  See the above copyright notices for more information.

=========================================================================*/


#include "QmitkRenderWindowMenu.h"

#include <QHBoxLayout>
#include <QSpacerItem>
#include <QSize>
#include <QPainter>

#include<QGroupBox>
#include<QRadioButton>
#include<QAction>
#include<QLine>

#include "QmitkStdMultiWidget.h"

#include"iconClose.xpm"
#include"iconFullScreen.xpm"
#include"iconCrosshairMode.xpm"
#include"iconHoriSplit.xpm"
#include"iconSettings.xpm"
#include"iconVertiSplit.xpm"
#include"iconLeaveFullScreen.xpm"

#include <math.h>


QmitkRenderWindowMenu::QmitkRenderWindowMenu(QWidget *parent, Qt::WindowFlags f )
:QWidget(parent,f),
m_Settings(NULL),
m_Layout(0),
m_LayoutDesign(0),
m_OldLayoutDesign(0),
m_FullScreenMode(false),
m_CrosshairMenu(NULL)
{
  //Create Menu Widget
  this->CreateMenuWidget();
  this->setMinimumWidth(61); //DIRTY.. If you add or remove a button, you need to change the size.
  this->setMaximumWidth(61);
  this->setAutoFillBackground( true );
  this->hide();
  
  //this->setAttribute( Qt::WA_NoSystemBackground  );
  //this->setBackgroundRole( QPalette::Dark );
  //this->update();

  //SetOpacity  --  its just posible if the widget is a window. 
  //Windows indicates that the widget is a window, usually with a window system frame and a title bar, 
  //irrespective of whether the widget has a parent or not.
  //this->setWindowFlags( Qt::Window | Qt::FramelessWindowHint);
  //this->setAttribute(Qt::WA_TranslucentBackground);
  //this->setWindowOpacity(0.75);
  //this->setFocusPolicy( Qt::NoFocus );
}

QmitkRenderWindowMenu::~QmitkRenderWindowMenu()
{}

void QmitkRenderWindowMenu::CreateMenuWidget()
{
  QHBoxLayout* layout = new QHBoxLayout(this);
  layout->setAlignment( Qt::AlignRight );
  layout->setContentsMargins(1,1,1,1);

  QSize size( 13, 13 );

  // button for changing rotation mode
  m_CrosshairModeButton = new QPushButton();
  m_CrosshairModeButton->setMaximumSize(15, 15);
  m_CrosshairModeButton->setIconSize(size);
  m_CrosshairModeButton->setFlat( true );
  m_CrosshairModeButton->setMenu( this->GetCrossHairMenu() );
  m_CrosshairModeButton->setIcon( QIcon( iconCrosshairMode_xpm ) );
  layout->addWidget( m_CrosshairModeButton );

  //fullScreenButton
  m_FullScreenButton = new QPushButton();
  m_FullScreenButton->setMaximumSize(15, 15);
  m_FullScreenButton->setIconSize(size);
  m_FullScreenButton->setFlat( true );
  m_FullScreenButton->setIcon( QIcon( iconFullScreen_xpm ));
  layout->addWidget( m_FullScreenButton );

  //settingsButton
  m_SettingsButton = new QPushButton();
  m_SettingsButton->setMaximumSize(15, 15);
  m_SettingsButton->setIconSize(size);
  m_SettingsButton->setFlat( true );
  m_SettingsButton->setIcon( QIcon( iconSettings_xpm ));
  layout->addWidget( m_SettingsButton );

  //Create Connections -- coming soon?
  connect( m_FullScreenButton, SIGNAL( clicked(bool) ), this, SLOT(OnFullScreenButton(bool)) );
  connect( m_SettingsButton, SIGNAL( clicked(bool) ), this, SLOT(OnSettingsButton(bool)) );
}

QMenu* QmitkRenderWindowMenu::CreateCrosshairMenu()
{
  QMenu* crosshairModesMenu = new QMenu(this);
  
  m_ShowHideCrosshairVisibilityAction = new QAction(crosshairModesMenu);
  m_ShowHideCrosshairVisibilityAction->setText("Show crosshair");
  m_ShowHideCrosshairVisibilityAction->setCheckable(true);
  m_ShowHideCrosshairVisibilityAction->setChecked(true); // TODO observe current status
  crosshairModesMenu->addAction( m_ShowHideCrosshairVisibilityAction );

  QAction* resetViewAction = new QAction(crosshairModesMenu);
  resetViewAction->setText("Reset view");
  crosshairModesMenu->addAction( resetViewAction );

  QAction* rotationGroupSeparator = new QAction(crosshairModesMenu);
  rotationGroupSeparator->setSeparator(true);
  rotationGroupSeparator->setText("Rotation mode");
  crosshairModesMenu->addAction( rotationGroupSeparator );

  QActionGroup* rotationModeActionGroup = new QActionGroup(crosshairModesMenu);
  rotationModeActionGroup->setExclusive(true);

  QAction* noCrosshairRotation = new QAction(crosshairModesMenu);
  noCrosshairRotation->setActionGroup(rotationModeActionGroup);
  noCrosshairRotation->setText("No crosshair rotation");
  noCrosshairRotation->setCheckable(true);
  noCrosshairRotation->setChecked(true); // TODO observe current status
  noCrosshairRotation->setData( QmitkStdMultiWidget::PLANE_MODE_SLICING );
  crosshairModesMenu->addAction( noCrosshairRotation );

  QAction* singleCrosshairRotation = new QAction(crosshairModesMenu);
  singleCrosshairRotation->setActionGroup(rotationModeActionGroup);
  singleCrosshairRotation->setText("Crosshair rotation");
  singleCrosshairRotation->setCheckable(true);
  singleCrosshairRotation->setData( QmitkStdMultiWidget::PLANE_MODE_ROTATION  );
  crosshairModesMenu->addAction( singleCrosshairRotation );

  QAction* coupledCrosshairRotation = new QAction(crosshairModesMenu);
  coupledCrosshairRotation->setActionGroup(rotationModeActionGroup);
  coupledCrosshairRotation->setText("Coupled crosshair rotation");
  coupledCrosshairRotation->setCheckable(true);
  coupledCrosshairRotation->setData( QmitkStdMultiWidget::PLANE_MODE_ROTATION );
  crosshairModesMenu->addAction( coupledCrosshairRotation );

  QAction* swivelMode = new QAction(crosshairModesMenu);
  swivelMode->setActionGroup(rotationModeActionGroup);
  swivelMode->setText("Swivel mode");
  swivelMode->setCheckable(true);
  swivelMode->setData( QmitkStdMultiWidget::PLANE_MODE_SWIVEL );
  crosshairModesMenu->addAction( swivelMode );
  
  connect( m_ShowHideCrosshairVisibilityAction, SIGNAL(toggled(bool)), this, SIGNAL(ShowCrosshair(bool)));
  connect( resetViewAction, SIGNAL(triggered()), this, SIGNAL(ResetView()));
  connect( rotationModeActionGroup, SIGNAL(triggered(QAction*)), this, SLOT(OnCrosshairRotationModeSelected(QAction*)) ); 

  return crosshairModesMenu;
}

QMenu* QmitkRenderWindowMenu::GetCrossHairMenu()
{
  if (!m_CrosshairMenu)
  {
    m_CrosshairMenu = this->CreateCrosshairMenu();
  }

  return m_CrosshairMenu;
}

void QmitkRenderWindowMenu::SetCrossHairMenu(QMenu* menu)
{
  if (m_CrosshairMenu != menu)
  {
    m_CrosshairMenu = menu;
  
    if (m_CrosshairModeButton)
    {
      m_CrosshairModeButton->setMenu( GetCrossHairMenu() );
    }
  }
}

void QmitkRenderWindowMenu::CreateSettingsWidget()
{
  m_Settings = new QMenu(this);

  m_DefaultLayoutAction = new QAction( "standard layout", m_Settings );
  m_DefaultLayoutAction->setDisabled( true );

  m_2DImagesUpLayoutAction = new QAction( "2D images top, 3D bottom", m_Settings );
  m_2DImagesUpLayoutAction->setDisabled( false );

  m_2DImagesLeftLayoutAction = new QAction( "2D images left, 3D right", m_Settings );
  m_2DImagesLeftLayoutAction->setDisabled( false );

  m_Big3DLayoutAction = new QAction( "Big 3D", m_Settings );
  m_Big3DLayoutAction->setDisabled( false );

  m_Widget1LayoutAction = new QAction( "Transversal plane", m_Settings );
  m_Widget1LayoutAction->setDisabled( false );

  m_Widget2LayoutAction = new QAction( "Sagittal plane", m_Settings );
  m_Widget2LayoutAction->setDisabled( false );

  m_Widget3LayoutAction = new QAction( "Coronal plane", m_Settings );
  m_Widget3LayoutAction->setDisabled( false );

  m_RowWidget3And4LayoutAction = new QAction( "Coronal top, 3D bottom", m_Settings );
  m_RowWidget3And4LayoutAction->setDisabled( false );

  m_ColumnWidget3And4LayoutAction = new QAction( "Coronal left, 3D right", m_Settings );
  m_ColumnWidget3And4LayoutAction->setDisabled( false );

  m_SmallUpperWidget2Big3and4LayoutAction = new QAction( "Sagittal top, Coronal n 3D bottom", m_Settings );
  m_SmallUpperWidget2Big3and4LayoutAction->setDisabled( false );

  m_2x2Dand3DWidgetLayoutAction = new QAction( "Transversal n Sagittal left, 3D right", m_Settings );
  m_2x2Dand3DWidgetLayoutAction->setDisabled( false );

  m_Left2Dand3DRight2DLayoutAction = new QAction( "Transversal n 3D left, Sagittal right", m_Settings );
  m_Left2Dand3DRight2DLayoutAction->setDisabled( false );

  m_Settings->addAction(m_DefaultLayoutAction);
  m_Settings->addAction(m_2DImagesUpLayoutAction);
  m_Settings->addAction(m_2DImagesLeftLayoutAction);
  m_Settings->addAction(m_Big3DLayoutAction);
  m_Settings->addAction(m_Widget1LayoutAction);
  m_Settings->addAction(m_Widget2LayoutAction);
  m_Settings->addAction(m_Widget3LayoutAction);
  m_Settings->addAction(m_RowWidget3And4LayoutAction);
  m_Settings->addAction(m_ColumnWidget3And4LayoutAction);
  m_Settings->addAction(m_SmallUpperWidget2Big3and4LayoutAction);
  m_Settings->addAction(m_2x2Dand3DWidgetLayoutAction);
  m_Settings->addAction(m_Left2Dand3DRight2DLayoutAction);

  m_Settings->setVisible( false );

  connect( m_DefaultLayoutAction, SIGNAL( triggered(bool) ), this, SLOT(OnChangeLayoutToDefault(bool)) );
  connect( m_2DImagesUpLayoutAction, SIGNAL( triggered(bool) ), this, SLOT(OnChangeLayoutTo2DImagesUp(bool)) );
  connect( m_2DImagesLeftLayoutAction, SIGNAL( triggered(bool) ), this, SLOT(OnChangeLayoutTo2DImagesLeft(bool)) );
  connect( m_Big3DLayoutAction, SIGNAL( triggered(bool) ), this, SLOT(OnChangeLayoutToBig3D(bool)) );
  connect( m_Widget1LayoutAction, SIGNAL( triggered(bool) ), this, SLOT(OnChangeLayoutToWidget1(bool)) );
  connect( m_Widget2LayoutAction, SIGNAL( triggered(bool) ), this, SLOT(OnChangeLayoutToWidget2(bool)) );
  connect( m_Widget3LayoutAction, SIGNAL( triggered(bool) ), this, SLOT(OnChangeLayoutToWidget3(bool)) );
  connect( m_RowWidget3And4LayoutAction, SIGNAL( triggered(bool) ), this, SLOT(OnChangeLayoutToRowWidget3And4(bool)) );
  connect( m_ColumnWidget3And4LayoutAction, SIGNAL( triggered(bool) ), this, SLOT(OnChangeLayoutToColumnWidget3And4(bool)) );
  connect( m_SmallUpperWidget2Big3and4LayoutAction, SIGNAL( triggered(bool) ), this, SLOT(OnChangeLayoutToSmallUpperWidget2Big3and4(bool)) );
  connect( m_2x2Dand3DWidgetLayoutAction, SIGNAL( triggered(bool) ), this, SLOT(OnChangeLayoutTo2x2Dand3DWidget(bool)) );
  connect( m_Left2Dand3DRight2DLayoutAction, SIGNAL( triggered(bool) ), this, SLOT(OnChangeLayoutToLeft2Dand3DRight2D(bool)) );

}

void QmitkRenderWindowMenu::paintEvent( QPaintEvent*  /*e*/ )
{
  QPainter painter(this);
  QColor semiTransparentColor = Qt::black;
  semiTransparentColor.setAlpha(255);
  painter.fillRect(rect(), semiTransparentColor);
}

void QmitkRenderWindowMenu::SetLayoutIndex( unsigned int layoutIndex )
{
  m_Layout = layoutIndex;
}

void QmitkRenderWindowMenu::ChangeFullScreenMode( bool state )
{
  this->OnFullScreenButton( state );
}
/// \brief
void QmitkRenderWindowMenu::OnFullScreenButton( bool  /*checked*/ )
{
  if( !m_FullScreenMode )
  {
    m_FullScreenMode = true;
    m_OldLayoutDesign = m_LayoutDesign;

    switch( m_Layout )
    {
    case TRANSVERSAL:
      {
        emit SignalChangeLayoutDesign( LAYOUT_TRANSVERSAL );
        break;
      }

    case SAGITTAL:
      {
        emit SignalChangeLayoutDesign( LAYOUT_SAGITTAL );
        break;
      }
    case CORONAL:
      {
        emit SignalChangeLayoutDesign( LAYOUT_CORONAL );
        break;
      }
    case THREE_D:
      {
        emit SignalChangeLayoutDesign( LAYOUT_BIG3D );
        break;
      }
    }

    //Move Widget and show again
    this->MoveWidgetToCorrectPos();

    //change icon
    this->ChangeFullScreenIcon();

  }
  else
  {
    m_FullScreenMode = false;
    emit SignalChangeLayoutDesign( m_OldLayoutDesign );

    //Move Widget and show again
    this->MoveWidgetToCorrectPos();

    //change icon
    this->ChangeFullScreenIcon();
  }
}


/// \brief
void QmitkRenderWindowMenu::OnSettingsButton( bool  /*checked*/ )
{
  if( m_Settings == NULL )
    this->CreateSettingsWidget();

  QPoint point = this->mapToGlobal( m_SettingsButton->geometry().topLeft() );
  m_Settings->setVisible( true );
  m_Settings->exec( point );
}

void QmitkRenderWindowMenu::OnChangeLayoutTo2DImagesUp(bool)
{
  //set Full Screen Mode to false, if Layout Design was changed by the LayoutDesign_List
  m_FullScreenMode = false;
  this->ChangeFullScreenIcon();

  m_LayoutDesign = LAYOUT_2DIMAGEUP;
  emit SignalChangeLayoutDesign( LAYOUT_2DIMAGEUP );
}
void QmitkRenderWindowMenu::OnChangeLayoutTo2DImagesLeft(bool)
{
  //set Full Screen Mode to false, if Layout Design was changed by the LayoutDesign_List
  m_FullScreenMode = false;
  this->ChangeFullScreenIcon();

  m_LayoutDesign = LAYOUT_2DIMAGELEFT;
  emit SignalChangeLayoutDesign( LAYOUT_2DIMAGELEFT );
}
void QmitkRenderWindowMenu::OnChangeLayoutToDefault(bool)
{
  //set Full Screen Mode to false, if Layout Design was changed by the LayoutDesign_List
  m_FullScreenMode = false;
  this->ChangeFullScreenIcon();

  m_LayoutDesign = LAYOUT_DEFAULT;
  emit SignalChangeLayoutDesign( LAYOUT_DEFAULT );
}
void QmitkRenderWindowMenu::OnChangeLayoutToBig3D(bool)
{ 
  //set Full Screen Mode to false, if Layout Design was changed by the LayoutDesign_List
  m_FullScreenMode = false;
  this->ChangeFullScreenIcon();

  m_LayoutDesign = LAYOUT_BIG3D;
  emit SignalChangeLayoutDesign( LAYOUT_BIG3D );
}
void QmitkRenderWindowMenu::OnChangeLayoutToWidget1(bool)
{
  //set Full Screen Mode to false, if Layout Design was changed by the LayoutDesign_List
  m_FullScreenMode = false;
  this->ChangeFullScreenIcon();

  m_LayoutDesign = LAYOUT_TRANSVERSAL;
  emit SignalChangeLayoutDesign( LAYOUT_TRANSVERSAL );
}
void QmitkRenderWindowMenu::OnChangeLayoutToWidget2(bool)
{
  //set Full Screen Mode to false, if Layout Design was changed by the LayoutDesign_List
  m_FullScreenMode = false;
  this->ChangeFullScreenIcon();

  m_LayoutDesign = LAYOUT_SAGITTAL;
  emit SignalChangeLayoutDesign( LAYOUT_SAGITTAL );
}
void QmitkRenderWindowMenu::OnChangeLayoutToWidget3(bool)
{
  //set Full Screen Mode to false, if Layout Design was changed by the LayoutDesign_List
  m_FullScreenMode = false;
  this->ChangeFullScreenIcon();

  m_LayoutDesign = LAYOUT_CORONAL;
  emit SignalChangeLayoutDesign( LAYOUT_CORONAL );
}
void QmitkRenderWindowMenu::OnChangeLayoutToRowWidget3And4(bool)
{
  //set Full Screen Mode to false, if Layout Design was changed by the LayoutDesign_List
  m_FullScreenMode = false;
  this->ChangeFullScreenIcon();

  m_LayoutDesign = LAYOUT_ROWWIDGET3AND4;
  emit SignalChangeLayoutDesign( LAYOUT_ROWWIDGET3AND4 );
}
void QmitkRenderWindowMenu::OnChangeLayoutToColumnWidget3And4(bool)
{
  //set Full Screen Mode to false, if Layout Design was changed by the LayoutDesign_List
  m_FullScreenMode = false;
  this->ChangeFullScreenIcon();

  m_LayoutDesign = LAYOUT_COLUMNWIDGET3AND4;
  emit SignalChangeLayoutDesign( LAYOUT_COLUMNWIDGET3AND4 );
}

void QmitkRenderWindowMenu::OnChangeLayoutToSmallUpperWidget2Big3and4(bool)
{
  //set Full Screen Mode to false, if Layout Design was changed by the LayoutDesign_List
  m_FullScreenMode = false;
  this->ChangeFullScreenIcon();

  m_LayoutDesign = LAYOUT_SMALLUPPERWIDGET2BIGAND4;
  emit SignalChangeLayoutDesign( LAYOUT_SMALLUPPERWIDGET2BIGAND4 );
}
void QmitkRenderWindowMenu::OnChangeLayoutTo2x2Dand3DWidget(bool)
{
  //set Full Screen Mode to false, if Layout Design was changed by the LayoutDesign_List
  m_FullScreenMode = false;
  this->ChangeFullScreenIcon();

  m_LayoutDesign = LAYOUT_2X2DAND3DWIDGET;
  emit SignalChangeLayoutDesign( LAYOUT_2X2DAND3DWIDGET );
}
void QmitkRenderWindowMenu::OnChangeLayoutToLeft2Dand3DRight2D(bool)
{
  //set Full Screen Mode to false, if Layout Design was changed by the LayoutDesign_List
  m_FullScreenMode = false;
  this->ChangeFullScreenIcon();

  m_LayoutDesign = LAYOUT_LEFT2DAND3DRIGHT2D;
  emit SignalChangeLayoutDesign( LAYOUT_LEFT2DAND3DRIGHT2D );
}

void QmitkRenderWindowMenu::UpdateLayoutDesignList( int layoutDesignIndex )
{
  m_LayoutDesign = layoutDesignIndex;

  if( m_Settings == NULL )
    this->CreateSettingsWidget();

  switch( m_LayoutDesign )
  {
  case LAYOUT_DEFAULT:
    {
      m_DefaultLayoutAction->setEnabled(false);
      m_2DImagesUpLayoutAction->setEnabled(true);
      m_2DImagesLeftLayoutAction->setEnabled(true);
      m_Big3DLayoutAction->setEnabled(true);
      m_Widget1LayoutAction->setEnabled(true);
      m_Widget2LayoutAction->setEnabled(true);
      m_Widget3LayoutAction->setEnabled(true);
      m_RowWidget3And4LayoutAction->setEnabled(true);
      m_ColumnWidget3And4LayoutAction->setEnabled(true);
      m_SmallUpperWidget2Big3and4LayoutAction->setEnabled(true);
      m_2x2Dand3DWidgetLayoutAction->setEnabled(true);
      m_Left2Dand3DRight2DLayoutAction->setEnabled(true);
      break;
    }

  case LAYOUT_2DIMAGEUP:
    {
      m_DefaultLayoutAction->setEnabled(true);
      m_2DImagesUpLayoutAction->setEnabled(false);
      m_2DImagesLeftLayoutAction->setEnabled(true);
      m_Big3DLayoutAction->setEnabled(true);
      m_Widget1LayoutAction->setEnabled(true);
      m_Widget2LayoutAction->setEnabled(true);
      m_Widget3LayoutAction->setEnabled(true);
      m_RowWidget3And4LayoutAction->setEnabled(true);
      m_ColumnWidget3And4LayoutAction->setEnabled(true);
      m_SmallUpperWidget2Big3and4LayoutAction->setEnabled(true);
      m_2x2Dand3DWidgetLayoutAction->setEnabled(true);
      m_Left2Dand3DRight2DLayoutAction->setEnabled(true);
      break;
    }
  case LAYOUT_2DIMAGELEFT:
    {
      m_DefaultLayoutAction->setEnabled(true);
      m_2DImagesUpLayoutAction->setEnabled(true);
      m_2DImagesLeftLayoutAction->setEnabled(false);
      m_Big3DLayoutAction->setEnabled(true);
      m_Widget1LayoutAction->setEnabled(true);
      m_Widget2LayoutAction->setEnabled(true);
      m_Widget3LayoutAction->setEnabled(true);
      m_RowWidget3And4LayoutAction->setEnabled(true);
      m_ColumnWidget3And4LayoutAction->setEnabled(true);
      m_SmallUpperWidget2Big3and4LayoutAction->setEnabled(true);
      m_2x2Dand3DWidgetLayoutAction->setEnabled(true);
      m_Left2Dand3DRight2DLayoutAction->setEnabled(true);
      break;
    }
  case LAYOUT_BIG3D:
    {
      m_DefaultLayoutAction->setEnabled(true);
      m_2DImagesUpLayoutAction->setEnabled(true);
      m_2DImagesLeftLayoutAction->setEnabled(true);
      m_Big3DLayoutAction->setEnabled(false);
      m_Widget1LayoutAction->setEnabled(true);
      m_Widget2LayoutAction->setEnabled(true);
      m_Widget3LayoutAction->setEnabled(true);
      m_RowWidget3And4LayoutAction->setEnabled(true);
      m_ColumnWidget3And4LayoutAction->setEnabled(true);
      m_SmallUpperWidget2Big3and4LayoutAction->setEnabled(true);
      m_2x2Dand3DWidgetLayoutAction->setEnabled(true);
      m_Left2Dand3DRight2DLayoutAction->setEnabled(true);
      break;
    }
  case LAYOUT_TRANSVERSAL:
    {
      m_DefaultLayoutAction->setEnabled(true);
      m_2DImagesUpLayoutAction->setEnabled(true);
      m_2DImagesLeftLayoutAction->setEnabled(true);
      m_Big3DLayoutAction->setEnabled(true);
      m_Widget1LayoutAction->setEnabled(false);
      m_Widget2LayoutAction->setEnabled(true);
      m_Widget3LayoutAction->setEnabled(true);
      m_RowWidget3And4LayoutAction->setEnabled(true);
      m_ColumnWidget3And4LayoutAction->setEnabled(true);
      m_SmallUpperWidget2Big3and4LayoutAction->setEnabled(true);
      m_2x2Dand3DWidgetLayoutAction->setEnabled(true);
      m_Left2Dand3DRight2DLayoutAction->setEnabled(true);
      break;
    }
  case LAYOUT_SAGITTAL:
    {
      m_DefaultLayoutAction->setEnabled(true);
      m_2DImagesUpLayoutAction->setEnabled(true);
      m_2DImagesLeftLayoutAction->setEnabled(true);
      m_Big3DLayoutAction->setEnabled(true);
      m_Widget1LayoutAction->setEnabled(true);
      m_Widget2LayoutAction->setEnabled(false);
      m_Widget3LayoutAction->setEnabled(true);
      m_RowWidget3And4LayoutAction->setEnabled(true);
      m_ColumnWidget3And4LayoutAction->setEnabled(true);
      m_SmallUpperWidget2Big3and4LayoutAction->setEnabled(true);
      m_2x2Dand3DWidgetLayoutAction->setEnabled(true);
      m_Left2Dand3DRight2DLayoutAction->setEnabled(true);
      break;
    }
  case LAYOUT_CORONAL:
    {
      m_DefaultLayoutAction->setEnabled(true);
      m_2DImagesUpLayoutAction->setEnabled(true);
      m_2DImagesLeftLayoutAction->setEnabled(true);
      m_Big3DLayoutAction->setEnabled(true);
      m_Widget1LayoutAction->setEnabled(true);
      m_Widget2LayoutAction->setEnabled(true);
      m_Widget3LayoutAction->setEnabled(false);
      m_RowWidget3And4LayoutAction->setEnabled(true);
      m_ColumnWidget3And4LayoutAction->setEnabled(true);
      m_SmallUpperWidget2Big3and4LayoutAction->setEnabled(true);
      m_2x2Dand3DWidgetLayoutAction->setEnabled(true);
      m_Left2Dand3DRight2DLayoutAction->setEnabled(true);
      break;
    }
  case LAYOUT_2X2DAND3DWIDGET:
    { 
      m_DefaultLayoutAction->setEnabled(true);
      m_2DImagesUpLayoutAction->setEnabled(true);
      m_2DImagesLeftLayoutAction->setEnabled(true);
      m_Big3DLayoutAction->setEnabled(true);
      m_Widget1LayoutAction->setEnabled(true);
      m_Widget2LayoutAction->setEnabled(true);
      m_Widget3LayoutAction->setEnabled(true);
      m_RowWidget3And4LayoutAction->setEnabled(true);
      m_ColumnWidget3And4LayoutAction->setEnabled(true);
      m_SmallUpperWidget2Big3and4LayoutAction->setEnabled(true);
      m_2x2Dand3DWidgetLayoutAction->setEnabled(false);
      m_Left2Dand3DRight2DLayoutAction->setEnabled(true);
      break;
    }
  case LAYOUT_ROWWIDGET3AND4:
    {
      m_DefaultLayoutAction->setEnabled(true);
      m_2DImagesUpLayoutAction->setEnabled(true);
      m_2DImagesLeftLayoutAction->setEnabled(true);
      m_Big3DLayoutAction->setEnabled(true);
      m_Widget1LayoutAction->setEnabled(true);
      m_Widget2LayoutAction->setEnabled(true);
      m_Widget3LayoutAction->setEnabled(true);
      m_RowWidget3And4LayoutAction->setEnabled(false);
      m_ColumnWidget3And4LayoutAction->setEnabled(true);
      m_SmallUpperWidget2Big3and4LayoutAction->setEnabled(true);
      m_2x2Dand3DWidgetLayoutAction->setEnabled(true);
      m_Left2Dand3DRight2DLayoutAction->setEnabled(true);
      break;
    }
  case LAYOUT_COLUMNWIDGET3AND4:
    {
      m_DefaultLayoutAction->setEnabled(true);
      m_2DImagesUpLayoutAction->setEnabled(true);
      m_2DImagesLeftLayoutAction->setEnabled(true);
      m_Big3DLayoutAction->setEnabled(true);
      m_Widget1LayoutAction->setEnabled(true);
      m_Widget2LayoutAction->setEnabled(true);
      m_Widget3LayoutAction->setEnabled(true);
      m_RowWidget3And4LayoutAction->setEnabled(true);
      m_ColumnWidget3And4LayoutAction->setEnabled(false);
      m_SmallUpperWidget2Big3and4LayoutAction->setEnabled(true);
      m_2x2Dand3DWidgetLayoutAction->setEnabled(true);
      m_Left2Dand3DRight2DLayoutAction->setEnabled(true);
      break;
    }
  case LAYOUT_SMALLUPPERWIDGET2BIGAND4:
    {
      m_DefaultLayoutAction->setEnabled(true);
      m_2DImagesUpLayoutAction->setEnabled(true);
      m_2DImagesLeftLayoutAction->setEnabled(true);
      m_Big3DLayoutAction->setEnabled(true);
      m_Widget1LayoutAction->setEnabled(true);
      m_Widget2LayoutAction->setEnabled(true);
      m_Widget3LayoutAction->setEnabled(true);
      m_RowWidget3And4LayoutAction->setEnabled(true);
      m_ColumnWidget3And4LayoutAction->setEnabled(true);
      m_SmallUpperWidget2Big3and4LayoutAction->setEnabled(false);
      m_2x2Dand3DWidgetLayoutAction->setEnabled(true);
      m_Left2Dand3DRight2DLayoutAction->setEnabled(true);
      break;
    }
  case LAYOUT_LEFT2DAND3DRIGHT2D:
    {
      m_DefaultLayoutAction->setEnabled(true);
      m_2DImagesUpLayoutAction->setEnabled(true);
      m_2DImagesLeftLayoutAction->setEnabled(true);
      m_Big3DLayoutAction->setEnabled(true);
      m_Widget1LayoutAction->setEnabled(true);
      m_Widget2LayoutAction->setEnabled(true);
      m_Widget3LayoutAction->setEnabled(true);
      m_RowWidget3And4LayoutAction->setEnabled(true);
      m_ColumnWidget3And4LayoutAction->setEnabled(true);
      m_SmallUpperWidget2Big3and4LayoutAction->setEnabled(true);
      m_2x2Dand3DWidgetLayoutAction->setEnabled(true);
      m_Left2Dand3DRight2DLayoutAction->setEnabled(false);
      break;
    }
  }
}

void QmitkRenderWindowMenu::MoveWidgetToCorrectPos()
{
  int moveX= floor( double(this->parentWidget()->width() - this->width() - 4.0) );
  this->move( moveX, 3 );
  this->show();
}

void QmitkRenderWindowMenu::ChangeFullScreenIcon()
{

 if( m_FullScreenMode )
 {   
   const QIcon icon( iconLeaveFullScreen_xpm );
   m_FullScreenButton->setIcon(icon);
 }
  else
  {
    const QIcon icon( iconFullScreen_xpm );
    m_FullScreenButton->setIcon(icon);
 }
}

void QmitkRenderWindowMenu::OnCrosshairRotationModeSelected(QAction* action)
{
  emit ChangeCrosshairRotationMode( action->data().toInt() );
  emit SetCrosshairRotationLinked( action->text().contains("Coupled") );
}

void QmitkRenderWindowMenu::UpdateCrosshairState( bool state )
{
  if( m_ShowHideCrosshairVisibilityAction->isChecked() != state )
    m_ShowHideCrosshairVisibilityAction->setChecked(state);
}

