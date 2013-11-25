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

#include "QmitkDilateTool3DGUI.h"

#include "QmitkNewSegmentationDialog.h"
#include <QApplication.h>

MITK_TOOL_GUI_MACRO(SegmentationUI_EXPORT, QmitkDilateTool3DGUI, "")

QmitkDilateTool3DGUI::QmitkDilateTool3DGUI() : QmitkToolGUI()
{
  m_Controls.setupUi(this);
  m_Controls.m_InformationWidget->hide();
  m_Controls.m_AdvancedControlsWidget->hide();

  connect( m_Controls.m_pbRun, SIGNAL(clicked()), this, SLOT(OnRun()) );
  connect( m_Controls.m_pbAcceptPreview, SIGNAL(clicked()), this, SLOT(OnAcceptPreview()) );
  connect( m_Controls.m_pbDifference, SIGNAL(clicked()), this, SLOT(OnCalculateDifference()) );
  connect( m_Controls.m_pbNewLabel, SIGNAL(clicked()), this, SLOT(OnNewLabel()) );
  connect( m_Controls.m_sbKernelSize, SIGNAL(valueChanged(int)), this, SLOT(OnKernelSizeChanged(int)) );
  connect( m_Controls.m_cbShowInformation, SIGNAL(toggled(bool)), this, SLOT(OnShowInformation(bool)) );
  connect( m_Controls.m_cbShowAdvancedControls, SIGNAL(toggled(bool)), this, SLOT(OnShowAdvancedControls(bool)) );
  connect( this, SIGNAL(NewToolAssociated(mitk::Tool*)), this, SLOT(OnNewToolAssociated(mitk::Tool*)) );
}

QmitkDilateTool3DGUI::~QmitkDilateTool3DGUI()
{
  if (m_DilateTool3D.IsNotNull())
  {
    m_DilateTool3D->CurrentlyBusy -= mitk::MessageDelegate1<QmitkDilateTool3DGUI, bool>( this, &QmitkDilateTool3DGUI::BusyStateChanged );
  }
}

void QmitkDilateTool3DGUI::OnNewToolAssociated(mitk::Tool* tool)
{
  if (m_DilateTool3D.IsNotNull())
  {
    m_DilateTool3D->CurrentlyBusy -= mitk::MessageDelegate1<QmitkDilateTool3DGUI, bool>( this, &QmitkDilateTool3DGUI::BusyStateChanged );
  }

  m_DilateTool3D = dynamic_cast<mitk::DilateTool3D*>( tool );

  if (m_DilateTool3D.IsNotNull())
  {
    m_DilateTool3D->CurrentlyBusy += mitk::MessageDelegate1<QmitkDilateTool3DGUI, bool>( this, &QmitkDilateTool3DGUI::BusyStateChanged );
  }
}

void QmitkDilateTool3DGUI::OnRun()
{
  if (m_DilateTool3D.IsNotNull())
  {
    m_DilateTool3D->Run();
  }
}

void QmitkDilateTool3DGUI::OnAcceptPreview()
{
  if (m_DilateTool3D.IsNotNull())
  {
    m_DilateTool3D->AcceptPreview();
  }
}

void QmitkDilateTool3DGUI::OnCalculateDifference()
{
  if (m_DilateTool3D.IsNotNull())
  {
    m_DilateTool3D->CalculateDifference();
  }
}

void QmitkDilateTool3DGUI::OnNewLabel()
{
  if (m_DilateTool3D.IsNotNull())
  {
    QmitkNewSegmentationDialog* dialog = new QmitkNewSegmentationDialog( this );
//    dialog->SetSuggestionList( m_OrganColors );
    dialog->setWindowTitle("New Label");

    int dialogReturnValue = dialog->exec();

    if ( dialogReturnValue == QDialog::Rejected ) return;

    mitk::Color color = dialog->GetColor();

    m_DilateTool3D->CreateNewLabel(dialog->GetSegmentationName().toStdString(), color);
  }
}

void QmitkDilateTool3DGUI::BusyStateChanged(bool value)
{
  if (value)
    QApplication::setOverrideCursor( QCursor(Qt::BusyCursor) );
  else
    QApplication::restoreOverrideCursor();
}

void QmitkDilateTool3DGUI::OnKernelSizeChanged(int value)
{
  if (m_DilateTool3D.IsNotNull())
  {
    m_DilateTool3D->SetRadius(value);
  }
}

void QmitkDilateTool3DGUI::OnShowInformation(bool on)
{
  if (on)
    m_Controls.m_InformationWidget->show();
  else
    m_Controls.m_InformationWidget->hide();
}

void QmitkDilateTool3DGUI::OnShowAdvancedControls( bool on )
{
  if (on)
    m_Controls.m_AdvancedControlsWidget->show();
  else
    m_Controls.m_AdvancedControlsWidget->hide();
}
