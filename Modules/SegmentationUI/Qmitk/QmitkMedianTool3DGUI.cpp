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

#include "QmitkMedianTool3DGUI.h"

#include <QApplication.h>

MITK_TOOL_GUI_MACRO(SegmentationUI_EXPORT, QmitkMedianTool3DGUI, "")

QmitkMedianTool3DGUI::QmitkMedianTool3DGUI() : QmitkToolGUI()
{
  m_Controls.setupUi(this);
  m_Controls.m_InformationWidget->hide();
  m_Controls.m_AdvancedControlsWidget->hide();

  connect( m_Controls.m_pbRun, SIGNAL(clicked()), this, SLOT(OnRun()) );
  connect( m_Controls.m_pbCancel, SIGNAL(clicked()), this, SLOT(OnCancel()) );
  connect( m_Controls.m_pbAcceptPreview, SIGNAL(clicked()), this, SLOT(OnAcceptPreview()) );
  connect( m_Controls.m_pbDifference, SIGNAL(clicked()), this, SLOT(OnCalculateDifference()) );
  connect( m_Controls.m_pbUnion, SIGNAL(clicked()), this, SLOT(OnCalculateUnion()) );
  connect( m_Controls.m_sbKernelSize, SIGNAL(valueChanged(int)), this, SLOT(OnKernelSizeChanged(int)) );
  connect( m_Controls.m_cbShowInformation, SIGNAL(toggled(bool)), this, SLOT(OnShowInformation(bool)) );
  connect( m_Controls.m_cbShowAdvancedControls, SIGNAL(toggled(bool)), this, SLOT(OnShowAdvancedControls(bool)) );
  connect( this, SIGNAL(NewToolAssociated(mitk::Tool*)), this, SLOT(OnNewToolAssociated(mitk::Tool*)) );
}

QmitkMedianTool3DGUI::~QmitkMedianTool3DGUI()
{
  if (m_MedianTool3D.IsNotNull())
  {
    m_MedianTool3D->CurrentlyBusy -= mitk::MessageDelegate1<QmitkMedianTool3DGUI, bool>( this, &QmitkMedianTool3DGUI::BusyStateChanged );
  }
}

void QmitkMedianTool3DGUI::OnNewToolAssociated(mitk::Tool* tool)
{
  if (m_MedianTool3D.IsNotNull())
  {
    m_MedianTool3D->CurrentlyBusy -= mitk::MessageDelegate1<QmitkMedianTool3DGUI, bool>( this, &QmitkMedianTool3DGUI::BusyStateChanged );
  }

  m_MedianTool3D = dynamic_cast<mitk::MedianTool3D*>( tool );

  if (m_MedianTool3D.IsNotNull())
  {
    m_MedianTool3D->CurrentlyBusy += mitk::MessageDelegate1<QmitkMedianTool3DGUI, bool>( this, &QmitkMedianTool3DGUI::BusyStateChanged );
  }
}

void QmitkMedianTool3DGUI::OnRun()
{
  if (m_MedianTool3D.IsNotNull())
  {
    m_MedianTool3D->Run();
  }
}

void QmitkMedianTool3DGUI::OnCancel()
{
  if (m_MedianTool3D.IsNotNull())
  {
    m_MedianTool3D->Cancel();
  }
}

void QmitkMedianTool3DGUI::OnAcceptPreview()
{
  if (m_MedianTool3D.IsNotNull())
  {
    m_MedianTool3D->AcceptPreview();
  }
}

void QmitkMedianTool3DGUI::OnCalculateDifference()
{
  if (m_MedianTool3D.IsNotNull())
  {
    m_MedianTool3D->CalculateDifference();
  }
}

void QmitkMedianTool3DGUI::OnCalculateUnion()
{
  if (m_MedianTool3D.IsNotNull())
  {
    m_MedianTool3D->CalculateUnion();
  }
}

void QmitkMedianTool3DGUI::OnKernelSizeChanged(int value)
{
  if (m_MedianTool3D.IsNotNull())
  {
    m_MedianTool3D->SetRadius(value);
  }
}

void QmitkMedianTool3DGUI::BusyStateChanged(bool value)
{
  if (value)
    QApplication::setOverrideCursor( QCursor(Qt::BusyCursor) );
  else
    QApplication::restoreOverrideCursor();
}

void QmitkMedianTool3DGUI::OnShowInformation( bool on )
{
  if (on)
    m_Controls.m_InformationWidget->show();
  else
    m_Controls.m_InformationWidget->hide();
}

void QmitkMedianTool3DGUI::OnShowAdvancedControls( bool on )
{
  if (on)
    m_Controls.m_AdvancedControlsWidget->show();
  else
    m_Controls.m_AdvancedControlsWidget->hide();
}
