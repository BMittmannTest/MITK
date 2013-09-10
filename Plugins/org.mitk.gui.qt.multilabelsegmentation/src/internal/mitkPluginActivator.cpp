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
#include "mitkPluginActivator.h"

#include "QmitkMultiLabelSegmentationView.h"
#include "QmitkMultiLabelSegmentationPreferencePage.h"

using namespace mitk;

void PluginActivator::start(ctkPluginContext *context)
{
  BERRY_REGISTER_EXTENSION_CLASS(QmitkMultiLabelSegmentationView, context)
  BERRY_REGISTER_EXTENSION_CLASS(QmitkMultiLabelSegmentationPreferencePage, context)
}

void PluginActivator::stop(ctkPluginContext *)
{
}

Q_EXPORT_PLUGIN2(org_mitk_gui_qt_multilabelsegmentation, mitk::PluginActivator)
