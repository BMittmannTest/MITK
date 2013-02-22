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

#include "mitkEventObserver.h"

void mitk::EventObserver::Notify(InteractionEvent::Pointer /*interactionEvent*/, bool /*isHandled*/)
{
}

mitk::EventObserver::EventObserver()
{
}

mitk::EventObserver::~EventObserver()
{
}

bool mitk::EventObserver::FilterEvents(InteractionEvent* /*interactionEvent*/, DataNode* /*dataNode*/)
{
  return true;
}
