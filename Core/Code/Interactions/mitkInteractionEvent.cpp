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

#include "mitkException.h"
#include "mitkInteractionEvent.h"

mitk::InteractionEvent::InteractionEvent(BaseRenderer* baseRenderer, std::string eventClass)
: m_Sender(baseRenderer)
, m_EventClass(eventClass)
{
}

void mitk::InteractionEvent::SetSender(mitk::BaseRenderer* sender)
{
  m_Sender = sender;
}

mitk::BaseRenderer* mitk::InteractionEvent::GetSender()
{
  return m_Sender;
}

bool mitk::InteractionEvent::MatchesTemplate(InteractionEvent::Pointer)
{
  return false;
}

mitk::InteractionEvent::~InteractionEvent()
{
}

bool mitk::InteractionEvent::IsSuperClassOf(InteractionEvent::Pointer baseClass)
{
  return (dynamic_cast<InteractionEvent*>(baseClass.GetPointer()) != NULL) ;
}

std::string mitk::InteractionEvent::GetEventClass() const
{
  return m_EventClass;
}
