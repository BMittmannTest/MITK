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

#include "mitkDispatcher.h"
#include "mitkInteractionEvent.h"
#include "mitkInternalEvent.h"

// MicroServices
#include "mitkGetModuleContext.h"
#include "mitkModule.h"
#include "mitkModuleRegistry.h"

#include "mitkInteractionEventObserver.h"


mitk::Dispatcher::Dispatcher() :
    m_ProcessingMode(REGULAR)
{
  m_EventObserverTracker = new mitk::ServiceTracker<InteractionEventObserver*>(GetModuleContext());
  m_EventObserverTracker->Open();
}

void mitk::Dispatcher::AddDataInteractor(const DataNode* dataNode)
{
  RemoveDataInteractor(dataNode);
  RemoveOrphanedInteractors();
  DataInteractor::Pointer dataInteractor = dataNode->GetDataInteractor();
  if (dataInteractor.IsNotNull())
  {
    m_Interactors.push_back(dataInteractor);
  }
}

/*
 * Note: One DataInteractor can only have one DataNode and vice versa,
 * BUT the m_Interactors list may contain another DataInteractor that is still connected to this DataNode,
 * in this case we have to remove >1 DataInteractor. (Some special case of switching DataNodes between DataInteractors and registering a
 * DataNode to a DataStorage after assigning it to an DataInteractor)
 */

void mitk::Dispatcher::RemoveDataInteractor(const DataNode* dataNode)
{
  for (ListInteractorType::iterator it = m_Interactors.begin(); it != m_Interactors.end();)
  {
    if ((*it)->GetDataNode() == dataNode)
    {
      it = m_Interactors.erase(it);
    }
    else
    {
      ++it;
    }
  }
}

size_t mitk::Dispatcher::GetNumberOfInteractors()
{
  return m_Interactors.size();
}

mitk::Dispatcher::~Dispatcher()
{
  m_EventObserverTracker->Close();
  delete m_EventObserverTracker;
}

bool mitk::Dispatcher::ProcessEvent(InteractionEvent* event)
{
  InteractionEvent::Pointer p = event;
  //MITK_INFO << event->GetEventClass();
  bool eventIsHandled = false;
  /* Filter out and handle Internal Events separately */
  InternalEvent* internalEvent = dynamic_cast<InternalEvent*>(event);
  if (internalEvent != NULL)
  {
    eventIsHandled = HandleInternalEvent(internalEvent);
    // InternalEvents that are handled are not sent to the listeners
    if (eventIsHandled)
    {
      return true;
    }
  }
  switch (m_ProcessingMode)
  {
  case CONNECTEDMOUSEACTION:
    // finished connected mouse action
    if (p->GetEventClass() == "MouseReleaseEvent")
    {
      m_ProcessingMode = REGULAR;
      eventIsHandled = m_SelectedInteractor->HandleEvent(event, m_SelectedInteractor->GetDataNode());
    }
    // give event to selected interactor
    if (eventIsHandled == false)
    {
      eventIsHandled = m_SelectedInteractor->HandleEvent(event, m_SelectedInteractor->GetDataNode());
    }
    break;

  case GRABINPUT:
    eventIsHandled = m_SelectedInteractor->HandleEvent(event, m_SelectedInteractor->GetDataNode());
    SetEventProcessingMode(m_SelectedInteractor);
    break;

  case PREFERINPUT:
    if (m_SelectedInteractor->HandleEvent(event, m_SelectedInteractor->GetDataNode()) == true)
    {
      SetEventProcessingMode(m_SelectedInteractor);
      eventIsHandled = true;
    }
    break;

  case REGULAR:
    break;
  }
  // Standard behavior. Is executed in STANDARD mode  and PREFERINPUT mode, if preferred interactor rejects event.
  if (m_ProcessingMode == REGULAR || (m_ProcessingMode == PREFERINPUT && eventIsHandled == false))
  {

    m_Interactors.sort(cmp()); // sorts interactors by layer (descending);
    for (std::list<DataInteractor::Pointer>::iterator it = m_Interactors.begin(); it != m_Interactors.end(); ++it)
    {
      // explicit copy of pointer because HandleEvent function causes the m_Interactors list to be updated,
      // which in turn invalidates the iterator.
      DataInteractor::Pointer dataInteractor = *it;
      if (dataInteractor->HandleEvent(event, dataInteractor->GetDataNode()))
      { // if an event is handled several properties are checked, in order to determine the processing mode of the dispatcher
        SetEventProcessingMode(dataInteractor);
        if (p->GetEventClass() == "MousePressEvent" && m_ProcessingMode == REGULAR)
        {
          m_SelectedInteractor = dataInteractor;
          m_ProcessingMode = CONNECTEDMOUSEACTION;
        }
        eventIsHandled = true;
        break;
      }
    }
  }

  /* Notify InteractionEventObserver  */
  std::list<mitk::ServiceReference> listEventObserver;
  m_EventObserverTracker->GetServiceReferences(listEventObserver);
  for (std::list<mitk::ServiceReference>::iterator it = listEventObserver.begin(); it != listEventObserver.end(); ++it)
  {
    InteractionEventObserver* interactionEventObserver = m_EventObserverTracker->GetService(*it);
    if (interactionEventObserver != NULL)
    {
      if (interactionEventObserver->IsEnabled())
      {
        interactionEventObserver->Notify(event, eventIsHandled);
      }
    }
  }

  // Process event queue
  if (!m_QueuedEvents.empty())
  {
    InteractionEvent::Pointer e = m_QueuedEvents.front();
    m_QueuedEvents.pop_front();
    ProcessEvent(e);
  }
  return eventIsHandled;
}

/*
 * Checks if DataNodes associated with DataInteractors point back to them.
 * If not remove the DataInteractors. (This can happen when s.o. tries to set DataNodes to multiple DataInteractors)
 */
void mitk::Dispatcher::RemoveOrphanedInteractors()
{
  for (ListInteractorType::iterator it = m_Interactors.begin(); it != m_Interactors.end();)
  {
    DataNode::Pointer dn = (*it)->GetDataNode();
    if (dn.IsNull())
    {
      it = m_Interactors.erase(it);
    }
    else
    {
      DataInteractor::Pointer interactor = dn->GetDataInteractor();
      if (interactor != it->GetPointer())
      {
        it = m_Interactors.erase(it);
      }
      else
      {
        ++it;
      }
    }
  }
}

void mitk::Dispatcher::QueueEvent(InteractionEvent* event)
{
  m_QueuedEvents.push_back(event);
}

void mitk::Dispatcher::SetEventProcessingMode(DataInteractor::Pointer dataInteractor)
{
  m_ProcessingMode = dataInteractor->GetMode();
  if (dataInteractor->GetMode() != REGULAR)
  {
    m_SelectedInteractor = dataInteractor;
  }
}

bool mitk::Dispatcher::HandleInternalEvent(InternalEvent* internalEvent)
{
  if (internalEvent->GetSignalName() == IntDeactivateMe && internalEvent->GetTargetInteractor() != NULL)
  {
    internalEvent->GetTargetInteractor()->GetDataNode()->SetDataInteractor(NULL);
    internalEvent->GetTargetInteractor()->SetDataNode(NULL);

    mitk::RenderingManager::GetInstance()->RequestUpdateAll();
    return true;
  }
  return false;
}
