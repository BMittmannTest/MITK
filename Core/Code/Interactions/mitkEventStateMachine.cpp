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

#include "mitkEventStateMachine.h"
#include "mitkStateMachineContainer.h"
#include "mitkInteractionEvent.h"
#include "mitkStateMachineAction.h"
#include "mitkStateMachineTransition.h"
#include "mitkStateMachineState.h"
// us
#include "mitkModule.h"
#include "mitkModuleResource.h"
#include "mitkModuleResourceStream.h"
#include "mitkModuleRegistry.h"

mitk::EventStateMachine::EventStateMachine() :
    m_StateMachineContainer(NULL), m_CurrentState(NULL)
{
}

bool mitk::EventStateMachine::LoadStateMachine(const std::string& filename, const Module* module)
{
  if (m_StateMachineContainer != NULL)
  {
    m_StateMachineContainer->Delete();
  }
  m_StateMachineContainer = StateMachineContainer::New();

  if (m_StateMachineContainer->LoadBehavior(filename, module))
  {
    m_CurrentState = m_StateMachineContainer->GetStartState();

    // clear actions map ,and connect all actions as declared in sub-class
    for(std::map<std::string, TActionFunctor*>::iterator i = m_ActionFunctionsMap.begin();
        i != m_ActionFunctionsMap.end(); ++i)
    {
      delete i->second;
    }
    m_ActionFunctionsMap.clear();
    for(ActionDelegatesMapType::iterator i = m_ActionDelegatesMap.begin();
        i != m_ActionDelegatesMap.end(); ++i)
    {
      delete i->second;
    }
    m_ActionDelegatesMap.clear();

    ConnectActionsAndFunctions();
    return true;
  }
  else
  {
    MITK_WARN<< "Unable to load StateMachine from file: " << filename;
    return false;
  }
}

mitk::EventStateMachine::~EventStateMachine()
{
  if (m_StateMachineContainer != NULL)
  {
    m_StateMachineContainer->Delete();
  }
}

void mitk::EventStateMachine::AddActionFunction(const std::string& action, mitk::TActionFunctor* functor)
{
  if (!functor)
    return;
// make sure double calls for same action won't cause memory leaks
  delete m_ActionFunctionsMap[action];
  ActionDelegatesMapType::iterator i = m_ActionDelegatesMap.find(action);
  if (i != m_ActionDelegatesMap.end())
  {
    delete i->second;
    m_ActionDelegatesMap.erase(i);
  }
  m_ActionFunctionsMap[action] = functor;
}

void mitk::EventStateMachine::AddActionFunction(const std::string& action, const ActionFunctionDelegate& delegate)
{
  std::map<std::string, TActionFunctor*>::iterator i = m_ActionFunctionsMap.find(action);
  if (i != m_ActionFunctionsMap.end())
  {
    delete i->second;
    m_ActionFunctionsMap.erase(i);
  }

  delete m_ActionDelegatesMap[action];
  m_ActionDelegatesMap[action] = delegate.Clone();
}

bool mitk::EventStateMachine::HandleEvent(InteractionEvent* event, DataNode* dataNode)
{
  if (!FilterEvents(event, dataNode))
  {
    return false;
  }
  // check if the current state holds a transition that works with the given event.
  StateMachineTransition::Pointer transition = m_CurrentState->GetTransition(event->GetNameOfClass(), MapToEventVariant(event));

  if (transition.IsNotNull())
  {
    // iterate over all actions in this transition and execute them
    ActionVectorType actions = transition->GetActions();
    bool success = false;
    for (ActionVectorType::iterator it = actions.begin(); it != actions.end(); ++it)
    {

      success |= ExecuteAction(*it, event); // treat an event as handled if at least one of the actions is executed successfully
    }
    if (success || actions.empty())  // an empty action list is always successful
    {
      // perform state change
      m_CurrentState = transition->GetNextState();
      //MITK_INFO<< "StateChange: " << m_CurrentState->GetName();
    }
    return success;
  }
  else
  {
    return false; // no transition found that matches event
  }
}

void mitk::EventStateMachine::ConnectActionsAndFunctions()
{
  MITK_WARN<< "ConnectActionsAndFunctions in DataInteractor not implemented.\n DataInteractor will not be able to process any events.";
}

bool mitk::EventStateMachine::ExecuteAction(StateMachineAction* action, InteractionEvent* event)
{

  if (action == NULL)
  {
    return false;
  }

  bool retVal = false;
  // Maps Action-Name to Functor and executes the Functor.
  ActionDelegatesMapType::iterator delegateIter = m_ActionDelegatesMap.find(action->GetActionName());
  if (delegateIter != m_ActionDelegatesMap.end())
  {
    retVal = delegateIter->second->Execute(action, event);
  }
  else
  {
    // try the legacy system
    std::map<std::string, TActionFunctor*>::iterator functionIter = m_ActionFunctionsMap.find(action->GetActionName());
    if (functionIter != m_ActionFunctionsMap.end())
    {
      retVal = functionIter->second->DoAction(action, event);
    }
  }
  return retVal;
}

mitk::StateMachineState* mitk::EventStateMachine::GetCurrentState() const
{
  return m_CurrentState.GetPointer();
}

bool mitk::EventStateMachine::FilterEvents(InteractionEvent* interactionEvent, DataNode* dataNode)
{
  if (dataNode == NULL)
  {
    MITK_WARN<< "EventStateMachine: Empty DataNode received along with this Event " << interactionEvent;
    return false;
  }
  bool visible = false;
  if (dataNode->GetPropertyList()->GetBoolProperty("visible", visible) == false)
  { //property doesn't exist
    return false;
  }
  return visible;
}
