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

#include "mitkEventConfig.h"
#include <vtkXMLDataElement.h>
#include <mitkStandardFileLocations.h>
#include <vtkObjectFactory.h>
#include <algorithm>
#include <sstream>
#include "mitkEventFactory.h"
#include "mitkInteractionEvent.h"
#include "mitkInternalEvent.h"
#include "mitkInteractionKeyEvent.h"
#include "mitkInteractionEventConst.h"
// us
#include "mitkModule.h"
#include "mitkModuleResource.h"
#include "mitkModuleResourceStream.h"
#include "mitkModuleRegistry.h"

namespace mitk
{
  vtkStandardNewMacro(EventConfig);
}

mitk::EventConfig::EventConfig() :
    m_Errors(false)
{
  if (m_PropertyList.IsNull())
  {
    m_PropertyList = PropertyList::New();
  }
}

mitk::EventConfig::~EventConfig()
{
}

void mitk::EventConfig::InsertMapping(EventMapping mapping)
{
  for (EventListType::iterator it = m_EventList.begin(); it != m_EventList.end(); ++it)
  {
    if ((*it).interactionEvent->MatchesTemplate(mapping.interactionEvent))
    {
      //MITK_INFO<< "Configuration overwritten:" << (*it).variantName;
      m_EventList.erase(it);
      break;
    }
  }
  m_EventList.push_back(mapping);
}

/**
 * @brief Loads the xml file filename and generates the necessary instances.
 **/
bool mitk::EventConfig::LoadConfig(std::string fileName, std::string moduleName)
{
  mitk::Module* module = mitk::ModuleRegistry::GetModule(moduleName);
  mitk::ModuleResource resource = module->GetResource("Interactions/" + fileName);
  if (!resource.IsValid())
  {
    mitkThrow()<< ("Resource not valid. State machine pattern not found:" + fileName);
  }
  mitk::ModuleResourceStream stream(resource);
  this->SetStream(&stream);
  return this->Parse() && !m_Errors;
}

void mitk::EventConfig::StartElement(const char* elementName, const char **atts)
{
  std::string name(elementName);

  if (name == xmlTagConfigRoot)
  {
    //
  }
  else if (name == xmlTagParam)
  {
    std::string name = ReadXMLStringAttribut(xmlParameterName, atts);
    std::string value = ReadXMLStringAttribut(xmlParameterValue, atts);
    m_PropertyList->SetStringProperty(name.c_str(), value.c_str());
  }
  else if (name == xmlTagEventVariant)
  {
    std::string eventClass = ReadXMLStringAttribut(xmlParameterEventClass, atts);
    std::string eventVariant = ReadXMLStringAttribut(xmlParameterName, atts);
    // New list in which all parameters are stored that are given within the <input/> tag
    m_EventPropertyList = PropertyList::New();
    m_EventPropertyList->SetStringProperty(xmlParameterEventClass.c_str(), eventClass.c_str());
    m_EventPropertyList->SetStringProperty(xmlParameterEventVariant.c_str(), eventVariant.c_str());
    m_CurrEventMapping.variantName = eventVariant;
  }
  else if (name == xmlTagAttribute)
  {
    // Attributes that describe an Input Event, such as which MouseButton triggered the event,or which modifier keys are pressed
    std::string name = ReadXMLStringAttribut(xmlParameterName, atts);
    std::string value = ReadXMLStringAttribut(xmlParameterValue, atts);
    m_EventPropertyList->SetStringProperty(name.c_str(), value.c_str());
  }
}

void mitk::EventConfig::EndElement(const char* elementName)
{
  std::string name(elementName);
  // At end of input section, all necessary infos are collected to created an interaction event.
  if (name == xmlTagEventVariant)
  {
    InteractionEvent::Pointer event = EventFactory::CreateEvent(m_EventPropertyList);
    if (event.IsNotNull())
    {
      m_CurrEventMapping.interactionEvent = event;
      InsertMapping(m_CurrEventMapping);
    }
    else
    {
      MITK_WARN<< "EventConfig: Unknown Event-Type in config. Entry skipped: " << name;
    }
  }
}

std::string mitk::EventConfig::ReadXMLStringAttribut(std::string name, const char** atts)
{
  if (atts)
  {
    const char** attsIter = atts;

    while (*attsIter)
    {
      if (name == *attsIter)
      {
        attsIter++;
        return *attsIter;
      }
      attsIter += 2;
    }
  }
  return std::string();
}

mitk::PropertyList::Pointer mitk::EventConfig::GetAttributes() const
{
  return m_PropertyList;
}

std::string mitk::EventConfig::GetMappedEvent(InteractionEvent::Pointer interactionEvent)
{
  // internal events are excluded from mapping
  if (interactionEvent->GetEventClass() == "InternalEvent")
  {
    InternalEvent* internalEvent = dynamic_cast<InternalEvent*>(interactionEvent.GetPointer());
    return internalEvent->GetSignalName();
  }

  for (EventListType::iterator it = m_EventList.begin(); it != m_EventList.end(); ++it)
  {
    if ((*it).interactionEvent->MatchesTemplate(interactionEvent))
    {
      return (*it).variantName;
    }
  }
  // if this part is reached, no mapping has been found,
  // so here we handle key events and map a key event to the string "Std" + letter/code
  // so "A" will be returned as "StdA"
  if (interactionEvent->GetEventClass() == "KeyEvent")
  {
    InteractionKeyEvent* keyEvent = dynamic_cast<InteractionKeyEvent*>(interactionEvent.GetPointer());
    return ("Std" + keyEvent->GetKey());
  }
  return "";
}

void mitk::EventConfig::ClearConfig()
{
  m_EventList.clear();
}

bool mitk::EventConfig::ReadXMLBooleanAttribut(std::string name, const char** atts)
{
  std::string s = ReadXMLStringAttribut(name, atts);
  std::transform(s.begin(), s.end(), s.begin(), ::toupper);

  return s == "TRUE";
}
