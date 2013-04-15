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

#include "mitkSimulation.h"
#include "mitkSimulationMapper3D.h"
#include "mitkSimulationModel.h"
#include "mitkSimulationObjectFactory.h"
#include "mitkSimulationTemplate.h"
#include <mitkCoreObjectFactory.h>
#include <sofa/component/init.h>
#include <sofa/core/ObjectFactory.h>
#include <sofa/simulation/common/xml/initXml.h>

static void InitializeSOFA()
{
  sofa::component::init();
  sofa::simulation::xml::initXml();

  int SimulationModelClass = sofa::core::RegisterObject("").add<mitk::SimulationModel>();
  sofa::core::ObjectFactory::AddAlias("OglModel", "SimulationModel", true);
  sofa::core::ObjectFactory::AddAlias("VisualModel", "SimulationModel", true);
}

mitk::SimulationObjectFactory::SimulationObjectFactory()
  : m_SimulationIOFactory(SimulationIOFactory::New()),
    m_SimulationTemplateIOFactory(SimulationTemplateIOFactory::New())
{
  itk::ObjectFactoryBase::RegisterFactory(m_SimulationIOFactory);
  itk::ObjectFactoryBase::RegisterFactory(m_SimulationTemplateIOFactory);

  std::string description = "SOFA Scene Files";
  m_FileExtensionsMap.insert(std::pair<std::string, std::string>("*.scn", description));
  m_FileExtensionsMap.insert(std::pair<std::string, std::string>("*.xml", description));

  description = "SOFA Scene File Templates";
  m_FileExtensionsMap.insert(std::pair<std::string, std::string>("*.scn.template", description));
  m_FileExtensionsMap.insert(std::pair<std::string, std::string>("*.xml.template", description));

  InitializeSOFA();
}

mitk::SimulationObjectFactory::~SimulationObjectFactory()
{
  itk::ObjectFactoryBase::UnRegisterFactory(m_SimulationTemplateIOFactory);
  itk::ObjectFactoryBase::UnRegisterFactory(m_SimulationIOFactory);
}

mitk::Mapper::Pointer mitk::SimulationObjectFactory::CreateMapper(mitk::DataNode* node, MapperSlotId slotId)
{
  mitk::Mapper::Pointer mapper;

  if (dynamic_cast<Simulation*>(node->GetData()) != NULL)
  {
    if (slotId == mitk::BaseRenderer::Standard3D)
      mapper = mitk::SimulationMapper3D::New();

    if (mapper.IsNotNull())
      mapper->SetDataNode(node);
  }

  return mapper;
}

const char* mitk::SimulationObjectFactory::GetDescription() const
{
  return "mitk::SimulationObjectFactory";
}

const char* mitk::SimulationObjectFactory::GetFileExtensions()
{
  std::string fileExtensions;
  this->CreateFileExtensions(m_FileExtensionsMap, fileExtensions);
  return fileExtensions.c_str();
}

mitk::CoreObjectFactoryBase::MultimapType mitk::SimulationObjectFactory::GetFileExtensionsMap()
{
  return m_FileExtensionsMap;
}

const char* mitk::SimulationObjectFactory::GetITKSourceVersion() const
{
  return ITK_SOURCE_VERSION;
}

const char* mitk::SimulationObjectFactory::GetSaveFileExtensions()
{
  std::string saveFileExtensions;
  this->CreateFileExtensions(m_FileExtensionsMap, saveFileExtensions);
  return saveFileExtensions.c_str();
}

mitk::CoreObjectFactoryBase::MultimapType mitk::SimulationObjectFactory::GetSaveFileExtensionsMap()
{
  return m_SaveFileExtensionsMap;
}

void mitk::SimulationObjectFactory::SetDefaultProperties(mitk::DataNode* node)
{
  if (node != NULL)
  {
    if (dynamic_cast<Simulation*>(node->GetData()) != NULL)
    {
      SimulationMapper3D::SetDefaultProperties(node);
    }
    else if (dynamic_cast<SimulationTemplate*>(node->GetData()) != NULL)
    {
      SimulationTemplate* simulationTemplate = static_cast<SimulationTemplate*>(node->GetData());
      simulationTemplate->SetProperties(node);
    }
  }
}

void mitk::RegisterSimulationObjectFactory()
{
  static bool alreadyRegistered = false;

  if (!alreadyRegistered)
  {
    mitk::CoreObjectFactory::GetInstance()->RegisterExtraFactory(mitk::SimulationObjectFactory::New());
    alreadyRegistered = true;
  }
}
