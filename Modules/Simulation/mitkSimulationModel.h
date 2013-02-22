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

#ifndef mitkSimulationModel_h
#define mitkSimulationModel_h

#include <SimulationExports.h>
#include <sofa/component/visualmodel/VisualModelImpl.h>
#include <vtkActor.h>
#include <vtkSmartPointer.h>
#include <vtkTexture.h>

namespace mitk
{
  class Simulation_EXPORT SimulationModel : public sofa::component::visualmodel::VisualModelImpl
  {
  public:
    SOFA_CLASS(SimulationModel, sofa::component::visualmodel::VisualModelImpl);

    std::vector<vtkSmartPointer<vtkActor> > GetActors() const;

    bool loadTexture(const std::string& filename);
    bool loadTextures();

  protected:
    void internalDraw(const sofa::core::visual::VisualParams* vparams, bool transparent);

  private:
    SimulationModel();
    ~SimulationModel();

    SimulationModel(const MyType&);
    MyType& operator=(const MyType&);

    void DrawGroup(int ig, const sofa::core::visual::VisualParams* vparams, bool transparent);

    std::vector<vtkSmartPointer<vtkActor> > m_Actors;
    std::map<unsigned int, vtkSmartPointer<vtkTexture> > m_Textures;

    double m_LastTime;
    bool m_LastShowNormals;
    bool m_LastShowVisualModels;
    bool m_LastShowWireFrame;
  };
}

#endif
