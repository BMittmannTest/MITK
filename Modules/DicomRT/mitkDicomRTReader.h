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


#ifndef mitkDicomRTReader_h
#define mitkDicomRTReader_h

#include <itkObject.h>
#include <itkObjectFactory.h>
#include <mitkCommon.h>
#include <mitkDicomRTExports.h>
#include "mitkContourModel.h"
#include "mitkContourElement.h"

#include "dcmtk/config/osconfig.h"
#include "dcmtk/ofstd/ofconapp.h"

#include "dcmtk/ofstd/ofcond.h"

#include "dcmtk/dcmrt/drtdose.h"
#include "dcmtk/dcmrt/drtimage.h"
#include "dcmtk/dcmrt/drtplan.h"
#include "dcmtk/dcmrt/drttreat.h"
#include "dcmtk/dcmrt/drtionpl.h"
#include "dcmtk/dcmrt/drtiontr.h"

#include "dcmtk/dcmrt/drtstrct.h"

#include <vtkObjectFactory.h>
#include <vtkPolyData.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkRibbonFilter.h>
#include <vtkPolyDataNormals.h>
#include <vtkSmartPointer.h>
#include <vtkCleanPolyData.h>

class vtkPolyData;
class DcmDataset;
class OFString;
class DRTContourSequence;
class DRTStructureSetIOD;

typedef std::vector<mitk::ContourModel::Pointer> ContourModelVector;

namespace mitk
{

  class mitkDicomRT_EXPORT DicomRTReader: public itk::Object
  {
    class BeamEntry
    {
    public:
      BeamEntry()
      {
        Number=-1;
        SrcAxisDistance=0.0;
        GantryAngle=0.0;
        PatientSupportAngle=0.0;
        BeamLimitingDeviceAngle=0.0;
        LeafJawPositions[0][0]=0.0;
        LeafJawPositions[0][1]=0.0;
        LeafJawPositions[1][0]=0.0;
        LeafJawPositions[1][1]=0.0;
      }
      unsigned int Number;
      std::string Name;
      std::string Type;
      std::string Description;
      double SrcAxisDistance;
      double GantryAngle;
      double PatientSupportAngle;
      double BeamLimitingDeviceAngle;
      double LeafJawPositions[2][2];
    };

    class RoiEntry
    {
    public:
      RoiEntry();
      virtual ~RoiEntry();
      RoiEntry(const RoiEntry& src);
      RoiEntry &operator=(const RoiEntry &src);

      void SetPolyData(vtkPolyData* roiPolyData);

      unsigned int Number;
      std::string  Name;
      std::string  Description;
      double       DisplayColor[3];
      vtkPolyData* PolyData;
    };

  public:

    mitkClassMacro( DicomRTReader, itk::Object );

    itkNewMacro( Self );

    ContourModelVector ReadDicomFile(char* filename);
    ContourModelVector ReadStructureSet(DcmDataset* dataset);
    int LoadRTPlan(DcmDataset* dataset);
    int LoadRTDose(DcmDataset* dataset);
    size_t GetNumberOfRois();
    RoiEntry* FindRoiByNumber(int roiNumber);
    OFString GetReferencedFrameOfReferenceSOPInstanceUID(DRTStructureSetIOD &structSetObject);

    bool Equals(mitk::ContourModel::Pointer first, mitk::ContourModel::Pointer second);

    /**
    * Virtual destructor.
    */
    virtual ~DicomRTReader();

  protected:

    std::vector<RoiEntry> RoiSequenceVector;
    std::vector<BeamEntry> BeamSequenceVector;
    ContourModelVector contourVector;

    /**
    * Constructor.
    */
    DicomRTReader();

  };



}

#endif