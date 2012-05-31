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


#ifndef MITKUSDevice_H_HEADER_INCLUDED_
#define MITKUSDevice_H_HEADER_INCLUDED_

#include <vector>
#include "mitkUSProbe.h"
#include "mitkUSImageMetadata.h"
#include "mitkUSImage.h"
#include <MitkUSExports.h>
#include <mitkCommon.h>
#include <mitkImageSource.h>
#include <itkObjectFactory.h>


namespace mitk {

    /**Documentation
    * \brief A device holds information about it's model, make and the connected probes. It is the
    * common super class for all devices and acts as an image source for mitkUSImages. It is the base class
    * for all US Devices, and every new device should extend it.
    * \ingroup US
    */
   class MitkUS_EXPORT USDevice : public mitk::ImageSource
    {
    public:
      mitkClassMacro(USDevice, mitk::ImageSource);
      /**
      * \brief Enforces minimal Metadata to be set. The isVideoOnly flag indicates that this class 
      *        only handles a videostream and does not receive Metadata from the physical device itself.
      */
      mitkNewMacro3Param(Self, std::string, std::string, bool);


      

      /**
      * \brief Add a probe to the device without connecting to it.
      *  This should usually be done before connecting to the probe.
      */
      virtual void AddProbe(mitk::USProbe::Pointer probe);

      /**
      * \brief Connect to a probe and activate it. The probe should be added first.
      *  Usually, a VideoDevice will simply add a probe it wants to connect to,
      *  but an SDK Device might require adding a probe first.
      */
      virtual void ActivateProbe(mitk::USProbe::Pointer probe);

      /**
      * \brief Deactivates the currently active probe.
      */
      virtual void DeactivateProbe();

      /**
      * \brief Removes a probe from the ist of currently added probes.
      */
      //virtual void removeProbe(mitk::USProbe::Pointer probe);
    
      /**
      *  \brief Returns a vector containing all connected probes.
      */
      std::vector<mitk::USProbe::Pointer> GetConnectedProbes();

      /**
      *\brief return the output (output with id 0) of the filter
      */
      USImage* GetOutput(void);

      /**
      *\brief return the output with id idx of the filter
      */
      USImage* GetOutput(unsigned int idx);

      
      /**
      *\brief Graft the specified DataObject onto this ProcessObject's output.
      *
      * See itk::ImageSource::GraftNthOutput for details
      */
      virtual void GraftNthOutput(unsigned int idx, itk::DataObject *graft);

      /**
      * \brief Graft the specified DataObject onto this ProcessObject's output.
      *
      * See itk::ImageSource::Graft Output for details
      */
      virtual void GraftOutput(itk::DataObject *graft);

      /**
      * \brief Make a DataObject of the correct type to used as the specified output.
      *
      * This method is automatically called when DataObject::DisconnectPipeline()
      * is called.  DataObject::DisconnectPipeline, disconnects a data object
      * from being an output of its current source.  When the data object
      * is disconnected, the ProcessObject needs to construct a replacement
      * output data object so that the ProcessObject is in a valid state.
      * Subclasses of USImageVideoSource that have outputs of different
      * data types must overwrite this method so that proper output objects
      * are created.
      */
      virtual DataObjectPointer MakeOutput(unsigned int idx);

      //########### GETTER & SETTER ##################//

      /**
      * \brief Returns the currently active probe or null, if none is active
      */
      itkGetMacro(ActiveProbe, mitk::USProbe::Pointer);
      std::string GetDeviceManufacturer();
      std::string GetDeviceModel();
      std::string GetDeviceComment();
      bool GetIsVideoOnly();

    protected:
      mitk::USProbe::Pointer m_ActiveProbe;
      std::vector<mitk::USProbe::Pointer> m_ConnectedProbes; 
      /**
      * \brief This metadata set is privately used to imprint USImages with Metadata later.
      *        At instantiation time, it only contains Information about the Device,
      *        At scan time, it integrates this data with the probe information and imprints it on 
      *        the produced images. This field is intentionally hidden from outside interference.
      */
      mitk::USImageMetadata::Pointer m_Metadata;

      
      /**
      * \brief Enforces minimal Metadata to be set. The isVideoOnly flag indicates that this class 
      *        only handles a videostream and does not recieve Metadata from the physical device itself.
      */
      USDevice(std::string manufacturer, std::string model, bool isVideoOnly);
      virtual ~USDevice();

      /**
      *  \brief Grabs the next frame from the Video input. This method is called internally, whenever Update() is invoked by an Output.
      */
       void GenerateData();

    };
} // namespace mitk
#endif
