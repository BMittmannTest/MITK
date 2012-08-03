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

//#define _USE_MATH_DEFINES
#include <QmitkUSDeviceManagerWidget.h>
#include <list>

//QT headers
#include <QColor>

//mitk headers



//itk headers

//microservices
#include <usModuleRegistry.h>
#include <usModule.h>
#include <usServiceTracker.h>


const std::string QmitkUSDeviceManagerWidget::VIEW_ID = "org.mitk.views.QmitkUSDeviceManagerWidget";

QmitkUSDeviceManagerWidget::QmitkUSDeviceManagerWidget(QWidget* parent, Qt::WindowFlags f): QWidget(parent, f) 
{
  m_Controls = NULL;
  CreateQtPartControl(this);

  // get ModuleContext
  mitk::Module* mitkUS = mitk::ModuleRegistry::GetModule("MitkUS");  
  m_MitkUSContext = mitkUS->GetModuleContext();

  //ServiceTracker<mitk::USDevice>* tracker = new ServiceTracker<mitk::USDevice>(m_MitkUSContext, this);

  // Register this Widget as a listener for Registry changes.
  // If devices are registered, unregistered or changed, notifications will go there
  std::string filter = "(&(";
  filter += mitk::ServiceConstants::OBJECTCLASS();
  filter += "=";
  //filter += us_service_interface_iid<mitk::USDevice>();
  filter += "org.mitk.services.UltrasoundDevice)(IsActive=false))";
   m_MitkUSContext->AddServiceListener(this, &QmitkUSDeviceManagerWidget::OnServiceEvent, filter);
}

QmitkUSDeviceManagerWidget::~QmitkUSDeviceManagerWidget()
{
}

//////////////////// INITIALIZATION /////////////////////

void QmitkUSDeviceManagerWidget::CreateQtPartControl(QWidget *parent)
{
  if (!m_Controls)
  {
    // create GUI widgets
    m_Controls = new Ui::QmitkUSDeviceManagerWidgetControls;
    m_Controls->setupUi(parent);
    this->CreateConnections();
  }
}

void QmitkUSDeviceManagerWidget::CreateConnections()
{
  if ( m_Controls )
  {
    connect( m_Controls->m_BtnActivate,   SIGNAL(clicked()), this, SLOT(OnClickedActivateDevice()) );
    connect( m_Controls->m_BtnDisconnect, SIGNAL(clicked()), this, SLOT(OnClickedDisconnectDevice()) );
    connect( m_Controls->m_ConnectedDevices, SIGNAL(currentItemChanged( QListWidgetItem *, QListWidgetItem *)), this, SLOT(OnDeviceSelectionChanged()) );
  }
}




///////////// Methods & Slots Handling Direct Interaction /////////////////

void QmitkUSDeviceManagerWidget::OnClickedActivateDevice()
{
  MITK_INFO << "Activated Device";
  mitk::USDevice::Pointer device = this->GetDeviceForListItem(this->m_Controls->m_ConnectedDevices->currentItem());
  if (device.IsNull()) return;
  if (device->GetIsActive()) device->Deactivate();
  else device->Activate();
}

void QmitkUSDeviceManagerWidget::OnClickedDisconnectDevice(){
  MITK_INFO << "Disconnected Device";
  mitk::USDevice::Pointer device = this->GetDeviceForListItem(this->m_Controls->m_ConnectedDevices->currentItem());
  if (device.IsNull()) return;
  device->Disconnect();
}

void QmitkUSDeviceManagerWidget::OnDeviceSelectionChanged(){
  mitk::USDevice::Pointer device = this->GetDeviceForListItem(this->m_Controls->m_ConnectedDevices->currentItem());
  if (device.IsNull()) return;
  if (device->GetIsActive()) m_Controls->m_BtnActivate->setText("Deactivate");
    else m_Controls->m_BtnActivate->setText("Activate");
}


///////////////// Methods & Slots Handling Logic //////////////////////////

void QmitkUSDeviceManagerWidget::OnServiceEvent(const mitk::ServiceEvent event){
  // Empty ListWidget
  this->m_ListContent.clear();
  m_Controls->m_ConnectedDevices->clear();
  

  
  // get Active Devices
  std::vector<mitk::USDevice::Pointer> devices = this->GetAllRegisteredDevices();
  // Transfer them to the List
  for(std::vector<mitk::USDevice::Pointer>::iterator it = devices.begin(); it != devices.end(); ++it)
  {
    QListWidgetItem *newItem = ConstructItemFromDevice(it->GetPointer());
    //Add new item to QListWidget
    m_Controls->m_ConnectedDevices->addItem(newItem);
    // Construct Link and add to internal List for reference
    QmitkUSDeviceManagerWidget::DeviceListLink link;
    link.device = it->GetPointer();
    link.item = newItem;
    m_ListContent.push_back(link);
  }
}


/////////////////////// HOUSEHOLDING CODE /////////////////////////////////

QListWidgetItem* QmitkUSDeviceManagerWidget::ConstructItemFromDevice(mitk::USDevice::Pointer device){
  QListWidgetItem *result = new QListWidgetItem;
  std::string text = device->GetDeviceManufacturer() + "|" + device->GetDeviceModel();

  if (device->GetIsActive())
  {
    result->foreground().setColor(Qt::blue);
    text += "|(ON)";
  } else text += "|(OFF)";

  result->setText(text.c_str());

  return result;
}


mitk::USDevice::Pointer QmitkUSDeviceManagerWidget::GetDeviceForListItem(QListWidgetItem* item)
{
  for(std::vector<QmitkUSDeviceManagerWidget::DeviceListLink>::iterator it = m_ListContent.begin(); it != m_ListContent.end(); ++it)
  {
    if (item == it->item) return it->device;
  }
  return 0;
}

std::vector <mitk::USDevice::Pointer> QmitkUSDeviceManagerWidget::GetAllRegisteredDevices(){
  
  //Get Service References
  std::list<mitk::ServiceReference> serviceRefs = m_MitkUSContext->GetServiceReferences<mitk::USDevice>();
  
  // Convert Service References to US Devices
  std::vector<mitk::USDevice::Pointer>* result = new std::vector<mitk::USDevice::Pointer>;
  std::list<mitk::ServiceReference>::const_iterator iterator;
  for (iterator = serviceRefs.begin(); iterator != serviceRefs.end(); ++iterator)
  {
    mitk::USDevice::Pointer device = m_MitkUSContext->GetService<mitk::USDevice>(*iterator);
    if (device) result->push_back(device);
  }

  return *result;
}