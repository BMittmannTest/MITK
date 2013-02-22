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

#ifndef MITKMOUSERELEASEEVENT_H_
#define MITKMOUSERELEASEEVENT_H_

#include "itkObject.h"
#include "itkObjectFactory.h"
#include "mitkCommon.h"
#include "mitkInteractionPositionEvent.h"
#include "mitkInteractionEventConst.h"
#include "mitkBaseRenderer.h"
#include "mitkInteractionEvent.h"

#include <MitkExports.h>

namespace mitk
{

  class MITK_CORE_EXPORT MouseReleaseEvent: public InteractionPositionEvent
  {

  public:
    mitkClassMacro(MouseReleaseEvent,InteractionPositionEvent);
    mitkNewMacro5Param(Self, BaseRenderer*, Point2D ,MouseButtons , ModifierKeys, MouseButtons) ;

    ModifierKeys GetModifiers() const;
    MouseButtons GetButtonStates() const;
    void SetModifiers(ModifierKeys modifiers);
    void SetButtonStates(MouseButtons buttons);
    MouseButtons GetEventButton() const;
    void SetEventButton(MouseButtons buttons);
    virtual bool MatchesTemplate(InteractionEvent::Pointer);

  protected:
    MouseReleaseEvent(BaseRenderer*, Point2D, MouseButtons buttonStates, ModifierKeys modifiers, MouseButtons eventButton);
    virtual ~MouseReleaseEvent();

  private:
    MouseButtons m_EventButton;
    MouseButtons m_ButtonStates;
    ModifierKeys m_Modifiers;
  };
} /* namespace mitk */

#endif /* MITKMOUSERELEASEEVENT_H_ */
