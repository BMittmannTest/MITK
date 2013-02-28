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

#ifndef MITKMOUSEPRESSEVENT_H_
#define MITKMOUSEPRESSEVENT_H_

#include "itkObject.h"
#include "itkObjectFactory.h"
#include "mitkCommon.h"
#include "mitkInteractionEventConst.h"
#include "mitkInteractionPositionEvent.h"
#include "mitkBaseRenderer.h"
#include "mitkInteractionEvent.h"

#include <MitkExports.h>

namespace mitk
{

  class MITK_CORE_EXPORT MousePressEvent: public InteractionPositionEvent
  {

  public:
    mitkClassMacro(MousePressEvent,InteractionPositionEvent)
    mitkNewMacro5Param(Self, BaseRenderer*, const Point2D , MouseButtons , ModifierKeys, MouseButtons)

    ModifierKeys GetModifiers() const;
    MouseButtons GetButtonStates() const;
    void SetModifiers(ModifierKeys modifiers);
    void SetButtonStates(MouseButtons buttons);
    MouseButtons GetEventButton() const;
    void SetEventButton(MouseButtons buttons);

    virtual bool MatchesTemplate(InteractionEvent::Pointer);
    virtual bool IsSuperClassOf(InteractionEvent::Pointer baseClass);

  protected:
    MousePressEvent(BaseRenderer*, const Point2D = Point2D(), mitk::MouseButtons buttonStates = NoButton, mitk::ModifierKeys modifiers =
        NoKey, mitk::MouseButtons eventButton = NoButton);
    virtual ~MousePressEvent();

  private:
    MouseButtons m_EventButton;
    MouseButtons m_ButtonStates;
    ModifierKeys m_Modifiers;
  };
} /* namespace mitk */

#endif /* MITKMOUSEPRESSEVENT_H_ */
