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

#include "mitkMousePressEvent.h"
#include "mitkMouseReleaseEvent.h"
#include "mitkMouseMoveEvent.h"
#include "mitkVtkPropRenderer.h"
#include "mitkInteractionEventConst.h"
#include "mitkTestingMacros.h"

int mitkInteractionEventTest(int /*argc*/, char* /*argv*/[])
{
  /*
   * Create different Events, fill them with data.
   * And check if isEqual method is implemented properly.
   */
  MITK_TEST_BEGIN("InteractionEvent")

  mitk::VtkPropRenderer::Pointer renderer = NULL;

  mitk::MouseButtons buttonStates =  mitk::LeftMouseButton | mitk::RightMouseButton;
  mitk::MouseButtons eventButton = mitk::LeftMouseButton;
  mitk::ModifierKeys modifiers = mitk::ControlKey | mitk::AltKey;

  mitk::Point2D point;
  point[0] = 17;
  point[1] = 170;

  // MousePress Events
  mitk::MousePressEvent::Pointer me1 = mitk::MousePressEvent::New(renderer,point, buttonStates, modifiers, eventButton);
  mitk::MousePressEvent::Pointer me2 = mitk::MousePressEvent::New(renderer,point, buttonStates, modifiers, eventButton);
  point[0] = 178;
  point[1] = 170;
  mitk::MousePressEvent::Pointer me3 = mitk::MousePressEvent::New(renderer,point, buttonStates, modifiers, eventButton);
  modifiers = mitk::ControlKey;
  mitk::MousePressEvent::Pointer me4 = mitk::MousePressEvent::New(renderer,point, buttonStates, modifiers, eventButton);


  MITK_TEST_CONDITION_REQUIRED(
      me1->MatchesTemplate(me2.GetPointer()) &&
      me1->MatchesTemplate(me3.GetPointer()) &&
      (me2->MatchesTemplate(me3.GetPointer())) &&
      !(me3->MatchesTemplate(me4.GetPointer()))
      , "Checking isEqual and Constructors of mitk::InteractionEvent, mitk::MousePressEvent");

  // MouseReleaseEvents
  mitk::MouseReleaseEvent::Pointer mr1 = mitk::MouseReleaseEvent::New(renderer,point, buttonStates, modifiers, eventButton);
  mitk::MouseReleaseEvent::Pointer mr2 = mitk::MouseReleaseEvent::New(renderer,point, buttonStates, modifiers, eventButton);
  point[0] = 178;
  point[1] = 170;
  mitk::MouseReleaseEvent::Pointer mr3 = mitk::MouseReleaseEvent::New(renderer,point, buttonStates, modifiers, eventButton);
  eventButton = mitk::RightMouseButton;
  mitk::MouseReleaseEvent::Pointer mr4 = mitk::MouseReleaseEvent::New(renderer,point, buttonStates, modifiers, eventButton);


  MITK_TEST_CONDITION_REQUIRED(
      mr1->MatchesTemplate(mr2.GetPointer()) &&
      mr1->MatchesTemplate(mr3.GetPointer()) &&
      (mr2->MatchesTemplate(mr3.GetPointer())) &&
      !(mr3->MatchesTemplate(mr4.GetPointer()))
      , "Checking isEqual and Constructors of mitk::InteractionEvent, mitk::MouseReleaseEvent");


  // MouseMoveEvents
  mitk::MouseMoveEvent::Pointer mm1 = mitk::MouseMoveEvent::New(renderer,point, buttonStates, modifiers);
    point[0] = 178;
  point[1] = 170;
  mitk::MouseMoveEvent::Pointer mm3 = mitk::MouseMoveEvent::New(renderer,point, buttonStates, modifiers);
  modifiers = mitk::AltKey;
  mitk::MouseMoveEvent::Pointer mm4 = mitk::MouseMoveEvent::New(renderer,point, buttonStates, modifiers);



  MITK_TEST_CONDITION_REQUIRED(
      mm1->MatchesTemplate(mm3.GetPointer()) &&
      !(mm3->MatchesTemplate(mm4.GetPointer()))
      , "Checking isEqual and Constructors of mitk::InteractionEvent, mitk::MouseMoveEvent");


  // always end with this!
  MITK_TEST_END()

}
