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

#ifndef mitkContourTool_h_Included
#define mitkContourTool_h_Included

#include "mitkCommon.h"
#include "SegmentationExports.h"
#include "mitkFeedbackContourTool.h"

namespace mitk
{

class Image;
class StateMachineAction;
class InteractionEvent;

/**
  \brief Simple contour filling tool.

  \sa FeedbackContourTool
  \sa ExtractImageFilter
  \sa OverwriteSliceImageFilter

  \ingroup Interaction
  \ingroup ToolManagerEtAl

  Fills a visible contour (from FeedbackContourTool) during mouse dragging. When the mouse button
  is released, ContourTool tries to extract a slice from the working image and fill in
  the (filled) contour as a binary image.

  The painting "color" is defined by m_PaintingPixelValue, which is set in the constructor
  (by sub-classes) or during some event (e.g. in OnInvertLogic - when CTRL is pressed).

  \warning Only to be instantiated by mitk::ToolManager.

  $Author$
*/
class Segmentation_EXPORT ContourTool : public FeedbackContourTool
{
  public:

    mitkClassMacro(ContourTool, FeedbackContourTool);

  protected:

    ContourTool(); // purposely hidden
    virtual ~ContourTool();

    virtual void Activated();
    virtual void Deactivated();

    /**
     Set the label under clicked position as the active one.
    */
   // virtual bool OnChangeActiveLabel(StateMachineAction*, InteractionEvent* interactionEvent );

    virtual bool OnMousePressed( StateMachineAction*, InteractionEvent* interactionEvent );
    virtual bool OnMouseMoved( StateMachineAction*, InteractionEvent* interactionEvent );
    virtual bool OnMouseReleased( StateMachineAction*, InteractionEvent* interactionEvent );
    virtual bool OnChangeActiveLabel (StateMachineAction*, InteractionEvent* interactionEvent);
//    virtual bool OnInvertLogic( StateMachineAction*, InteractionEvent* interactionEvent );

    virtual void ConnectActionsAndFunctions();

    int m_PaintingPixelValue;

    bool m_LogicInverted;
};

} // namespace

#endif


