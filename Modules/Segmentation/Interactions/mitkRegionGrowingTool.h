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

#ifndef mitkRegionGrowingTool_h_Included
#define mitkRegionGrowingTool_h_Included

#include "mitkFeedbackContourTool.h"
#include "mitkLegacyAdaptors.h"
#include "SegmentationExports.h"

struct mitkIpPicDescriptor;

namespace us {
class ModuleResource;
}

namespace mitk
{

/**
  \brief A slice based region growing tool.

  \sa FeedbackContourTool

  \ingroup Interaction
  \ingroup ToolManagerEtAl

  When the user presses the mouse button, RegionGrowingTool will use the gray values at that position
  to initialize a region growing algorithm (in the affected 2D slice).

  By moving the mouse up and down while the button is still pressed, the user can change the parameters
  of the region growing algorithm (selecting more or less of an object).
  The current result of region growing will always be shown as a contour to the user.

  After releasing the button, the current result of the region growing algorithm will be written to the
  working image of this tool's ToolManager.

  If the first click is <i>inside</i> a segmentation that was generated by region growing (recently),
  the tool will try to cut off a part of the segmentation. For this reason a skeletonization of the segmentation
  is generated and the optimal cut point is determined.

  \warning Only to be instantiated by mitk::ToolManager.

  $Author$
*/
class Segmentation_EXPORT RegionGrowingTool : public FeedbackContourTool
{
  public:

    mitkClassMacro(RegionGrowingTool, FeedbackContourTool);
    itkNewMacro(RegionGrowingTool);

    virtual const char** GetXPM() const;
    virtual us::ModuleResource GetCursorIconResource() const;
    us::ModuleResource GetIconResource() const;

    virtual const char* GetName() const;

  protected:

    RegionGrowingTool(); // purposely hidden
    virtual ~RegionGrowingTool();

    void ConnectActionsAndFunctions();

    virtual void Activated();
    virtual void Deactivated();

    virtual bool OnMousePressed ( StateMachineAction*, InteractionEvent* interactionEvent );
    virtual bool OnMousePressedInside ( StateMachineAction*, InteractionEvent* interactionEvent, mitkIpPicDescriptor* workingPicSlice, int initialWorkingOffset);
    virtual bool OnMousePressedOutside ( StateMachineAction*, InteractionEvent* interactionEvent );
    virtual bool OnMouseMoved   ( StateMachineAction*, InteractionEvent* interactionEvent );
    virtual bool OnMouseReleased( StateMachineAction*, InteractionEvent* interactionEvent );
//virtual bool OnChangeActiveLabel(Action*, const StateEvent*);
    mitkIpPicDescriptor* PerformRegionGrowingAndUpdateContour(int timestep=0);

    Image::Pointer m_ReferenceSlice;
    Image::Pointer m_WorkingSlice;

    ScalarType m_LowerThreshold;
    ScalarType m_UpperThreshold;
    ScalarType m_InitialLowerThreshold;
    ScalarType m_InitialUpperThreshold;

    Point2I m_LastScreenPosition;
    int m_ScreenYDifference;

  private:

    /**
      Smoothes a binary ipPic image with a 5x5 mask. The image borders (some first and last rows) are treated differently.
    */
    mitkIpPicDescriptor* SmoothIPPicBinaryImage( mitkIpPicDescriptor* image, int &contourOfs, mitkIpPicDescriptor* dest = NULL );

    /**
      Helper method for SmoothIPPicBinaryImage. Smoothes a given part of and image.

      \param sourceImage The original binary image.
      \param dest The smoothed image (will be written without bounds checking).
      \param contourOfs One offset of the contour. Is updated if a pixel is changed (which might change the contour).
      \param maskOffsets Memory offsets that describe the smoothing mask.
      \param maskSize Entries of the mask.
      \param startOffset First pixel that should be smoothed using this mask.
      \param endOffset Last pixel that should be smoothed using this mask.
    */
    void SmoothIPPicBinaryImageHelperForRows( mitkIpPicDescriptor* source, mitkIpPicDescriptor* dest, int &contourOfs, int* maskOffsets, int maskSize, int startOffset, int endOffset );

    mitkIpPicDescriptor* m_OriginalPicSlice;
    int m_SeedPointMemoryOffset;

    ScalarType m_VisibleWindow;
    ScalarType m_DefaultWindow;
    ScalarType m_MouseDistanceScaleFactor;

    int m_PaintingPixelValue;
    int m_LastWorkingSeed;

    bool m_FillFeedbackContour;
};

} // namespace

#endif


