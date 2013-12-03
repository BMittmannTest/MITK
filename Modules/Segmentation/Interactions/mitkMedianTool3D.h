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

#ifndef mitkMedianTool3D_h_Included
#define mitkMedianTool3D_h_Included

#include "mitkAutoSegmentationTool.h"
#include "mitkLegacyAdaptors.h"
#include "SegmentationExports.h"
#include "mitkDataNode.h"
#include "mitkLabelSetImage.h"
#include "mitkStateEvent.h"

#include "itkImage.h"


namespace us {
class ModuleResource;
}


namespace mitk
{

/**
  \brief Segmentation smoothing tool.

  This tool applies a median smoothing filter on the active label
*/
class Segmentation_EXPORT MedianTool3D : public AutoSegmentationTool
{
  public:

    mitkClassMacro(MedianTool3D, AutoSegmentationTool)
    itkNewMacro(MedianTool3D)

    /* icon stuff */
    virtual const char** GetXPM() const;
    virtual const char* GetName() const;
    us::ModuleResource GetIconResource() const;

    void Run();

    void SetRadius(int);
    int GetRadius();

  protected:

    MedianTool3D();
    virtual ~MedianTool3D();

    int m_Radius;

    template < typename TPixel, unsigned int VDimension >
    void InternalProcessing( itk::Image< TPixel, VDimension>* input );
};

} // namespace

#endif
