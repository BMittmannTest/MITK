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

#ifndef _mitkLabelSetImageToSurfaceThreadedFilter_H_
#define _mitkLabelSetImageToSurfaceThreadedFilter_H_

#include "mitkSegmentationSink.h"
#include "mitkSurface.h"
#include "mitkLabelSetImageToSurfaceFilter.h"
#include "mitkColorSequenceRainbow.h"

namespace mitk
{

class Segmentation_EXPORT LabelSetImageToSurfaceThreadedFilter : public SegmentationSink
{
  public:

    mitkClassMacro( LabelSetImageToSurfaceThreadedFilter, SegmentationSink )
    mitkAlgorithmNewMacro( LabelSetImageToSurfaceThreadedFilter );

  protected:

    LabelSetImageToSurfaceThreadedFilter();  // use smart pointers
    virtual ~LabelSetImageToSurfaceThreadedFilter();

    virtual void Initialize(const NonBlockingAlgorithm* other = NULL);
    virtual bool ReadyToRun();

    virtual bool ThreadedUpdateFunction(); // will be called from a thread after calling StartAlgorithm

    virtual void ThreadedUpdateSuccessful(); // will be called from a thread after calling StartAlgorithm

  private:

     mitk::ColorSequenceRainbow::Pointer m_ColorSequenceRainbow;
     int m_RequestedLabel;
     Surface::Pointer m_ResultSurface;
   //  std::stringstream m_ResultNodeName;
};

} // namespace

#endif // _mitkLabelSetImageToSurfaceThreadedFilter_H_
