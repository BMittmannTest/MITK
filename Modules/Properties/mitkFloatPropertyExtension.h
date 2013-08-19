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

#ifndef mitkFloatPropertyExtension_h
#define mitkFloatPropertyExtension_h

#include "mitkPropertyExtension.h"
#include <PropertiesExports.h>

namespace mitk
{
  class Properties_EXPORT FloatPropertyExtension : public PropertyExtension
  {
  public:
    FloatPropertyExtension();
    FloatPropertyExtension(float minimum, float maximum, float singleStep = 0.1f, int decimals = 2);
    ~FloatPropertyExtension();

    int GetDecimals() const;
    void SetDecimals(int decimals);

    float GetMaximum() const;
    void SetMaximum(float maximum);

    float GetMinimum() const;
    void SetMinimum(float minimum);

    float GetSingleStep() const;
    void SetSingleStep(float singleStep);

  private:
    struct Impl;
    Impl* m_Impl;
  };
}

#endif
