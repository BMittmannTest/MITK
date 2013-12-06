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
#ifndef MITKLookupTable_H_HEADER_INCLUDED_C1EBD53D
#define MITKLookupTable_H_HEADER_INCLUDED_C1EBD53D

#include <MitkExports.h>
#include <mitkCommon.h>
#include "vtkLookupTable.h"
#include "vtkSmartPointer.h"
#include <itkDataObject.h>
#include <itkObjectFactory.h>

class vtkColorTransferFunction;
class vtkPiecewiseFunction;

namespace mitk
{
/**
 * @brief The LookupTable class mitk wrapper for a vtkLookupTable
 * @ingroup DataManagement
 *
 * This class can be used to color images with a LookupTable, such as the
 * vtkLookupTable.
 * @note If you want to use this as a property for an mitk::Image, make sure
 * to use the mitk::LookupTableProperty and set the mitk::RenderingModeProperty
 * to a mode which supports lookup tables (e.g. LOOKUPTABLE_COLOR). Make
 * sure to check the documentation of the mitk::RenderingModeProperty. For a
 * code example how to use the mitk::LookupTable check the
 * mitkImageVtkMapper2DLookupTableTest.cpp in Core\Code\Testing.
 */
class MITK_CORE_EXPORT LookupTable : public itk::DataObject
{
public:
    /**
     *@brief Some convenient typedefs.
     */
    typedef unsigned char RawLookupTableType;

    mitkClassMacro( LookupTable, itk::DataObject );

    itkNewMacro( Self );

    /**
     * @returns the associated vtkLookupTable
     */
    virtual vtkLookupTable* GetVtkLookupTable() const;

    virtual RawLookupTableType * GetRawLookupTable() const;

    virtual void SetVtkLookupTable( vtkLookupTable* lut );

    virtual void ChangeOpacityForAll( float opacity );

    virtual void ChangeOpacity(int index, float opacity );

    virtual void GetColor(int, double rgb[3]);

    virtual void GetTableValue(int, double rgba[4]);

    virtual void SetTableValue(int, double rgba[4]);

    virtual void SetActiveColormap(int index);

    itkSetMacro(Window, float);
    itkSetMacro(Level, float);
    itkSetMacro(Opacity, float);

    /*!
    * \brief equality operator implementation
    */
    virtual bool operator==( const mitk::LookupTable& LookupTable ) const;

    /*!
    * \brief non equality operator implementation
    */
    virtual bool operator!=( const LookupTable& LookupTable ) const;

    /*!
    * \brief implementation necessary because operator made
    * private in itk::Object
    */
    virtual LookupTable& operator=( const LookupTable& LookupTable );

    /**
     * Updates the output information of the current object by calling
     * updateOutputInformation of the data objects source object.
     */
    virtual void UpdateOutputInformation( );

    /**
     * Sets the requested Region to the largest possible region.
     * This method is not implemented, since this is the default
     * behaviour of the itk pipeline and we do not support the
     * requested-region mechanism for lookup-tables
     */
    virtual void SetRequestedRegionToLargestPossibleRegion( );

    /**
     * Checks, if the requested region lies outside of the buffered region by
     * calling verifyRequestedRegion().
     */
    virtual bool RequestedRegionIsOutsideOfTheBufferedRegion( );

    /**
     * Checks if the requested region is completely contained in
     * the buffered region. Since we always want to process the lookup
     * table as a whole, this method always returns true
     */
    virtual bool VerifyRequestedRegion( );

    /**
     * This method has no effect for lookup tables, since we do
     * not support the region-mechanism
     */
    virtual void SetRequestedRegion(const itk::DataObject *data );

    LookupTable();
    virtual ~LookupTable();

    void CreateColorTransferFunction(vtkColorTransferFunction*& colorFunction);
    void CreateOpacityTransferFunction(vtkPiecewiseFunction*& opacityFunction);
    void CreateGradientTransferFunction(vtkPiecewiseFunction*& gradientFunction);

protected:

    void PrintSelf(std::ostream &os, itk::Indent indent) const;

    LookupTable(const LookupTable& other);

    virtual void BuildGrayScaleLookupTable();
    virtual void BuildLegacyBinaryLookupTable();
    virtual void BuildInverseGrayScaleLookupTable();
    virtual void BuildRedBlackAlphaLookupTable();
    virtual void BuildGreenBlackAlphaLookupTable();
    virtual void BuildBlueBlackAlphaLookupTable();
    virtual void BuildVRMusclesBonesLookupTable();
    virtual void BuildVRRedVesselsLookupTable();
    virtual void BuildSternLookupTable();
    virtual void BuildHotGreenLookupTable();
    virtual void BuildVRBonesLookupTable();
    virtual void BuildCardiacLookupTable();
    virtual void BuildGrayRainbowLookupTable();
    virtual void BuildRainbowLookupTable();
    virtual void BuildHotMetalLookupTable();
    virtual void BuildHotIronLookupTable();
    virtual void BuildSpectrumLookupTable();
    virtual void BuildGEColorLookupTable();
    virtual void BuildFlowLookupTable();
    virtual void BuildLONILookupTable();
    virtual void BuildLONI2LookupTable();
    virtual void BuildAsymmetryLookupTable();
    virtual void BuildPValueLookupTable();
    virtual void BuildMultiLabelLookupTable();

    vtkSmartPointer<vtkLookupTable> m_LookupTable;

    float m_Window;

    float m_Level;

    float m_Opacity;

private:

    virtual itk::LightObject::Pointer InternalClone() const;
};
} // namespace mitk

#endif /* LookupTable_H_HEADER_INCLUDED_C1EBD53D */
