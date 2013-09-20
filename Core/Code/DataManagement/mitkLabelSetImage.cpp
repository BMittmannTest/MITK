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

#include "mitkLabelSetImage.h"

#include "mitkImageAccessByItk.h"
#include "mitkInteractionConst.h"
#include "mitkRenderingManager.h"
#include "mitkColormapProperty.h"
#include "mitkImageCast.h"
#include "mitkImageReadAccessor.h"
#include "mitkLookupTableProperty.h"
#include "mitkPadImageFilter.h"

#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkCell.h>

#include <itkImageRegionIterator.h>
#include <itkQuadEdgeMesh.h>
#include <itkTriangleMeshToBinaryImageFilter.h>

mitk::LabelSetImage::LabelSetImage() : mitk::Image()
{
  this->CreateDefaultLabelSet();
}

mitk::LabelSetImage::LabelSetImage(Image::Pointer image) : mitk::Image()
{
  try
  {
    this->Initialize(image);
    mitk::ImageReadAccessor imgA(image, image->GetVolumeData(0));
    this->SetVolume(imgA.GetData());
    this->CreateDefaultLabelSet();
  }
  catch(mitk::Exception e)
  {
    mitkReThrow(e) << "Cannot access image data while constructing labelset image";
  }
}

mitk::LabelSetImage::~LabelSetImage()
{
}

void mitk::LabelSetImage::CreateDefaultLabelSet()
{
  m_LabelSet = mitk::LabelSet::New();
  mitk::Color color;
  color.SetRed(0.0);
  color.SetGreen(0.0);
  color.SetBlue(0.0);
  mitk::Label::Pointer label = mitk::Label::New();
  label->SetColor(color);
  label->SetName("Exterior");
  label->SetExterior(true);
  label->SetOpacity(0.0);
  label->SetLocked(false);
  m_LabelSet->AddLabel(*label);

  mitk::LookupTable::Pointer lut = mitk::LookupTable::New();
  lut->SetActiveColormap(mitk::ColormapProperty::CM_MULTILABEL);

  mitk::LookupTableProperty::Pointer lutProp = mitk::LookupTableProperty::New();
  lutProp->SetLookupTable(lut);

  this->GetPropertyList()->SetProperty( "LookupTable", lutProp );
}

void mitk::LabelSetImage::SetName(const std::string& name)
{
  this->m_LabelSet->SetName(name);
}

std::string mitk::LabelSetImage::GetLabelSetName()
{
  return this->m_LabelSet->GetName();
}

void mitk::LabelSetImage::SetLabelSetLastModified(const std::string& name)
{
  this->m_LabelSet->SetLastModified(name);
}

std::string mitk::LabelSetImage::GetLabelSetLastModified()
{
  return this->m_LabelSet->GetLastModified();
}

void mitk::LabelSetImage::CalculateLabelVolume(int index)
{
// todo
}

void mitk::LabelSetImage::EraseLabel(int index, bool reorder)
{
  AccessByItk_2(this, EraseLabelProcessing, index, reorder);
  this->Modified();
}

bool mitk::LabelSetImage::Concatenate(mitk::LabelSetImage* other)
{
  const unsigned int* otherDims = other->GetDimensions();
  const unsigned int* thisDims = this->GetDimensions();
  if ( (otherDims[0] != thisDims[0]) || (otherDims[1] != thisDims[1]) || (otherDims[2] != thisDims[2]) )
      return false;

  AccessByItk_1(this, ConcatenateProcessing, other);

  mitk::LabelSet::ConstPointer ls = other->GetConstLabelSet();
  for( int i=1; i<ls->GetNumberOfLabels(); i++) // skip exterior
  {
      this->AddLabel( *ls->GetLabel(i) );
  }

  this->Modified();

  return true;
}

void mitk::LabelSetImage::ClearBuffer()
{
  AccessByItk(this, ClearBufferProcessing);
  this->Modified();
}

template < typename LabelSetImageType >
void mitk::LabelSetImage::CalculateCenterOfMassProcessing(LabelSetImageType* itkImage, int index)
{
  // for now, we just retrieve the voxel in the middle
  typename typedef itk::ImageRegionConstIterator< LabelSetImageType > IteratorType;
  IteratorType iter( itkImage, itkImage->GetLargestPossibleRegion() );
  iter.GoToBegin();

  typename std::vector< LabelSetImageType::IndexType > indexVector;

  while ( !iter.IsAtEnd() )
  {
    if ( iter.Get() == static_cast<int>(index) )
    {
      indexVector.push_back(iter.GetIndex());
    }
    ++iter;
  }

  mitk::Point3D pos;
  pos.Fill(0.0);

  if (!indexVector.empty())
  {
   typename itk::ImageRegionConstIteratorWithIndex< LabelSetImageType >::IndexType centerIndex;
   centerIndex = indexVector.at(indexVector.size()/2);
   pos[0] = centerIndex[0];
   pos[1] = centerIndex[1];
   pos[2] = centerIndex[2];
  }

  this->m_LabelSet->SetLabelCenterOfMassIndex(index, pos);
  this->GetSlicedGeometry()->IndexToWorld(pos, pos);
  this->m_LabelSet->SetLabelCenterOfMassCoordinates(index, pos);
}

template < typename LabelSetImageType >
void mitk::LabelSetImage::ClearBufferProcessing(LabelSetImageType* itkImage)
{
  itkImage->FillBuffer(0);
}

template < typename LabelSetImageType >
void mitk::LabelSetImage::ConcatenateProcessing(LabelSetImageType* itkTarget, mitk::LabelSetImage* other)
{
  typename LabelSetImageType::Pointer itkSource = LabelSetImageType::New();
  mitk::CastToItkImage( other, itkSource );

  typedef itk::ImageRegionConstIterator< LabelSetImageType > ConstIteratorType;
  typedef itk::ImageRegionIterator< LabelSetImageType >      IteratorType;

  ConstIteratorType sourceIter( itkSource, itkSource->GetLargestPossibleRegion() );
  IteratorType      targetIter( itkTarget, itkTarget->GetLargestPossibleRegion() );

  int numberOfTargetLabels = this->GetNumberOfLabels() - 1; // skip exterior
  sourceIter.GoToBegin();
  targetIter.GoToBegin();

  while (!sourceIter.IsAtEnd())
  {
    int sourceValue = static_cast<int>(sourceIter.Get());
    int targetValue =  static_cast<int>(targetIter.Get());
    if ( (sourceValue != 0) && !this->GetLabelLocked(targetValue) ) // skip exterior and locked labels
    {
        targetIter.Set( sourceValue + numberOfTargetLabels );
    }
    ++sourceIter;
    ++targetIter;
  }
}

template < typename LabelSetImageType >
void mitk::LabelSetImage::EraseLabelProcessing(LabelSetImageType* itkImage, int index, bool reorder)
{
 typedef itk::ImageRegionIterator< LabelSetImageType > IteratorType;

 IteratorType iter( itkImage, itkImage->GetLargestPossibleRegion() );
 iter.GoToBegin();

 int nLabels = this->GetNumberOfLabels();

 while (!iter.IsAtEnd())
 {
  int value = static_cast<int>(iter.Get());

  if (value == index)
  {
    iter.Set( 0 );
  }
  else if (reorder)
  {
    for (int i=index+1; i<nLabels; ++i) // for each index up-to the end
    {
      if (value == i)
      {
        iter.Set( i-1 );
      }
    }
  }
  ++iter;
 }
}

template < typename LabelSetImageType >
void mitk::LabelSetImage::MergeLabelsProcessing(LabelSetImageType* itkImage, int index)
{
 typedef itk::ImageRegionIterator< LabelSetImageType > IteratorType;

 const mitk::Label* activeLabel = this->m_LabelSet->GetActiveLabel();

 IteratorType iter( itkImage, itkImage->GetLargestPossibleRegion() );
 iter.GoToBegin();

 while (!iter.IsAtEnd())
 {
   if (static_cast<int>(iter.Get()) == index)
   {
     iter.Set( activeLabel->GetIndex() );
   }
   ++iter;
 }
}

void mitk::LabelSetImage::AddLabel(const mitk::Label& label)
{
  this->m_LabelSet->AddLabel(label);
  mitk::LookupTableProperty::Pointer lutProp = dynamic_cast<mitk::LookupTableProperty*>(this->GetPropertyList()->GetProperty( "LookupTable" ));
  double rgba[4];
  const mitk::Color& color = label.GetColor();
  rgba[0] = color.GetRed();
  rgba[1] = color.GetGreen();
  rgba[2] = color.GetBlue();
  rgba[3] = label.GetOpacity();
  // the active label is the one we just added
  lutProp->GetLookupTable()->SetTableValue( m_LabelSet->GetActiveLabelIndex(), rgba );
  AddLabelEvent.Send();
}

void mitk::LabelSetImage::AddLabel(const std::string& name, const mitk::Color& color)
{
  this->m_LabelSet->AddLabel(name, color);
  mitk::LookupTableProperty::Pointer lutProp = dynamic_cast<mitk::LookupTableProperty*>(this->GetPropertyList()->GetProperty( "LookupTable" ));
  double rgba[4];
  rgba[0] = color.GetRed();
  rgba[1] = color.GetGreen();
  rgba[2] = color.GetBlue();
  rgba[3] = m_LabelSet->GetActiveLabel()->GetOpacity();
  lutProp->GetLookupTable()->SetTableValue( m_LabelSet->GetActiveLabelIndex(), rgba );
  AddLabelEvent.Send();
}

const mitk::Color& mitk::LabelSetImage::GetLabelColor(int index)
{
   return this->m_LabelSet->GetLabelColor(index);
}

void mitk::LabelSetImage::SetLabelColor(int index, const mitk::Color& color)
{
  mitk::LookupTableProperty::Pointer lutProp = dynamic_cast<mitk::LookupTableProperty*>(this->GetPropertyList()->GetProperty( "LookupTable" ));
  double rgba[4];
  lutProp->GetLookupTable()->GetTableValue( index, rgba );
  rgba[0] = color.GetRed();
  rgba[1] = color.GetGreen();
  rgba[2] = color.GetBlue();
  rgba[3] = this->m_LabelSet->GetLabelOpacity(index);
  lutProp->GetLookupTable()->SetTableValue( index, rgba );
  this->m_LabelSet->SetLabelColor(index, color);
  ModifyLabelEvent.Send(index);
}

float mitk::LabelSetImage::GetLabelOpacity(int index)
{
  return this->m_LabelSet->GetLabelOpacity(index);
}

float mitk::LabelSetImage::GetLabelVolume(int index)
{
  return this->m_LabelSet->GetLabelVolume(index);
}

void mitk::LabelSetImage::SetLabelOpacity(int index, float value)
{
  mitk::LookupTableProperty::Pointer lutProp = dynamic_cast<mitk::LookupTableProperty*>(this->GetPropertyList()->GetProperty( "LookupTable" ));
  double rgba[4];
  lutProp->GetLookupTable()->GetTableValue(index,rgba);
  rgba[3] = value;
  lutProp->GetLookupTable()->SetTableValue(index,rgba);
  m_LabelSet->SetLabelOpacity(index,value);
  ModifyLabelEvent.Send(index);
}

void mitk::LabelSetImage::SetLabelVolume(int index, float value)
{
  m_LabelSet->SetLabelVolume(index,value);
  ModifyLabelEvent.Send(index);
}

void mitk::LabelSetImage::RenameLabel(int index, const std::string& name, const mitk::Color& color)
{
  m_LabelSet->SetLabelName(index, name);
  m_LabelSet->SetLabelColor(index, color);
  mitk::LookupTableProperty::Pointer lutProp = dynamic_cast<mitk::LookupTableProperty*>(this->GetPropertyList()->GetProperty( "LookupTable" ));
  double rgba[4];
  lutProp->GetLookupTable()->GetTableValue(index,rgba);
  rgba[0] = color.GetRed();
  rgba[1] = color.GetGreen();
  rgba[2] = color.GetBlue();
  lutProp->GetLookupTable()->SetTableValue(index,rgba);
  ModifyLabelEvent.Send(index);
}

void mitk::LabelSetImage::SetLabelName(int index, const std::string& name)
{
  this->m_LabelSet->SetLabelName(index, name);
  ModifyLabelEvent.Send(index);
}

void mitk::LabelSetImage::SetAllLabelsLocked(bool value)
{
  this->m_LabelSet->SetAllLabelsLocked(value);
  AllLabelsModifiedEvent.Send();
}

void mitk::LabelSetImage::SetAllLabelsVisible(bool value)
{
  this->m_LabelSet->SetAllLabelsVisible(value);
  mitk::LookupTableProperty::Pointer lutProp = dynamic_cast<mitk::LookupTableProperty*>(this->GetPropertyList()->GetProperty( "LookupTable" ));
  mitk::LabelSet::LabelContainerConstIteratorType _it = this->m_LabelSet->IteratorBegin();
  mitk::LabelSet::LabelContainerConstIteratorType _end = this->m_LabelSet->IteratorEnd();
  double rgba[4];
  for(; _it!=_end; ++_it)
  {
      lutProp->GetLookupTable()->GetTableValue((*_it)->GetIndex(),rgba);
      rgba[3] = value ? (*_it)->GetOpacity() : 0.0;
      lutProp->GetLookupTable()->SetTableValue((*_it)->GetIndex(),rgba);
  }
  AllLabelsModifiedEvent.Send();
}

void mitk::LabelSetImage::RemoveLabel(int index)
{
  this->EraseLabel(index, true);
  this->m_LabelSet->RemoveLabel(index);
  this->ResetLabels();
  RemoveLabelEvent.Send();
}

void mitk::LabelSetImage::ResetLabels()
{
  this->m_LabelSet->ResetLabels();
  mitk::LookupTableProperty::Pointer lutProp = dynamic_cast<mitk::LookupTableProperty*>(this->GetPropertyList()->GetProperty( "LookupTable" ));
  mitk::LabelSet::LabelContainerConstIteratorType _it = this->m_LabelSet->IteratorBegin();
  mitk::LabelSet::LabelContainerConstIteratorType _end = this->m_LabelSet->IteratorEnd();
  double rgba[4];
  for(; _it!=_end; ++_it)
  {
    const mitk::Color& color = (*_it)->GetColor();
    rgba[0] = color.GetRed();
    rgba[1] = color.GetGreen();
    rgba[2] = color.GetBlue();
    rgba[3] = (*_it)->GetOpacity();
    lutProp->GetLookupTable()->SetTableValue((*_it)->GetIndex(),rgba);
  }
}

void mitk::LabelSetImage::MergeLabels(std::vector<int>& indexes, int index)
{
  this->SetActiveLabel( index, false );
  for (int idx=0; idx<indexes.size(); idx++)
  {
    AccessByItk_1(this, MergeLabelsProcessing, indexes[idx]);
  }

  RemoveLabelEvent.Send();
}

void mitk::LabelSetImage::RemoveLabels(std::vector<int>& indexes)
{
  std::sort(indexes.begin(), indexes.end());

  for (int i=0; i<indexes.size(); i++)
  {
    this->EraseLabel(indexes[i]-i, true);
    this->m_LabelSet->RemoveLabel(indexes[i]-i);
    this->ResetLabels();
  }

  RemoveLabelEvent.Send();
}

void mitk::LabelSetImage::EraseLabels(std::vector<int>& indexes)
{
  std::sort(indexes.begin(), indexes.end());

  for (int i=0; i<indexes.size(); i++)
  {
    this->EraseLabel(indexes[i]-i, false);
  }
}

std::string mitk::LabelSetImage::GetLabelName(int index)
{
  return this->m_LabelSet->GetLabelName(index);
}

unsigned int mitk::LabelSetImage::GetActiveLabelIndex() const
{
  return this->m_LabelSet->GetActiveLabel()->GetIndex();
}

unsigned int mitk::LabelSetImage::GetLabelComponent(int index) const
{
  return this->m_LabelSet->GetLabelComponent(index);
}

unsigned int mitk::LabelSetImage::GetActiveLabelComponent() const
{
  return this->m_LabelSet->GetActiveLabel()->GetComponent();
}

const mitk::Label* mitk::LabelSetImage::GetActiveLabel() const
{
  return this->m_LabelSet->GetActiveLabel();
}

const mitk::Color& mitk::LabelSetImage::GetActiveLabelColor() const
{
  return this->m_LabelSet->GetActiveLabel()->GetColor();
}

bool mitk::LabelSetImage::GetLabelLocked(int index) const
{
  return this->m_LabelSet->GetLabelLocked(index);
}

bool mitk::LabelSetImage::GetLabelSelected(int index) const
{
  return this->m_LabelSet->GetLabelSelected(index);
}

mitk::LabelSet::ConstPointer mitk::LabelSetImage::GetConstLabelSet() const
{
  return m_LabelSet.GetPointer();
}

mitk::LabelSet* mitk::LabelSetImage::GetLabelSet()
{
  return m_LabelSet.GetPointer();
}

void mitk::LabelSetImage::SetLabelSet(mitk::LabelSet* labelset)
{
  if (m_LabelSet != labelset)
  {
    m_LabelSet = labelset;
    this->ResetLabels();
  }
}

bool mitk::LabelSetImage::IsLabelSelected(mitk::Label::Pointer label)
{
  return label->GetSelected();
}

bool mitk::LabelSetImage::GetLabelVisible(int index) const
{
  return this->m_LabelSet->GetLabelVisible(index);
}

void mitk::LabelSetImage::SetLabelVisible(int index, bool value)
{
  this->m_LabelSet->SetLabelVisible(index, value);
  mitk::LookupTableProperty::Pointer lutProp = dynamic_cast<mitk::LookupTableProperty*>(this->GetPropertyList()->GetProperty( "LookupTable" ));
  double rgba[4];
  lutProp->GetLookupTable()->GetTableValue(index,rgba);
  rgba[3] = value ? this->m_LabelSet->GetLabelOpacity(index) : 0.0;
  lutProp->GetLookupTable()->SetTableValue(index,rgba);
  ModifyLabelEvent.Send(index);
}

const mitk::Point3D& mitk::LabelSetImage::GetLabelCenterOfMassIndex(int index, bool update)
{
  if (update)
    AccessByItk_1( this, CalculateCenterOfMassProcessing, index );
  return this->m_LabelSet->GetLabelCenterOfMassIndex(index);
}

const mitk::Point3D& mitk::LabelSetImage::GetLabelCenterOfMassCoordinates(int index, bool update)
{
  if (update)
    AccessByItk_1( this, CalculateCenterOfMassProcessing, index );
  return this->m_LabelSet->GetLabelCenterOfMassCoordinates(index);
}

void mitk::LabelSetImage::SetActiveLabel(int index, bool sendEvent)
{
  this->m_LabelSet->SetActiveLabel(index);
  if (sendEvent)
      ModifyLabelEvent.Send(index);
  this->Modified();
}

void mitk::LabelSetImage::SetLabelLocked(int index, bool value)
{
  this->m_LabelSet->SetLabelLocked(index, value);
}

void mitk::LabelSetImage::SetLabelSelected(int index, bool value)
{
  this->m_LabelSet->SetLabelSelected(index, value);
}

void mitk::LabelSetImage::RemoveAllLabels()
{
  this->m_LabelSet->RemoveAllLabels();
  this->ClearBuffer();
  this->Modified();
  RemoveLabelEvent.Send();
}

int mitk::LabelSetImage::GetNumberOfLabels()
{
  return this->m_LabelSet->GetNumberOfLabels();
}

void mitk::LabelSetImage::PasteMaskOnActiveLabel(mitk::Image* mask, bool forceOverwrite)
{
  if (!mask) return;

  mitk::PadImageFilter::Pointer padImageFilter = mitk::PadImageFilter::New();
  padImageFilter->SetInput(0, mask);
  padImageFilter->SetInput(1, this);
  padImageFilter->SetPadConstant(0);
  padImageFilter->SetBinaryFilter(false);
  padImageFilter->SetLowerThreshold(0);
  padImageFilter->SetUpperThreshold(1);

  padImageFilter->Update();

  mitk::Image::Pointer paddedMask = padImageFilter->GetOutput();

  if (paddedMask.IsNull()) return;

  AccessByItk_2(this, PasteMaskOnActiveLabelProcessing, paddedMask, forceOverwrite);
}

template < typename LabelSetImageType >
void mitk::LabelSetImage::PasteMaskOnActiveLabelProcessing(LabelSetImageType* itkImage, mitk::Image* mask, bool forceOverwrite)
{
  typedef itk::ImageRegionIterator< LabelSetImageType > IteratorType;

  typename LabelSetImageType::Pointer itkMask;
  mitk::CastToItkImage(mask, itkMask);

  typedef itk::ImageRegionConstIterator< LabelSetImageType > SourceIteratorType;
  typedef itk::ImageRegionIterator< LabelSetImageType > TargetIteratorType;

  SourceIteratorType sourceIter( itkMask, itkMask->GetLargestPossibleRegion() );
  sourceIter.GoToBegin();

  TargetIteratorType targetIter( itkImage, itkImage->GetLargestPossibleRegion() );
  targetIter.GoToBegin();

  int activeLabel = this->GetActiveLabelIndex();

  while ( !sourceIter.IsAtEnd() )
  {
    int sourceValue = static_cast<int>(sourceIter.Get());
    int targetValue =  static_cast<int>(targetIter.Get());

    if ( (sourceValue != 0) && (forceOverwrite || !this->GetLabelLocked(targetValue)) ) // skip exterior and locked labels
    {
      targetIter.Set( activeLabel );
    }
    ++sourceIter;
    ++targetIter;
  }

  this->Modified();
}

void mitk::LabelSetImage::PasteSurfaceOnActiveLabel(mitk::Surface* surface, bool forceOverwrite)
{
  if (!surface) return;

  typedef itk::Image< unsigned char, 3 > LabelSetImageType;

  LabelSetImageType::Pointer itkImage;
  mitk::CastToItkImage(this, itkImage);

  vtkPolyData *polydata = surface->GetVtkPolyData();

  vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
  transform->SetMatrix( surface->GetGeometry()->GetVtkTransform()->GetMatrix() );
  transform->Update();

  vtkSmartPointer<vtkTransformPolyDataFilter> transformer = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
  transformer->SetInput(polydata);
  transformer->SetTransform(transform);
  transformer->Update();

  typedef double Coord;
  typedef itk::QuadEdgeMesh< Coord, 3 > MeshType;

  MeshType::Pointer mesh = MeshType::New();
  mesh->SetCellsAllocationMethod( MeshType::CellsAllocatedDynamicallyCellByCell );
  unsigned int numberOfPoints = polydata->GetNumberOfPoints();
  mesh->GetPoints()->Reserve( numberOfPoints );

  vtkPoints* points = polydata->GetPoints();

  MeshType::PointType point;
  for( int i=0; i < numberOfPoints; i++ )
  {
    double* aux = points->GetPoint(i);
    point[0] = aux[0];
    point[1] = aux[1];
    point[2] = aux[2];
    mesh->SetPoint( i, point );
  }

  // Load the polygons into the itk::Mesh
  typedef MeshType::CellAutoPointer     CellAutoPointerType;
  typedef MeshType::CellType            CellType;
  typedef itk::TriangleCell< CellType > TriangleCellType;
  typedef MeshType::PointIdentifier     PointIdentifierType;
  typedef MeshType::CellIdentifier      CellIdentifierType;

  // Read the number of polygons
  CellIdentifierType numberOfPolygons = 0;
  numberOfPolygons = polydata->GetNumberOfPolys();
  vtkCellArray* polys = polydata->GetPolys();

  PointIdentifierType numberOfCellPoints = 3;
  CellIdentifierType i = 0;

  for (i=0; i<numberOfPolygons; i++)
  {
    vtkIdList *cellIds;
    vtkCell *vcell = polydata->GetCell(i);
    cellIds = vcell->GetPointIds();

    CellAutoPointerType cell;
    TriangleCellType * triangleCell = new TriangleCellType;
    PointIdentifierType k;
    for( k = 0; k < numberOfCellPoints; k++ )
    {
      triangleCell->SetPointId( k, cellIds->GetId(k) );
    }

    cell.TakeOwnership( triangleCell );
    mesh->SetCell( i, cell );
  }

  typedef itk::TriangleMeshToBinaryImageFilter<MeshType, LabelSetImageType> TriangleMeshToBinaryImageFilterType;

  TriangleMeshToBinaryImageFilterType::Pointer filter = TriangleMeshToBinaryImageFilterType::New();
  filter->SetInput(mesh);
  filter->SetInfoImage(itkImage);
  filter->SetInsideValue(1);
  filter->SetOutsideValue(0);

  try
  {
    filter->Update();
  }
  catch(mitk::Exception excep)
  {
    mitkReThrow(excep) << "Could not stamp the provided surface";
  }

  LabelSetImageType::Pointer resultImage = filter->GetOutput();
  resultImage->DisconnectPipeline();

// for now, we just retrieve the voxel in the middle
  typedef itk::ImageRegionConstIterator< LabelSetImageType > SourceIteratorType;
  typedef itk::ImageRegionIterator< LabelSetImageType > TargetIteratorType;

  SourceIteratorType sourceIter( resultImage, resultImage->GetLargestPossibleRegion() );
  sourceIter.GoToBegin();

  TargetIteratorType targetIter( itkImage, itkImage->GetLargestPossibleRegion() );
  targetIter.GoToBegin();

  int activeLabel = this->GetActiveLabelIndex();

  while ( !sourceIter.IsAtEnd() )
  {
    int sourceValue = static_cast<int>(sourceIter.Get());
    int targetValue =  static_cast<int>(targetIter.Get());

    if ( (sourceValue != 0) && (forceOverwrite || !this->GetLabelLocked(targetValue)) ) // skip exterior and locked labels
    {
      targetIter.Set( activeLabel );
    }
    ++sourceIter;
    ++targetIter;
  }

  this->Modified();
}
