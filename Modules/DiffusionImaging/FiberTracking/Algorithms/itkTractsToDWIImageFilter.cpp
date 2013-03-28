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

#include "itkTractsToDWIImageFilter.h"
#include <boost/progress.hpp>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkPoints.h>
#include <vtkPolyLine.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <mitkGibbsRingingArtifact.h>
#include <itkResampleImageFilter.h>
#include <itkNearestNeighborInterpolateImageFunction.h>
#include <itkBSplineInterpolateImageFunction.h>
#include <itkCastImageFilter.h>
#include <itkImageFileWriter.h>
#include <itkRescaleIntensityImageFilter.h>
#include <itkWindowedSincInterpolateImageFunction.h>
#include <itkResampleDwiImageFilter.h>

namespace itk
{
TractsToDWIImageFilter::TractsToDWIImageFilter()
    : m_CircleDummy(false)
    , m_VolumeAccuracy(10)
    , m_Upsampling(1)
    , m_NumberOfRepetitions(1)
    , m_EnforcePureFiberVoxels(false)
    , m_InterpolationShrink(10)
    , m_FiberRadius(20)
    , m_SignalScale(300)
{
    m_Spacing.Fill(2.5); m_Origin.Fill(0.0);
    m_DirectionMatrix.SetIdentity();
    m_ImageRegion.SetSize(0, 10);
    m_ImageRegion.SetSize(1, 10);
    m_ImageRegion.SetSize(2, 10);
}

TractsToDWIImageFilter::~TractsToDWIImageFilter()
{

}

std::vector< TractsToDWIImageFilter::DoubleDwiType::Pointer > TractsToDWIImageFilter::AddKspaceArtifacts( std::vector< DoubleDwiType::Pointer >& images )
{
    // create slice object
    SliceType::Pointer slice = SliceType::New();
    ImageRegion<2> region;
    region.SetSize(0, m_UpsampledImageRegion.GetSize()[0]);
    region.SetSize(1, m_UpsampledImageRegion.GetSize()[1]);
    slice->SetLargestPossibleRegion( region );
    slice->SetBufferedRegion( region );
    slice->SetRequestedRegion( region );
    slice->Allocate();

    boost::progress_display disp(images.size()*images[0]->GetVectorLength()*images[0]->GetLargestPossibleRegion().GetSize(2));
    std::vector< DoubleDwiType::Pointer > outImages;
    for (int i=0; i<images.size(); i++)
    {
        DoubleDwiType::Pointer image = images.at(i);
        DoubleDwiType::Pointer newImage = DoubleDwiType::New();
        newImage->SetSpacing( m_Spacing );
        newImage->SetOrigin( m_Origin );
        newImage->SetDirection( m_DirectionMatrix );
        newImage->SetLargestPossibleRegion( m_ImageRegion );
        newImage->SetBufferedRegion( m_ImageRegion );
        newImage->SetRequestedRegion( m_ImageRegion );
        newImage->SetVectorLength( image->GetVectorLength() );
        newImage->Allocate();

        DiffusionSignalModel<double>* signalModel;
        if (i<m_FiberModels.size())
            signalModel = m_FiberModels.at(i);
        else
            signalModel = m_NonFiberModels.at(i-m_FiberModels.size());

        for (int g=0; g<image->GetVectorLength(); g++)
            for (int z=0; z<image->GetLargestPossibleRegion().GetSize(2); z++)
            {
                ++disp;

                // extract slice from channel g
                for (int y=0; y<image->GetLargestPossibleRegion().GetSize(1); y++)
                    for (int x=0; x<image->GetLargestPossibleRegion().GetSize(0); x++)
                    {
                        SliceType::IndexType index2D;
                        index2D[0]=x; index2D[1]=y;
                        DoubleDwiType::IndexType index3D;
                        index3D[0]=x; index3D[1]=y; index3D[2]=z;

                        SliceType::PixelType pix2D = image->GetPixel(index3D)[g];
                        slice->SetPixel(index2D, pix2D);
                    }

                // fourier transform slice
                itk::FFTRealToComplexConjugateImageFilter< SliceType::PixelType, 2 >::Pointer fft = itk::FFTRealToComplexConjugateImageFilter< SliceType::PixelType, 2 >::New();
                fft->SetInput(slice);
                fft->Update();
                ComplexSliceType::Pointer fSlice = fft->GetOutput();
                fSlice = RearrangeSlice(fSlice);

                // add artifacts
                for (int a=0; a<m_KspaceArtifacts.size(); a++)
                {
                    m_KspaceArtifacts.at(a)->SetT2(signalModel->GetT2());
                    fSlice = m_KspaceArtifacts.at(a)->AddArtifact(fSlice);
                }

                // save k-space slice of s0 image
                if (g==m_FiberModels.at(0)->GetFirstBaselineIndex())
                    for (int y=0; y<fSlice->GetLargestPossibleRegion().GetSize(1); y++)
                        for (int x=0; x<fSlice->GetLargestPossibleRegion().GetSize(0); x++)
                        {
                            DoubleDwiType::IndexType index3D;
                            index3D[0]=x; index3D[1]=y; index3D[2]=z;
                            SliceType::IndexType index2D;
                            index2D[0]=x; index2D[1]=y;
                            double kpix = sqrt(fSlice->GetPixel(index2D).real()*fSlice->GetPixel(index2D).real()+fSlice->GetPixel(index2D).imag()*fSlice->GetPixel(index2D).imag());
                            m_KspaceImage->SetPixel(index3D, m_KspaceImage->GetPixel(index3D)+kpix);
                        }

                // inverse fourier transform slice
                SliceType::Pointer newSlice;
                itk::FFTComplexConjugateToRealImageFilter< SliceType::PixelType, 2 >::Pointer ifft = itk::FFTComplexConjugateToRealImageFilter< SliceType::PixelType, 2 >::New();
                ifft->SetInput(fSlice);
                ifft->Update();
                newSlice = ifft->GetOutput();

                // put slice back into channel g
                for (int y=0; y<fSlice->GetLargestPossibleRegion().GetSize(1); y++)
                    for (int x=0; x<fSlice->GetLargestPossibleRegion().GetSize(0); x++)
                    {
                        DoubleDwiType::IndexType index3D;
                        index3D[0]=x; index3D[1]=y; index3D[2]=z;
                        DoubleDwiType::PixelType pix3D = newImage->GetPixel(index3D);
                        SliceType::IndexType index2D;
                        index2D[0]=x; index2D[1]=y;

                        pix3D[g] = newSlice->GetPixel(index2D);
                        newImage->SetPixel(index3D, pix3D);
                    }
            }
        outImages.push_back(newImage);
    }
    return outImages;
}

TractsToDWIImageFilter::ComplexSliceType::Pointer TractsToDWIImageFilter::RearrangeSlice(ComplexSliceType::Pointer slice)
{
    ImageRegion<2> region = slice->GetLargestPossibleRegion();
    ComplexSliceType::Pointer rearrangedSlice = ComplexSliceType::New();
    rearrangedSlice->SetLargestPossibleRegion( region );
    rearrangedSlice->SetBufferedRegion( region );
    rearrangedSlice->SetRequestedRegion( region );
    rearrangedSlice->Allocate();

    int xHalf = region.GetSize(0)/2;
    int yHalf = region.GetSize(1)/2;

    for (int y=0; y<region.GetSize(1); y++)
        for (int x=0; x<region.GetSize(0); x++)
        {
            SliceType::IndexType idx;
            idx[0]=x; idx[1]=y;
            vcl_complex< double > pix = slice->GetPixel(idx);

            if( idx[0] <  xHalf )
                idx[0] = idx[0] + xHalf;
            else
                idx[0] = idx[0] - xHalf;

            if( idx[1] <  yHalf )
                idx[1] = idx[1] + yHalf;
            else
                idx[1] = idx[1] - yHalf;

            rearrangedSlice->SetPixel(idx, pix);
        }

    return rearrangedSlice;
}

void TractsToDWIImageFilter::GenerateData()
{
    // check input data
    if (m_FiberBundle.IsNull())
        itkExceptionMacro("Input fiber bundle is NULL!");

    int numFibers = m_FiberBundle->GetNumFibers();
    if (numFibers<=0)
        itkExceptionMacro("Input fiber bundle contains no fibers!");

    if (m_FiberModels.empty())
        itkExceptionMacro("No diffusion model for fiber compartments defined!");

    if (m_NonFiberModels.empty())
        itkExceptionMacro("No diffusion model for non-fiber compartments defined!");

    int baselineIndex = m_FiberModels[0]->GetFirstBaselineIndex();
    if (baselineIndex<0)
        itkExceptionMacro("No baseline index found!");

    // determine k-space undersampling
    for (int i=0; i<m_KspaceArtifacts.size(); i++)
        if ( dynamic_cast<mitk::GibbsRingingArtifact<double>*>(m_KspaceArtifacts.at(i)) )
            m_Upsampling = dynamic_cast<mitk::GibbsRingingArtifact<double>*>(m_KspaceArtifacts.at(i))->GetKspaceCropping();
    if (m_Upsampling<1)
        m_Upsampling = 1;

    if (m_TissueMask.IsNotNull())
    {
        // use input tissue mask
        m_Spacing = m_TissueMask->GetSpacing();
        m_Origin = m_TissueMask->GetOrigin();
        m_DirectionMatrix = m_TissueMask->GetDirection();
        m_ImageRegion = m_TissueMask->GetLargestPossibleRegion();

        if (m_Upsampling>1)
        {
            ImageRegion<3> region = m_ImageRegion;
            region.SetSize(0, m_ImageRegion.GetSize(0)*m_Upsampling);
            region.SetSize(1, m_ImageRegion.GetSize(1)*m_Upsampling);
            mitk::Vector3D spacing = m_Spacing;
            spacing[0] /= m_Upsampling;
            spacing[1] /= m_Upsampling;
            itk::RescaleIntensityImageFilter<ItkUcharImgType,ItkUcharImgType>::Pointer rescaler = itk::RescaleIntensityImageFilter<ItkUcharImgType,ItkUcharImgType>::New();
            rescaler->SetInput(0,m_TissueMask);
            rescaler->SetOutputMaximum(100);
            rescaler->SetOutputMinimum(0);
            rescaler->Update();

            itk::ResampleImageFilter<ItkUcharImgType, ItkUcharImgType>::Pointer resampler = itk::ResampleImageFilter<ItkUcharImgType, ItkUcharImgType>::New();
            resampler->SetInput(rescaler->GetOutput());
            resampler->SetOutputParametersFromImage(m_TissueMask);
            resampler->SetSize(region.GetSize());
            resampler->SetOutputSpacing(spacing);
            resampler->Update();
            m_TissueMask = resampler->GetOutput();
        }
        MITK_INFO << "Using tissue mask";
    }

    // initialize output dwi image
    OutputImageType::Pointer outImage = OutputImageType::New();
    outImage->SetSpacing( m_Spacing );
    outImage->SetOrigin( m_Origin );
    outImage->SetDirection( m_DirectionMatrix );
    outImage->SetLargestPossibleRegion( m_ImageRegion );
    outImage->SetBufferedRegion( m_ImageRegion );
    outImage->SetRequestedRegion( m_ImageRegion );
    outImage->SetVectorLength( m_FiberModels[0]->GetNumGradients() );
    outImage->Allocate();
    OutputImageType::PixelType temp;
    temp.SetSize(m_FiberModels[0]->GetNumGradients());
    temp.Fill(0.0);
    outImage->FillBuffer(temp);

    // is input slize size a power of two?
    int x=2; int y=2;
    while (x<m_ImageRegion.GetSize(0))
        x *= 2;
    while (y<m_ImageRegion.GetSize(1))
        y *= 2;

    // if not, adjust size and dimension (needed for FFT); zero-padding
    if (x!=m_ImageRegion.GetSize(0))
    {
        MITK_INFO << "Adjusting image width: " << m_ImageRegion.GetSize(0) << " --> " << x << " --> " << x*m_Upsampling;
        m_ImageRegion.SetSize(0, x);
    }
    if (y!=m_ImageRegion.GetSize(1))
    {
        MITK_INFO << "Adjusting image height: " << m_ImageRegion.GetSize(1) << " --> " << y << " --> " << y*m_Upsampling;
        m_ImageRegion.SetSize(1, y);
    }

    // initialize k-space image
    m_KspaceImage = ItkDoubleImgType::New();
    m_KspaceImage->SetSpacing( m_Spacing );
    m_KspaceImage->SetOrigin( m_Origin );
    m_KspaceImage->SetDirection( m_DirectionMatrix );
    m_KspaceImage->SetLargestPossibleRegion( m_ImageRegion );
    m_KspaceImage->SetBufferedRegion( m_ImageRegion );
    m_KspaceImage->SetRequestedRegion( m_ImageRegion );
    m_KspaceImage->Allocate();
    m_KspaceImage->FillBuffer(0);

    // apply undersampling to image parameters
    m_UpsampledSpacing = m_Spacing;
    m_UpsampledImageRegion = m_ImageRegion;
    m_UpsampledSpacing[0] /= m_Upsampling;
    m_UpsampledSpacing[1] /= m_Upsampling;
    m_UpsampledImageRegion.SetSize(0, m_ImageRegion.GetSize()[0]*m_Upsampling);
    m_UpsampledImageRegion.SetSize(1, m_ImageRegion.GetSize()[1]*m_Upsampling);

    // everything from here on is using the upsampled image parameters!!!
    if (m_TissueMask.IsNull())
    {
        m_TissueMask = ItkUcharImgType::New();
        m_TissueMask->SetSpacing( m_UpsampledSpacing );
        m_TissueMask->SetOrigin( m_Origin );
        m_TissueMask->SetDirection( m_DirectionMatrix );
        m_TissueMask->SetLargestPossibleRegion( m_UpsampledImageRegion );
        m_TissueMask->SetBufferedRegion( m_UpsampledImageRegion );
        m_TissueMask->SetRequestedRegion( m_UpsampledImageRegion );
        m_TissueMask->Allocate();
        m_TissueMask->FillBuffer(1);
    }

    // resample fiber bundle for sufficient voxel coverage
    double segmentVolume = 0.0001;
    float minSpacing = 1;
    if(m_UpsampledSpacing[0]<m_UpsampledSpacing[1] && m_UpsampledSpacing[0]<m_UpsampledSpacing[2])
        minSpacing = m_UpsampledSpacing[0];
    else if (m_UpsampledSpacing[1] < m_UpsampledSpacing[2])
        minSpacing = m_UpsampledSpacing[1];
    else
        minSpacing = m_UpsampledSpacing[2];
    FiberBundleType fiberBundle = m_FiberBundle->GetDeepCopy();
    fiberBundle->ResampleFibers(minSpacing/m_VolumeAccuracy);
    double mmRadius = m_FiberRadius/1000;
    if (mmRadius>0)
        segmentVolume = M_PI*mmRadius*mmRadius*minSpacing/m_VolumeAccuracy;

    // generate double images to wokr with because we don't want to lose precision
    // we use a separate image for each compartment model
    std::vector< DoubleDwiType::Pointer > compartments;
    for (int i=0; i<m_FiberModels.size()+m_NonFiberModels.size(); i++)
    {
        DoubleDwiType::Pointer doubleDwi = DoubleDwiType::New();
        doubleDwi->SetSpacing( m_UpsampledSpacing );
        doubleDwi->SetOrigin( m_Origin );
        doubleDwi->SetDirection( m_DirectionMatrix );
        doubleDwi->SetLargestPossibleRegion( m_UpsampledImageRegion );
        doubleDwi->SetBufferedRegion( m_UpsampledImageRegion );
        doubleDwi->SetRequestedRegion( m_UpsampledImageRegion );
        doubleDwi->SetVectorLength( m_FiberModels[0]->GetNumGradients() );
        doubleDwi->Allocate();
        DoubleDwiType::PixelType pix;
        pix.SetSize(m_FiberModels[0]->GetNumGradients());
        pix.Fill(0.0);
        doubleDwi->FillBuffer(pix);
        compartments.push_back(doubleDwi);
    }

    double interpFact = 2*atan(-0.5*m_InterpolationShrink);
    double maxVolume = 0;

    vtkSmartPointer<vtkPolyData> fiberPolyData = fiberBundle->GetFiberPolyData();
    vtkSmartPointer<vtkCellArray> vLines = fiberPolyData->GetLines();
    vLines->InitTraversal();

    MITK_INFO << "Generating signal of " << m_FiberModels.size() << " fiber compartments";
    boost::progress_display disp(numFibers);
    for( int i=0; i<numFibers; i++ )
    {
        ++disp;
        vtkIdType   numPoints(0);
        vtkIdType*  points(NULL);
        vLines->GetNextCell ( numPoints, points );
        if (numPoints<2)
            continue;

        for( int j=0; j<numPoints; j++)
        {
            double* temp = fiberPolyData->GetPoint(points[j]);
            itk::Point<float, 3> vertex = GetItkPoint(temp);
            itk::Vector<double> v = GetItkVector(temp);

            itk::Vector<double, 3> dir(3);
            if (j<numPoints-1)
                dir = GetItkVector(fiberPolyData->GetPoint(points[j+1]))-v;
            else
                dir = v-GetItkVector(fiberPolyData->GetPoint(points[j-1]));

            itk::Index<3> idx;
            itk::ContinuousIndex<float, 3> contIndex;
            m_TissueMask->TransformPhysicalPointToIndex(vertex, idx);
            m_TissueMask->TransformPhysicalPointToContinuousIndex(vertex, contIndex);

            double frac_x = contIndex[0] - idx[0];
            double frac_y = contIndex[1] - idx[1];
            double frac_z = contIndex[2] - idx[2];
            if (frac_x<0)
            {
                idx[0] -= 1;
                frac_x += 1;
            }
            if (frac_y<0)
            {
                idx[1] -= 1;
                frac_y += 1;
            }
            if (frac_z<0)
            {
                idx[2] -= 1;
                frac_z += 1;
            }

            frac_x = atan((0.5-frac_x)*m_InterpolationShrink)/interpFact + 0.5;
            frac_y = atan((0.5-frac_y)*m_InterpolationShrink)/interpFact + 0.5;
            frac_z = atan((0.5-frac_z)*m_InterpolationShrink)/interpFact + 0.5;

            // use trilinear interpolation
            itk::Index<3> newIdx;
            for (int x=0; x<2; x++)
            {
                frac_x = 1-frac_x;
                for (int y=0; y<2; y++)
                {
                    frac_y = 1-frac_y;
                    for (int z=0; z<2; z++)
                    {
                        frac_z = 1-frac_z;

                        newIdx[0] = idx[0]+x;
                        newIdx[1] = idx[1]+y;
                        newIdx[2] = idx[2]+z;

                        double frac = frac_x*frac_y*frac_z;

                        // is position valid?
                        if (!m_TissueMask->GetLargestPossibleRegion().IsInside(newIdx) || m_TissueMask->GetPixel(newIdx)<=0)
                            continue;

                        // generate signal for each fiber compartment
                        for (int k=0; k<m_FiberModels.size(); k++)
                        {
                            DoubleDwiType::Pointer doubleDwi = compartments.at(k);
                            m_FiberModels[k]->SetFiberDirection(dir);
                            DoubleDwiType::PixelType pix = doubleDwi->GetPixel(newIdx);
                            pix += segmentVolume*frac*m_FiberModels[k]->SimulateMeasurement();
                            doubleDwi->SetPixel(newIdx, pix );
                            if (pix[baselineIndex]>maxVolume)
                                maxVolume = pix[baselineIndex];
                        }
                    }
                }
            }
        }
    }

    MITK_INFO << "Generating signal of " << m_NonFiberModels.size() << " non-fiber compartments";
    ImageRegionIterator<ItkUcharImgType> it3(m_TissueMask, m_TissueMask->GetLargestPossibleRegion());
    boost::progress_display disp3(m_TissueMask->GetLargestPossibleRegion().GetNumberOfPixels());
    double voxelVolume = m_UpsampledSpacing[0]*m_UpsampledSpacing[1]*m_UpsampledSpacing[2];

    double fact = 1;
    if (m_FiberRadius<0.0001)
        fact = voxelVolume/maxVolume;

    while(!it3.IsAtEnd())
    {
        ++disp3;
        DoubleDwiType::IndexType index = it3.GetIndex();

        if (it3.Get()>0)
        {
            // get fiber volume fraction
            DoubleDwiType::Pointer fiberDwi = compartments.at(0);
            DoubleDwiType::PixelType fiberPix = fiberDwi->GetPixel(index); // intra axonal compartment
            if (fact>1) // auto scale intra-axonal if no fiber radius is specified
            {
                fiberPix *= fact;
                fiberDwi->SetPixel(index, fiberPix);
            }
            double f = fiberPix[baselineIndex];

            if (f>voxelVolume || f>0 && m_EnforcePureFiberVoxels)  // more fiber than space in voxel?
            {
                fiberDwi->SetPixel(index, fiberPix*voxelVolume/f);

                for (int i=1; i<m_FiberModels.size(); i++)
                {
                    DoubleDwiType::PixelType pix; pix.Fill(0.0);
                    compartments.at(i)->SetPixel(index, pix);
                }
            }
            else
            {
                double nonf = voxelVolume-f;    // non-fiber volume
                double inter = 0;
                if (m_FiberModels.size()>1)
                    inter = nonf * f;           // intra-axonal fraction of non fiber compartment scales linearly with f
                double other = nonf - inter;    // rest of compartment

                // adjust non-fiber and intra-axonal signal
                for (int i=1; i<m_FiberModels.size(); i++)
                {
                    DoubleDwiType::Pointer doubleDwi = compartments.at(i);
                    DoubleDwiType::PixelType pix = doubleDwi->GetPixel(index);
                    if (pix[baselineIndex]>0)
                        pix /= pix[baselineIndex];
                    pix *= inter;
                    doubleDwi->SetPixel(index, pix);
                }
                for (int i=0; i<m_NonFiberModels.size(); i++)
                {
                    DoubleDwiType::Pointer doubleDwi = compartments.at(i+m_FiberModels.size());
                    DoubleDwiType::PixelType pix = doubleDwi->GetPixel(index) + m_NonFiberModels[i]->SimulateMeasurement()*other*m_NonFiberModels[i]->GetWeight();
                    doubleDwi->SetPixel(index, pix);
                }
            }
        }
        ++it3;
    }

    // do k-space stuff
    MITK_INFO << "Adjusting complex signal";
    compartments = AddKspaceArtifacts(compartments);

    MITK_INFO << "Summing compartments and adding noise";
    unsigned int window = 0;
    unsigned int min = itk::NumericTraits<unsigned int>::max();
    ImageRegionIterator<DWIImageType> it4 (outImage, outImage->GetLargestPossibleRegion());
    DoubleDwiType::PixelType signal; signal.SetSize(m_FiberModels[0]->GetNumGradients());
    boost::progress_display disp4(outImage->GetLargestPossibleRegion().GetNumberOfPixels());
    while(!it4.IsAtEnd())
    {
        ++disp4;
        DWIImageType::IndexType index = it4.GetIndex();
        signal.Fill(0.0);

        // adjust fiber signal
        for (int i=0; i<m_FiberModels.size(); i++)
            signal += compartments.at(i)->GetPixel(index)*m_SignalScale;

        // adjust non-fiber signal
        for (int i=0; i<m_NonFiberModels.size(); i++)
            signal += compartments.at(m_FiberModels.size()+i)->GetPixel(index)*m_SignalScale;

        DoubleDwiType::PixelType accu = signal; accu.Fill(0.0);
        for (int i=0; i<m_NumberOfRepetitions; i++)
        {
            DoubleDwiType::PixelType temp = signal;
            m_NoiseModel->AddNoise(temp);
            accu += temp;
        }
        signal = accu/m_NumberOfRepetitions;
        for (int i=0; i<signal.Size(); i++)
        {
            if (signal[i]>0)
                signal[i] = floor(signal[i]+0.5);
            else
                signal[i] = ceil(signal[i]-0.5);

            if (!m_FiberModels.at(0)->IsBaselineIndex(i) && signal[i]>window)
                window = signal[i];
            if (!m_FiberModels.at(0)->IsBaselineIndex(i) && signal[i]<min)
                min = signal[i];
        }
        it4.Set(signal);
        ++it4;
    }
    window -= min;
    unsigned int level = window/2 + min;
    m_LevelWindow.SetLevelWindow(level, window);
    this->SetNthOutput(0, outImage);
}

itk::Point<float, 3> TractsToDWIImageFilter::GetItkPoint(double point[3])
{
    itk::Point<float, 3> itkPoint;
    itkPoint[0] = point[0];
    itkPoint[1] = point[1];
    itkPoint[2] = point[2];
    return itkPoint;
}

itk::Vector<double, 3> TractsToDWIImageFilter::GetItkVector(double point[3])
{
    itk::Vector<double, 3> itkVector;
    itkVector[0] = point[0];
    itkVector[1] = point[1];
    itkVector[2] = point[2];
    return itkVector;
}

vnl_vector_fixed<double, 3> TractsToDWIImageFilter::GetVnlVector(double point[3])
{
    vnl_vector_fixed<double, 3> vnlVector;
    vnlVector[0] = point[0];
    vnlVector[1] = point[1];
    vnlVector[2] = point[2];
    return vnlVector;
}


vnl_vector_fixed<double, 3> TractsToDWIImageFilter::GetVnlVector(Vector<float,3>& vector)
{
    vnl_vector_fixed<double, 3> vnlVector;
    vnlVector[0] = vector[0];
    vnlVector[1] = vector[1];
    vnlVector[2] = vector[2];
    return vnlVector;
}

}
