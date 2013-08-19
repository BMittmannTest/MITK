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

#include "mitkTestingMacros.h"
#include <mitkTestingConfig.h>
#include <mitkIOUtil.h>
#include <mitkImageGenerator.h>

int mitkIOUtilTest(int  argc , char* argv[])
{
    // always start with this!
    MITK_TEST_BEGIN("mitkIOUtilTest");
    MITK_TEST_CONDITION_REQUIRED( argc >= 4, "Testing if input parameters are set.");

    std::string pathToImage = argv[1];
    std::string pathToPointSet = argv[2];
    std::string pathToSurface = argv[3];

    mitk::Image::Pointer img1 = mitk::IOUtil::LoadImage(pathToImage);
    MITK_TEST_CONDITION( img1.IsNotNull(), "Testing if image 1 could be loaded");
    mitk::PointSet::Pointer pointset = mitk::IOUtil::LoadPointSet(pathToPointSet);
    MITK_TEST_CONDITION( pointset.IsNotNull(), "Testing if pointset could be loaded");
    mitk::Surface::Pointer surface = mitk::IOUtil::LoadSurface(pathToSurface);
    MITK_TEST_CONDITION( surface.IsNotNull(), "Testing if surface could be loaded");

    std::string outDir = MITK_TEST_OUTPUT_DIR;
    std::string imagePath = outDir+"/diffpic3d.nrrd";
    std::string imagePath2 = outDir+"/diffpic3d.nii.gz";
    std::string pointSetPath = outDir + "/diffpointset.mps";
    std::string surfacePath = outDir + "/diffsurface.stl";
    std::string pointSetPathWithDefaultExtension = outDir + "/diffpointset2.mps";
    std::string pointSetPathWithoutDefaultExtension = outDir + "/diffpointset2.xXx";

    // the cases where no exception should be thrown
    MITK_TEST_CONDITION(mitk::IOUtil::SaveImage(img1, imagePath), "Testing if the image could be saved");
    MITK_TEST_CONDITION(mitk::IOUtil::SaveBaseData(img1.GetPointer(), imagePath2), "Testing if the image could be saved");
    MITK_TEST_CONDITION(mitk::IOUtil::SavePointSet(pointset, pointSetPath), "Testing if the pointset could be saved");
    MITK_TEST_CONDITION(mitk::IOUtil::SaveSurface(surface, surfacePath), "Testing if the surface could be saved");

    // test if defaultextension is inserted if no extension is present
    MITK_TEST_CONDITION(mitk::IOUtil::SavePointSet(pointset, pointSetPathWithoutDefaultExtension.c_str()), "Testing if the pointset could be saved");

    // test if exception is thrown as expected on unknown extsension
    MITK_TEST_FOR_EXCEPTION(mitk::Exception, mitk::IOUtil::SaveSurface(surface,"testSurface.xXx"));
    //load data which does not exist
    MITK_TEST_FOR_EXCEPTION(mitk::Exception, mitk::IOUtil::LoadImage("fileWhichDoesNotExist.nrrd"));

    //delete the files after the test is done
    remove(imagePath.c_str());
    remove(pointSetPath.c_str());
    remove(surfacePath.c_str());
    //remove the pointset with default extension and not the one without
    remove(pointSetPathWithDefaultExtension.c_str());

    mitk::Image::Pointer relativImage = mitk::ImageGenerator::GenerateGradientImage<float>(4,4,4,1);
    mitk::IOUtil::SaveImage(relativImage, "tempfile.nrrd");
    try
    {
      mitk::IOUtil::LoadImage("tempfile.nrrd");
      MITK_TEST_CONDITION(true, "Temporary image is in right place");
      remove("tempfile.nrrd");
    }
    catch (mitk::Exception &e)
    {
      MITK_TEST_CONDITION(false, "Temporary image is in right place");
    }

    MITK_TEST_END();
}