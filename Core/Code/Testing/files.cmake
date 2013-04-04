
# tests with no extra command line parameter
set(MODULE_TESTS
  mitkAccessByItkTest.cpp
  mitkCoreObjectFactoryTest.cpp
  mitkMaterialTest.cpp
  mitkActionTest.cpp
  mitkDispatcherTest.cpp
  mitkEnumerationPropertyTest.cpp
  mitkEventTest.cpp
  #mitkEventConfigTest.cpp ## needs to be re-written, test indirect since EventConfig is no longer exported as interface Bug 14529
  mitkFocusManagerTest.cpp
  mitkGenericPropertyTest.cpp
  mitkGeometry3DTest.cpp
  mitkGeometryDataToSurfaceFilterTest.cpp
  mitkGlobalInteractionTest.cpp
  mitkImageDataItemTest.cpp
  #mitkImageMapper2DTest.cpp
  mitkImageGeneratorTest.cpp
  mitkBaseDataTest.cpp
  #mitkImageToItkTest.cpp
  mitkInstantiateAccessFunctionTest.cpp
  mitkInteractorTest.cpp
  mitkInteractionEventTest.cpp
  mitkITKThreadingTest.cpp
  mitkLevelWindowTest.cpp
  mitkMessageTest.cpp
  #mitkPipelineSmartPointerCorrectnessTest.cpp
  mitkPixelTypeTest.cpp
  mitkPlaneGeometryTest.cpp
  mitkPointSetFileIOTest.cpp
  mitkPointSetTest.cpp
  mitkPointSetWriterTest.cpp
  mitkPointSetReaderTest.cpp
  mitkPointSetInteractorTest.cpp
  mitkPropertyTest.cpp
  mitkPropertyListTest.cpp
  #mitkRegistrationBaseTest.cpp
  #mitkSegmentationInterpolationTest.cpp
  mitkSlicedGeometry3DTest.cpp
  mitkSliceNavigationControllerTest.cpp
  mitkStateMachineTest.cpp
  ##mitkStateMachineContainerTest.cpp ## rewrite test, indirect since no longer exported Bug 14529
  mitkStateTest.cpp
  mitkSurfaceTest.cpp
  mitkSurfaceToSurfaceFilterTest.cpp
  mitkTimeSlicedGeometryTest.cpp
  mitkTransitionTest.cpp
  mitkUndoControllerTest.cpp
  mitkVtkWidgetRenderingTest.cpp
  mitkVerboseLimitedLinearUndoTest.cpp
  mitkWeakPointerTest.cpp
  mitkTransferFunctionTest.cpp
  #mitkAbstractTransformGeometryTest.cpp
  mitkStepperTest.cpp
  itkTotalVariationDenoisingImageFilterTest.cpp
  mitkRenderingManagerTest.cpp
  vtkMitkThickSlicesFilterTest.cpp
  mitkNodePredicateSourceTest.cpp
  mitkVectorTest.cpp
  mitkClippedSurfaceBoundsCalculatorTest.cpp
  #QmitkRenderingTestHelper.cpp
  mitkExceptionTest.cpp
  mitkExtractSliceFilterTest.cpp
  mitkLogTest.cpp
  mitkImageDimensionConverterTest.cpp
  mitkLoggingAdapterTest.cpp
  mitkUIDGeneratorTest.cpp
  mitkShaderRepositoryTest.cpp
)

# test with image filename as an extra command line parameter
set(MODULE_IMAGE_TESTS
  mitkPlanePositionManagerTest.cpp
  mitkSurfaceVtkWriterTest.cpp
  #mitkImageSliceSelectorTest.cpp
  mitkImageTimeSelectorTest.cpp
  # mitkVtkPropRendererTest.cpp
  mitkDataNodeFactoryTest.cpp
  #mitkSTLFileReaderTest.cpp
  mitkImageAccessorTest.cpp
)

# list of images for which the tests are run
set(MODULE_TESTIMAGES
 # Pic-Factory no more available in Core, test images now in .nrrd format
  US4DCyl.nrrd
  Pic3D.nrrd
  Pic2DplusT.nrrd
  BallBinary30x30x30.nrrd
  binary.stl
  ball.stl
)

set(MODULE_CUSTOM_TESTS
    #mitkLabeledImageToSurfaceFilterTest.cpp
    #mitkExternalToolsTest.cpp
    mitkDataStorageTest.cpp
    mitkDataNodeTest.cpp
    mitkDicomSeriesReaderTest.cpp
    mitkDICOMLocaleTest.cpp
    mitkEventMapperTest.cpp
    mitkNodeDependentPointSetInteractorTest.cpp
    mitkStateMachineFactoryTest.cpp
    mitkPointSetLocaleTest.cpp
    mitkImageTest.cpp
    mitkImageWriterTest.cpp
    mitkImageVtkMapper2DTest.cpp
    mitkImageVtkMapper2DLevelWindowTest.cpp
    mitkImageVtkMapper2DOpacityTest.cpp
    mitkImageVtkMapper2DColorTest.cpp
    mitkImageVtkMapper2DSwivelTest.cpp
    mitkImageVtkMapper2DTransferFunctionTest.cpp
    mitkIOUtilTest.cpp
    mitkSurfaceVtkMapper3DTest
    mitkSurfaceVtkMapper3DTexturedSphereTest.cpp
    mitkVolumeCalculatorTest.cpp
    mitkLevelWindowManagerTest.cpp
)

set(MODULE_RESOURCE_FILES
  Interactions/AddAndRemovePoints.xml
  Interactions/globalConfig.xml
  Interactions/StatemachineTest.xml
  Interactions/StatemachineConfigTest.xml
)

# Create an artificial module initializing class for
# the usServiceListenerTest.cpp
usFunctionGenerateModuleInit(testdriver_init_file
                             NAME ${MODULE_NAME}TestDriver
                             DEPENDS "Mitk"
                             VERSION "0.1.0"
                             EXECUTABLE
                            )

# Embed the resources
set(testdriver_resources )
usFunctionEmbedResources(testdriver_resources
                         EXECUTABLE_NAME ${MODULE_NAME}TestDriver
                         ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Resources
                         FILES ${MODULE_RESOURCE_FILES}
                        )

set(TEST_CPP_FILES ${testdriver_init_file} ${testdriver_resources})
