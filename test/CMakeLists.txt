itk_module_test()

set(ITKIOMINC2Tests
   itkMINC2ImageIOTest.cxx
)

# set(ITKIOTransformMincTests
#    itkIOTransformMincTest.cxx
# )

CreateTestDriver(ITKIOMINC2 "${ITKIOMINC2-Test_LIBRARIES}" "${ITKIOMINC2Tests}")

itk_add_test(NAME itkMINC2ImageIOTest
  COMMAND ITKIOMINC2TestDriver itkMINC2ImageIOTest ${ITK_TEST_OUTPUT_DIR} )
 
# CreateTestDriver(ITKIOTransformMINC "${ITKIOMINC-Test_LIBRARIES}" "${ITKIOTransformMincTests}")
# TARGET_LINK_LIBRARIES(ITKIOTransformMINCTestDriver ITKIOTransformBase )
# 
# itk_add_test(NAME itkIOTransformMINCTest
#       COMMAND ITKIOTransformMINCTestDriver itkIOTransformMincTest)

