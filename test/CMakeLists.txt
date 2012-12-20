itk_module_test()

IF(NOT ITK_USE_SYSTEM_HDF5)
  include_directories(
    "${ITKHDF5_SOURCE_DIR}/src/itkhdf5"
    "${ITKHDF5_BINARY_DIR}/src/itkhdf5"
    )
ENDIF()

set(ITKIOMINCTests
   itkMINCImageIOTest.cxx
   itkMINCImageIOTest2.cxx
   itkMINCImageIOTest3.cxx
  )

# set(ITKIOTransformMincTests
#    itkIOTransformMincTest.cxx
# )

CreateTestDriver(ITKIOMINC "${ITKIOMINC-Test_LIBRARIES}" "${ITKIOMINCTests}")

itk_add_test(NAME itkMINCImageIOTest
  COMMAND ITKIOMINCTestDriver itkMINCImageIOTest ${ITK_TEST_OUTPUT_DIR} )

itk_add_test(NAME itkMINCImageIOTest2
  COMMAND ITKIOMINCTestDriver itkMINCImageIOTest2
      DATA{Input/pd_z+_byte_cor.mnc} ${ITK_TEST_OUTPUT_DIR}/pd_z+_byte_cor.mnc)

itk_add_test(NAME itkMINCImageIOTest3
  COMMAND ITKIOMINCTestDriver itkMINCImageIOTest3
      DATA{Input/pd_z+_byte_cor.mnc} ${ITK_TEST_OUTPUT_DIR}/pd_z+_byte_cor_2.mnc)

# CreateTestDriver(ITKIOTransformMINC "${ITKIOMINC-Test_LIBRARIES}" "${ITKIOTransformMincTests}")
# TARGET_LINK_LIBRARIES(ITKIOTransformMINCTestDriver ITKIOTransformBase )
#
# itk_add_test(NAME itkIOTransformMINCTest
#       COMMAND ITKIOTransformMINCTestDriver itkIOTransformMincTest)

