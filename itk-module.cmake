set(DOCUMENTATION "This module contains classes for reading and writing image
files in the MINC2 file Format (mnc).")

itk_module(ITKIOMINC2
  DEPENDS
    ITKMINC2
    ITKIOImageBase
  TEST_DEPENDS
    ITKTestKernel
  DESCRIPTION
    "${DOCUMENTATION}"
)
