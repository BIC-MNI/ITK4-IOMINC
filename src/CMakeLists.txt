set(ITKIOMINC_SRC
  itkMINCImageIO.cxx
  itkMINCImageIOFactory.cxx
)

add_library(ITKIOMINC ${ITKIOMINC_SRC})
target_link_libraries(ITKIOMINC  ${ITKMINC_LIBRARIES} ${ITKIOImageBase_LIBRARIES})
itk_module_target(ITKIOMINC)
