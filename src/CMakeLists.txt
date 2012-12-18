set(ITKIOMINC2_SRC
  itkMINC2ImageIO.cxx
  itkMINC2ImageIOFactory.cxx
)

add_library(ITKIOMINC2 ${ITKIOMINC2_SRC})
target_link_libraries(ITKIOMINC2  ${ITKMINC2_LIBRARIES} ${ITKIOImageBase_LIBRARIES})
itk_module_target(ITKIOMINC2)
