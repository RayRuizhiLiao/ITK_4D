project(ITKOpenJPEG)
set(ITKOpenJPEG_THIRD_PARTY 1)

set(ITKOpenJPEG_INCLUDE_DIRS
  ${ITKOpenJPEG_BINARY_DIR}/src/openjpeg
  ${ITKOpenJPEG_SOURCE_DIR}/src/openjpeg
  )
set(ITKOpenJPEG_LIBRARIES itkopenjpeg)
set(ITK_LIBRARY_PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS FALSE)

itk_module_impl()

install(FILES ${ITKOpenJPEG_SOURCE_DIR}/src/itk_openjpeg.h
  DESTINATION ${ITKOpenJPEG_INSTALL_INCLUDE_DIR}
  COMPONENT Development
  )
