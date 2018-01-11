project(ITKv4ImageRegistration)


add_executable(TemporalRegistration_mask_BSpline   TemporalRegistration_mask_BSpline.cxx )
target_link_libraries(TemporalRegistration_mask_BSpline ${ITK_LIBRARIES})


add_executable(itkInverseDisplacementFieldImageFilterTest   itkInverseDisplacementFieldImageFilterTest.cxx )
target_link_libraries(itkInverseDisplacementFieldImageFilterTest ${ITK_LIBRARIES})

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
