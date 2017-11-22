project(ITKv4ImageRegistration)


add_executable(TemporalRegistration_mask_rigid_test4   TemporalRegistration_mask_rigid_test4.cxx )
target_link_libraries(TemporalRegistration_mask_rigid_test4 ${ITK_LIBRARIES})

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
