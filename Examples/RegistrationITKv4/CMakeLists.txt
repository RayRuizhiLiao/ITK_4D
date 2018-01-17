project(ITKv4ImageRegistration)

add_executable(TemporalRegistration_mask_rigid_maskcenter   TemporalRegistration_mask_rigid_maskcenter.cxx )
target_link_libraries(TemporalRegistration_mask_rigid_maskcenter ${ITK_LIBRARIES})

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
