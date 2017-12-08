project(ITKv4ImageRegistration)


add_executable(TemporalRegistration_mask_rigid_maskcenter_noOptimization   TemporalRegistration_mask_rigid_maskcenter_noOptimization.cxx )
target_link_libraries(TemporalRegistration_mask_rigid_maskcenter_noOptimization ${ITK_LIBRARIES})

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
