project(ITKv4ImageRegistration)


add_executable(TemporalRegistration_mask_rigid   TemporalRegistration_mask_rigid.cxx )
target_link_libraries(TemporalRegistration_mask_rigid ${ITK_LIBRARIES})

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
