project(ITKv4ImageRegistration)


add_executable(TemporalRegistration5_mask_rigid   TemporalRegistration5_mask_rigid.cxx )
target_link_libraries(TemporalRegistration5_mask_rigid ${ITK_LIBRARIES})

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
