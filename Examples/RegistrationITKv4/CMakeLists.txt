project(ITKv4ImageRegistration)

add_executable(TemporalRegistration_Rigid_Final_TemporalRegRotation  TemporalRegistration_Rigid_Final_TemporalRegRotation.cxx )
target_link_libraries(TemporalRegistration_Rigid_Final_TemporalRegRotation ${ITK_LIBRARIES})

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
