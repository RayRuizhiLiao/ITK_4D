project(ITKv4ImageRegistration)

add_executable(TemporalRegistration_Rigid_Final_TemporalReg  TemporalRegistration_Rigid_Final_TemporalReg.cxx )
target_link_libraries(TemporalRegistration_Rigid_Final_TemporalReg ${ITK_LIBRARIES})

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
