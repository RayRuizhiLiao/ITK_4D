project(ITKv4ImageRegistration)

add_executable(TemporalRegistration_Rigid_Final_noInitialization   TemporalRegistration_Rigid_Final_noInitialization.cxx )
target_link_libraries(TemporalRegistration_Rigid_Final_noInitialization ${ITK_LIBRARIES})

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
