project(ITKv4ImageRegistration)

add_executable(TemporalRegistration_Rigid_Final   TemporalRegistration_Rigid_Final.cxx )
target_link_libraries(TemporalRegistration_Rigid_Final ${ITK_LIBRARIES})

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
