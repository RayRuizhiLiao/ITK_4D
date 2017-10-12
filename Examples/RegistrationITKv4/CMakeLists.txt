project(ITKv4ImageRegistration)


add_executable(TemporalRegistration2   TemporalRegistration2.cxx )
target_link_libraries(TemporalRegistration2 ${ITK_LIBRARIES})

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
