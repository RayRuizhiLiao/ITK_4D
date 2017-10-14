project(ITKv4ImageRegistration)


add_executable(TemporalRegistration4   TemporalRegistration4.cxx )
target_link_libraries(TemporalRegistration4 ${ITK_LIBRARIES})

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
