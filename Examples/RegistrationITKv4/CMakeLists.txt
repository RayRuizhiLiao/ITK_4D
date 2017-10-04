project(ITKv4ImageRegistration)


add_executable(TemporalRegistration   TemporalRegistration.cxx )
target_link_libraries(TemporalRegistration ${ITK_LIBRARIES})

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
