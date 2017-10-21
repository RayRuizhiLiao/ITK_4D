project(ITKv4ImageRegistration)


add_executable(TemporalRegistration5_newMetric   TemporalRegistration5_newMetric.cxx )
target_link_libraries(TemporalRegistration5_newMetric ${ITK_LIBRARIES})

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
