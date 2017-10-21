project(ITKv4ImageRegistration)


add_executable(TemporalRegistration5   TemporalRegistration5.cxx )
target_link_libraries(TemporalRegistration5 ${ITK_LIBRARIES})

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
