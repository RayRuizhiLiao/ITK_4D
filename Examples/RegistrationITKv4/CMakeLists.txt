project(ITKv4ImageRegistration)

add_executable(TemporalRegistration_BSpline_Final_noChange  TemporalRegistration_BSpline_Final_noChange.cxx )
target_link_libraries(TemporalRegistration_BSpline_Final_noChange ${ITK_LIBRARIES})

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
