project(ITKv4ImageRegistration)

add_executable(TemporalRegistration_BSpline_Final_fixTemplate   TemporalRegistration_BSpline_Final_fixTemplate.cxx )
target_link_libraries(TemporalRegistration_BSpline_Final_fixTemplate ${ITK_LIBRARIES})

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
