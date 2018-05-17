project(ITKv4ImageRegistration)

add_executable(TemporalRegistration_BSpline_Final_TemporalReg  TemporalRegistration_BSpline_Final_TemporalReg.cxx )
target_link_libraries(TemporalRegistration_BSpline_Final_TemporalReg ${ITK_LIBRARIES})

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
