project(ITKv4ImageRegistration)

add_executable(TemporalRegistration_BSpline_Final_TemporalReg_TwoLevels  TemporalRegistration_BSpline_Final_TemporalReg_TwoLevels.cxx )
target_link_libraries(TemporalRegistration_BSpline_Final_TemporalReg_TwoLevels ${ITK_LIBRARIES})

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
