project(TemporalRegistration)

add_executable(TemporalRegistration_Rigid  TemporalRegistration_Rigid.cxx )
target_link_libraries(TemporalRegistration_Rigid ${ITK_LIBRARIES})

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
