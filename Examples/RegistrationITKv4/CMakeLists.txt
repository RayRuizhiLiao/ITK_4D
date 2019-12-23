project(ITKv4ImageRegistration)

add_executable(DeformableRegistration7_12GridNodes  DeformableRegistration7_12GridNodes.cxx )
target_link_libraries(DeformableRegistration7_12GridNodes ${ITK_LIBRARIES})

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
