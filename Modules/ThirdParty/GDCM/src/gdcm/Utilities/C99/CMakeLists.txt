# C99 install rules:
#if(WIN32 AND NOT CYGWIN)

if(NOT GDCM_INSTALL_NO_DEVELOPMENT)
  install(FILES stdint.h
    DESTINATION ${GDCM_INSTALL_INCLUDE_DIR} COMPONENT Headers
  )
endif()
