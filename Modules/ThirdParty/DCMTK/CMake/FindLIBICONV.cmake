# Find iconv library
#
# Released under BSD license
#
#  LIBICONV_INCLUDE_DIRS - where to find iconv.h, etc
#  LIBICONV_LIBRARIES    - Lists of libraries when using iconv
#  LIBICONV_FOUND        - True if iconv found

INCLUDE(FindPackageHandleStandardArgs)

# Look for the header file
if(NOT LIBICONV_INCLUDE_DIR)
  FIND_PATH(LIBICONV_INCLUDE_DIR NAMES iconv.h)
endif()
MARK_AS_ADVANCED(LIBICONV_INCLUDE_DIR)

# Look for the library
SET(LIBICONV_LIBS iconv)

if(NOT LIBICONV_LIBRARY)
  FIND_LIBRARY(LIBICONV_LIBRARY NAMES ${LIBICONV_LIBS})
endif()
MARK_AS_ADVANCED(LIBICONV_LIBRARY)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(LIBICONV REQUIRED_VARS LIBICONV_LIBRARY LIBICONV_INCLUDE_DIR)

# Copy the result to output variables
IF(LIBICONV_FOUND)
  SET(LIBICONV_LIBRARIES ${LIBICONV_LIBRARY})
  SET(LIBICONV_INCLUDE_DIRS ${LIBICONV_INCLUDE_DIR})
  include_directories(${LIBICONV_INCLUDE_DIR})
  link_directories(${LIBICONV_LIBDIR})
ELSE(LIBICONV_FOUND)
  SET(LIBICONV_LIBS)
  SET(LIBICONV_LIBRARY)
  SET(LIBICONV_LIBRARIES)
  SET(LIBICONV_INCLUDE_DIR)
  SET(LIBICONV_INCLUDE_DIRS)
ENDIF(LIBICONV_FOUND)
