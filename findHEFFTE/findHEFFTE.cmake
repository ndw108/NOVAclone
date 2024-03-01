# - Find heFFTe 


#find libs
find_library(
  HEFFTE_LIBRARIES
  NAMES "libheffte"
  PATHS  ${CMAKE_INSTALL_PREFIX} 
  PATH_SUFFIXES "lib" "lib64"
)

#find includes
find_path(
  HEFFTE_INCLUDES
  NAMES "heffte.h"
  PATHS ${CMAKE_INSTALL_PREFIX} 
  PATH_SUFFIXES "include"
)


set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_SAV} )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Heffte DEFAULT_MSG
                                  HEFFTE_INCLUDES HEFFTE_LIBRARIES)

mark_as_advanced(HEFFTE_INCLUDES HEFFTE_LIBRARIES)
