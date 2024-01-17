# - Find the PFFT library
#
# Usage:
#   find_package(PFFT [REQUIRED] [QUIET] )
#     
# It sets the following variables:
#   PFFT_FOUND               ... true if PFFT is found on the system
#   PFFT_LIBRARIES           ... full path to PFFT library
#   PFFT_INCLUDES            ... PFFT include directory
#
# The following variables will be checked by the function
#   PFFT_USE_STATIC_LIBS    ... if true, only static libraries are found
#   PFFT_ROOT               ... if set, the libraries are exclusively searched
#                               under this path
#   PFFT_LIBRARY            ... PFFT library to use
#   PFFT_INCLUDE_DIR        ... PFFT include directory
#

#If environment variable PFFTDIR is specified, it has same effect as PFFT_ROOT
if( NOT FFTMPI_ROOT AND ENV{FFTMPIDIR} )
  set( FFTMPI_ROOT $ENV{FFTMPIDIR} )
endif()

# Check if we can use PkgConfig
find_package(PkgConfig)

#Determine from PKG
if( PKG_CONFIG_FOUND AND NOT FFTMPI_ROOT )
  pkg_check_modules( PKG_FFTMPI "fftmpi" QUIET )
endif()

#Check whether to search static or dynamic libs
set( CMAKE_FIND_LIBRARY_SUFFIXES_SAV ${CMAKE_FIND_LIBRARY_SUFFIXES} )

if( ${FFTMPI_USE_STATIC_LIBS} )
  set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_STATIC_LIBRARY_SUFFIX} )
else()
  set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_SHARED_LIBRARY_SUFFIX} )
endif()

#find libs
find_library(
  FFTMPI_LIBRARIES
  NAMES "libfft3dmpi"
  PATHS ${FFTMPI_ROOT} ${PKG_FFTMPI_PREFIX} ${PKG_FFTMPI_LIBRARY_DIRS} ${LIB_INSTALL_DIR}
  PATH_SUFFIXES "lib" "lib64" "so"
)

#find includes
find_path(
  FFTMPI_INCLUDES
  NAMES "fft3d.h"
  PATHS ${FFTMPI_ROOT} ${PKG_FFTMPI_PREFIX} ${PKG_FFTMPI_INCLUDE_DIRS} ${INCLUDE_INSTALL_DIR}
  PATH_SUFFIXES "include"
)

set( CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES_SAV} )

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTMPI DEFAULT_MSG
                                  FFTMPI_INCLUDES FFTMPI_LIBRARIES)

mark_as_advanced(FFTMPI_INCLUDES FFTMPI_LIBRARIES)

