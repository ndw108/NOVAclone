cmake_minimum_required (VERSION 3.6)

IF(DEFINED CMAKE_BUILD_TYPE)
   SET(CMAKE_BUILD_TYPE ${CMAKE_BUILD_TYPE} CACHE STRING "Choose the type of
build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug
Release RelWithDebInfo MinSizeRel.")
ELSE()
   SET(CMAKE_BUILD_TYPE RELEASE CACHE STRING "Choose the type of build,
options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release
RelWithDebInfo MinSizeRel.")
ENDIF()

set(CMAKE_CXX_FLAGS_RELEASE "-O3")

project (NOVA)

if (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set (CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}/NOVA" CACHE PATH "default install path" FORCE )
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/findFFTW")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/findFFTMPI")
set(LIBS)

set(CMAKE_CXX_STANDARD 20) # C++20...
set(CMAKE_CXX_STANDARD_REQUIRED ON) #...is required...
set(CMAKE_CXX_EXTENSIONS OFF) #...without compiler extensions like gnu++11

set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost 1.61 COMPONENTS mpi timer system) 



include_directories(${Boost_INCLUDE_DIRS}) 

option(USE_MPI "Use MPI" ON)
option(USE_OMP "Use OpenMP" OFF)
option(USE_FFTMPI "Use FFTMPI" ON)

###mpi
if(USE_MPI)
  add_definitions(-DHAVE_MPI)
  find_package(MPI REQUIRED)
  if (MPI_FOUND)
    include_directories(SYSTEM ${MPI_INCLUDE_PATH})
  else (MPI_FOUND)
    message(SEND_ERROR "MPI not found")
  endif (MPI_FOUND)
endif()

###mui
option(USE_MUI "Use MUI" OFF)
if(USE_MUI)
  add_definitions(-DHAVE_MUI)
  include_directories(../)
endif()

###openmp
if(USE_OMP)
  find_package(OpenMP)
  message(SEND_ERROR "OMP not yet supported / experimental.")
endif()

if(USE_FFTMPI)
  find_package(FFTW REQUIRED)
  find_package(FFTMPI REQUIRED)
  add_definitions(-DHAVE_FFTMPI)
  include_directories(${FFTW_INCLUDE_DIRS})
  include_directories(${FFTMPI_INCLUDES})
  set(LIBS ${LIBS} "${FFTMPI_LIBRARIES}" "${FFTW_LIBRARIES}" )
endif()

option(PROFILE "Use gperftools" OFF)
if(PROFILE)
  add_definitions(-DPROF)
  find_package(Gperftools)
  include_directories(${GPERFTOOLS_INCLUDE_DIR})
  set(LIBS ${LIBS} "${GPERFTOOLS_LIBRARIES}")
endif()

include_directories(src)

add_subdirectory (src)
add_subdirectory (app)

