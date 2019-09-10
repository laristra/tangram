#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

project(tangram)

cinch_minimum_required(VERSION 1.0)


# SEMANTIC VERSION NUMBERS - UPDATE DILIGENTLY
# As soon as a change with a new version number is merged into the master,
# tag the central repository.

set(TANGRAM_VERSION_MAJOR 0)
set(TANGRAM_VERSION_MINOR 9)
set(TANGRAM_VERSION_PATCH 6)


# If a C++14 compiler is available, then set the appropriate flags
include(cxx14)
check_for_cxx14_compiler(CXX14_COMPILER)
if(CXX14_COMPILER)
  enable_cxx14()
else()
  message(STATUS "C++14 compatible compiler not found")
endif()

# If we couldn't find a C++14 compiler, try to see if a C++11
# compiler is available, then set the appropriate flags
if (NOT CXX14_COMPILER)
  include(cxx11)
  check_for_cxx11_compiler(CXX11_COMPILER)
  if(CXX11_COMPILER)
    enable_cxx11()
  else()
    message(FATAL_ERROR "C++11 compatible compiler not found")
  endif()
endif()

# cinch extras

cinch_load_extras()

set(CINCH_HEADER_SUFFIXES "\\.h")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")

# set the name of the Portage library

set(TANGRAM_LIBRARY "tangram" CACHE STRING "Name of the tangram library")


#-----------------------------------------------------------------------------
# Gather all the third party libraries needed for Tangram
#-----------------------------------------------------------------------------

set(ARCHOS ${CMAKE_SYSTEM_PROCESSOR}_${CMAKE_SYSTEM_NAME})

set(TANGRAM_EXTRA_LIBRARIES)

#-----------------------------------------------------------------------------
# Wonton
#-----------------------------------------------------------------------------
if (WONTON_DIR)

  # Link with an existing installation of Wonton, if provided. 
  find_package(WONTON REQUIRED)
  message(STATUS "WONTON_LIBRARIES=${WONTON_LIBRARIES}" )
  include_directories(${WONTON_INCLUDE_DIR})
  message(STATUS "WONTON_INCLUDE_DIRS=${WONTON_INCLUDE_DIR}")
 
  list(APPEND TANGRAM_EXTRA_LIBRARIES ${WONTON_LIBRARIES})

else (WONTON_DIR)

  # Build Wonton from a submodule
  file(GLOB _wonton_contents ${CMAKE_SOURCE_DIR}/wonton/*)
  if (_wonton_contents)
    if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
      # We are building tangram, and wonton is a subdirectory
      add_subdirectory(${CMAKE_SOURCE_DIR}/wonton)
    endif()
    include_directories(${CMAKE_SOURCE_DIR}/wonton)
    include_directories(${CMAKE_BINARY_DIR}/wonton)  # for wonton-config.h
    list(APPEND TANGRAM_EXTRA LIBRARIES wonton)
    set(WONTON_FOUND True)

    # If Wonton is included as a submodule, it will get installed alongside Portage
    set(WONTON_DIR ${CMAKE_INSTALL_PREFIX})
  endif()
endif (WONTON_DIR)
  
if (NOT WONTON_FOUND)
  message(FATAL_ERROR "WONTON_DIR is not specified and Wonton is not a subdirectory !")
endif() 


#------------------------------------------------------------------------------#
# If we are building with FleCSI, then we need a modern C++ compiler
#------------------------------------------------------------------------------#
if(ENABLE_FleCSI)
  # we already checked for CXX14 in project.cmake
  if(NOT CXX14_COMPILER)
    message(STATUS "C++14 compatible compiler not found")
  endif()
endif()

#------------------------------------------------------------------------------#
# Set up MPI builds
# (eventually most of this should be pushed down into cinch)
#------------------------------------------------------------------------------#
set(ENABLE_MPI OFF CACHE BOOL "")
if (ENABLE_MPI)
  find_package(MPI REQUIRED)
  set(TANGRAM_ENABLE_MPI True CACHE BOOL "Is Tangram compiled with MPI?")
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER} CACHE FILEPATH "C compiler to use" FORCE)
  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER} CACHE FILEPATH "C++ compiler to use" FORCE)
endif ()


set(ARCHOS ${CMAKE_SYSTEM_PROCESSOR}_${CMAKE_SYSTEM_NAME})

#-----------------------------------------------------------------------------
# FleCSI and FleCSI-SP location
#-----------------------------------------------------------------------------

set(ENABLE_FleCSI FALSE CACHE BOOL "Use FleCSI")
if (ENABLE_FleCSI)
 
 find_package(FleCSI REQUIRED)
 message(STATUS "FleCSI_LIBRARIES=${FleCSI_LIBRARIES}" )
 include_directories(${FleCSI_INCLUDE_DIR})
 message(STATUS "FleCSI_INCLUDE_DIRS=${FleCSI_INCLUDE_DIR}")
 list(APPEND TANGRAM_EXTRA_LIBRARIES ${FleCSI_LIBRARIES})

 find_package(FleCSISP REQUIRED)
 message(STATUS "FleCSISP_LIBRARIES=${FleCSISP_LIBRARIES}" )
 include_directories(${FleCSISP_INCLUDE_DIR})
 message(STATUS "FleCSISP_INCLUDE_DIRS=${FleCSISP_INCLUDE_DIR}")
 list(APPEND TANGRAM_EXTRA_LIBRARIES ${FleCSISP_LIBRARIES})

  ######################################################################
  # This is a placeholder for how we would do IO with FleCSI
  # There are still some issues with dumping the targetMesh data
  #
  # WARNING!!! THIS IS POTENTIALLY FRAGILE
  # it appears to work, but could cause conflicts with EXODUS and
  # other libraries used by Jali
  #
  # FOR NOW THIS IS DISABLED UNTIL WE CAN GET A PROPER WORKAROUND
  ######################################################################
  # STRING(REPLACE "flecsi" "flecsi-tpl" FLECSI_TPL_DIR ${FLECSI_INSTALL_DIR})
  # message(STATUS "FLECSI_TPL_DIR=${FLECSI_TPL_DIR}")
  # if(IS_DIRECTORY ${FLECSI_TPL_DIR})
  #   find_library(EXODUS_LIBRARY
  #     NAMES exodus
  #     PATHS ${FLECSI_TPL_DIR}
  #     PATH_SUFFIXES lib
  #     NO_DEFAULT_PATH)
  #   find_path(EXODUS_INCLUDE_DIR
  #     NAMES exodusII.h
  #     PATHS ${FLECSI_TPL_DIR}
  #     PATH_SUFFIXES include
  #     NO_DEFAULT_PATH)

  #   if(EXODUS_LIBRARY AND EXODUS_INCLUDE_DIR)
  #     set(FLECSI_LIBRARIES ${EXODUS_LIBRARY} ${FLECSI_LIBRARIES})
  #     include_directories(${EXODUS_INCLUDE_DIR})
  #     add_definitions(-DHAVE_EXODUS)
  #   endif(EXODUS_LIBRARY AND EXODUS_INCLUDE_DIR)

  # endif(IS_DIRECTORY ${FLECSI_TPL_DIR})
endif()



#------------------------------------------------------------------------------#
# Configure Jali
# (this includes the TPLs that Jali will need)
#------------------------------------------------------------------------------#

set(ENABLE_JALI OFF CACHE BOOL "")
if (JALI_DIR)  # forgive users for capitalization mistake
  set(Jali_DIR ${JALI_DIR})
endif (JALI_DIR)
if (Jali_DIR)

   # Look for the Jali package

   find_package(Jali REQUIRED
                HINTS ${Jali_DIR}/lib)

   message(STATUS "Located Jali")
   message(STATUS "Jali_DIR=${Jali_DIR}")

   # add full path to jali libs
   unset(_LIBS)
   foreach (_lib ${Jali_LIBRARIES})
      set(_LIBS ${_LIBS};${Jali_LIBRARY_DIRS}/lib${_lib}.a)
   endforeach()
   set(Jali_LIBRARIES ${_LIBS})

   include_directories(${Jali_INCLUDE_DIRS} ${Jali_TPL_INCLUDE_DIRS})

   add_definitions(-DENABLE_JALI)

   list(APPEND TANGRAM_EXTRA_LIBRARIES ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES})
endif (Jali_DIR)

#------------------------------------------------------------------------------#
# Configure XMOF2D
#------------------------------------------------------------------------------#

if (XMOF2D_DIR)

   # Look for the XMOF2D package

   find_package(XMOF2D REQUIRED
                HINTS ${XMOF2D_DIR}/lib)

   message(STATUS "Located XMOF2D")
   message(STATUS "XMOF2D_DIR=${XMOF2D_DIR}")

   # add full path to XMOF2D libs
   unset(_LIBS)
   foreach (_lib ${XMOF2D_LIBRARIES})
      set(_LIBS ${_LIBS};${XMOF2D_LIBRARY_DIRS}/lib${_lib}.a)
   endforeach()
   set(XMOF2D_LIBRARIES ${_LIBS})

   # message(STATUS "XMOF2D_INCLUDE_DIRS=${XMOF2D_INCLUDE_DIRS}")
   # message(STATUS "XMOF2D_LIBRARY_DIRS=${XMOF2D_LIBRARY_DIRS}")
   # message(STATUS "XMOF2D_LIBRARIES=${XMOF2D_LIBRARIES}")

   include_directories(${XMOF2D_INCLUDE_DIRS})
   list(APPEND PORTAGE_EXTRA_LIBRARIES ${XMOF2D_LIBRARIES})

endif (XMOF2D_DIR)

#-----------------------------------------------------------------------------
# General NGC include directory information
#-----------------------------------------------------------------------------
set(NGC_INCLUDE_DIR "$ENV{NGC_INCLUDE_DIR}" CACHE PATH "NGC include directory")
if(NGC_INCLUDE_DIR)
  message(STATUS "Using NGC_INCLUDE_DIR=${NGC_INCLUDE_DIR}")
endif(NGC_INCLUDE_DIR)

#-----------------------------------------------------------------------------
# Thrust information
#-----------------------------------------------------------------------------
set(ENABLE_THRUST FALSE CACHE BOOL "Use Thrust")
if(ENABLE_THRUST)
  message(STATUS "Enabling compilation with Thrust")

  set(TANGRAM_ENABLE_THRUST True CACHE BOOL "Is Tangram compiled with Thrust?")

  # allow the user to specify a THRUST_DIR, otherwise use ${NGC_INCLUDE_DIR}
  # NOTE: thrust internally uses include paths from the 'root' directory, e.g.
  #
  #       #include "thrust/device_vector.h"
  #
  #       so the path here should point to the directory that has thrust as
  #       a subdirectory.
  # Use THRUST_DIR directly if specified, otherwise try to build from NGC
  set(THRUST_DIR "${NGC_INCLUDE_DIR}" CACHE PATH "Thrust directory")
  message(STATUS "Using THRUST_DIR=${THRUST_DIR}")

  # Allow for swapping backends - should this be in CACHE?
  set(THRUST_BACKEND "THRUST_DEVICE_SYSTEM_OMP" CACHE STRING "Thrust backend")
  message(STATUS "Using ${THRUST_BACKEND} as Thrust backend.")
  include_directories(${THRUST_DIR})
  add_definitions(-DTHRUST)
  add_definitions(-DTHRUST_DEVICE_SYSTEM=${THRUST_BACKEND})

  if("${THRUST_BACKEND}" STREQUAL "THRUST_DEVICE_SYSTEM_OMP")
    FIND_PACKAGE( OpenMP REQUIRED)
    if(OPENMP_FOUND)
      set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif(OPENMP_FOUND)
  endif ()

  if("${THRUST_BACKEND}" STREQUAL "THRUST_DEVICE_SYSTEM_TBB")
    FIND_PACKAGE(TBB REQUIRED)
    if(TBB_FOUND)
      include_directories(${TBB_INCLUDE_DIRS})
      link_directories(${TBB_LIBRARY_DIRS})
      set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -ltbb")
    endif(TBB_FOUND)
  endif()

else (ENABLE_THRUST)

  find_package(Boost REQUIRED)
  if(Boost_FOUND)
   message(STATUS "Boost location: ${Boost_INCLUDE_DIRS}")
   include_directories( ${Boost_INCLUDE_DIRS} )
  endif()

endif(ENABLE_THRUST)


#-----------------------------------------------------------------------------
# Now add the source directories and library targets
#-----------------------------------------------------------------------------

# In addition to the include directories of the source set by cinch,
# we need to include the build directory to get the autogenerated
# wonton-config.h

include_directories(${CMAKE_BINARY_DIRECTORY})

# Apps and Libraries
cinch_add_application_directory(app)
cinch_add_library_target(tangram tangram)

# Add application tests
# May pull this logic into cinch at some future point
option(ENABLE_APP_TESTS "Enable testing of full app" OFF)
if(ENABLE_APP_TESTS)
  enable_testing()
endif()

#------------------------------------------------------------------------------#
#
#------------------------------------------------------------------------------#

# retrieve all the definitions we added for compiling
get_directory_property(TANGRAM_COMPILE_DEFINITIONS DIRECTORY ${CMAKE_SOURCE_DIR} COMPILE_DEFINITIONS)

# build the TANGRAM_LIBRARIES variable
set(TANGRAM_LIBRARIES ${TANGRAM_LIBRARY} ${TANGRAM_EXTRA_LIBRARIES} CACHE STRING "List of libraries to link with tangram")

############################################################################## 
# Write a configuration file from template replacing only variables enclosed
# by the @ sign. This will let other programs build on TANGRAM discover how
# TANGRAM was built and which TPLs it used
#############################################################################

configure_file(${PROJECT_SOURCE_DIR}/cmake/tangram_config.cmake.in 
               ${PROJECT_BINARY_DIR}/tangram_config.cmake @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/tangram_config.cmake 
        DESTINATION ${CMAKE_INSTALL_PREFIX}/share/cmake/)

configure_file(${PROJECT_SOURCE_DIR}/config/tangram-config.h.in
               ${PROJECT_BINARY_DIR}/tangram-config.h @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/tangram-config.h
        DESTINATION ${CMAKE_INSTALL_PREFIX}/include/)

