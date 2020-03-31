#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

cmake_minimum_required(VERSION 3.13)

project(tangram CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

cinch_minimum_required(VERSION 1.0)

if (CMAKE_VERSION_MAJOR GREATER_EQUAL 3.13)
  CMAKE_POLICY(SET CMP0079 NEW)  # allow target_link_libraries to reference
                                 # targets from other directories
endif()

cmake_policy(SET CMP0074 NEW)  # Don't ignore Pkg_ROOT variables



# SEMANTIC VERSION NUMBERS - UPDATE DILIGENTLY
# As soon as a change with a new version number is merged into the master,
# tag the central repository.

set(TANGRAM_VERSION_MAJOR 0)
set(TANGRAM_VERSION_MINOR 9)
set(TANGRAM_VERSION_PATCH 9)


# Top level target
add_library(tangram INTERFACE)

# Alias (Daniel Pfeiffer, Effective CMake) - this allows other
# projects that use Pkg as a subproject to find_package(Nmspc::Pkg)
# which does nothing because Pkg is already part of the project

add_library(tangram::tangram ALIAS tangram)
set(TANGRAM_LIBRARIES tangram::tangram)

# cinch extras

cinch_load_extras()

set(CINCH_HEADER_SUFFIXES "\\.h")

# Explicitly turn off module paths inherited from Cinch. When we get rid of
# Cinch we can delete this line
set(CMAKE_MODULE_PATH "")

# Find our modules first
if (CMAKE_VERSION GREATER_EQUAL 3.15)
  list(PREPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
else ()
  set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH}")
endif ()


#-----------------------------------------------------------------------------
# Find or compile Wonton
#-----------------------------------------------------------------------------
if (WONTON_ROOT)

  # Link with an existing installation of Wonton, if provided. 
  find_package(WONTON QUIET REQUIRED NAMES wonton)

  target_include_directories(tangram INTERFACE ${WONTON_INCLUDE_DIR})
  message(STATUS "WONTON_INCLUDE_DIR=${WONTON_INCLUDE_DIR}")

  target_link_libraries(tangram INTERFACE ${WONTON_LIBRARIES})
  message(STATUS "WONTON_LIBRARIES=${WONTON_LIBRARIES}" )

  if (ENABLE_THRUST AND NOT WONTON_ENABLE_THRUST)
    message(FATAL_ERROR "Thrust enabled for Tangram but Wonton is not built with Thrust")
  endif ()

  if (ENABLE_MPI AND NOT WONTON_ENABLE_MPI)
    message(FATAL_ERROR "MPI enabled for Tangram but Wonton is not compiled with MPI (WONTON_ENABLE_MPI=${WONTON_ENABLE_MPI})")
  endif ()

else ()

  # Build Wonton from a submodule
  file(GLOB _wonton_contents ${PROJECT_SOURCE_DIR}/wonton/*)
  if (_wonton_contents)
    if (CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
      # We are building tangram, and wonton is a subdirectory
      message(STATUS "Recursing down into wonton")
      add_subdirectory(${PROJECT_SOURCE_DIR}/wonton)
    endif()
    
    target_include_directories(tangram INTERFACE
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/wonton>
      $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/wonton>
      $<INSTALL_INTERFACE:include>)

    # set wonton_LIBRARIES which is needed for tangram_support library
    # but is set by the wonton build system only later (it would be
    # nice to avoid this)
#    set(WONTON_LIBRARIES wonton::wonton)
    target_link_libraries(tangram INTERFACE ${WONTON_LIBRARIES})
    set(WONTON_FOUND True)

    # If Wonton is included as a submodule, it will get installed alongside Portage
    set(WONTON_DIR ${CMAKE_INSTALL_PREFIX})
    set(WONTON_IS_SUBMODULE True)
  endif ()
  
endif ()
  
if (NOT WONTON_FOUND)
  message(FATAL_ERROR "WONTON_DIR is not specified and Wonton is not a subdirectory !")
endif() 


#-----------------------------------------------------------------------------
# Thrust information
#-----------------------------------------------------------------------------
set(ENABLE_THRUST FALSE CACHE BOOL "Use Thrust")  # create a cache var with default
if (ENABLE_THRUST)
  if (NOT WONTON_ENABLE_THRUST)
    message(FATAL_ERROR "Thrust enabled for Tangram but not for Wonton")
  endif ()
endif()


#------------------------------------------------------------------------------#
# Configure XMOF2D
#------------------------------------------------------------------------------#

if (ENABLE_XMOF2D)

  # Look for the XMOF2D package
  
  find_package(XMOF2D REQUIRED)
  
  message(STATUS "Located XMOF2D")
  
  target_include_directories(tangram INTERFACE ${XMOF2D_INCLUDE_DIRS})

  # XMOF2D_LIBRARIES doesn't contain a real target name. Until we
  # upgrade XMOF2D's cmake to export the correct target, we have to
  # also specify where to find the library
  target_link_directories(tangram INTERFACE ${XMOF2D_LIBRARY_DIR})
  target_link_libraries(tangram INTERFACE ${XMOF2D_LIBRARIES})
  
endif (ENABLE_XMOF2D)


#-----------------------------------------------------------------------------
# Recurse down the source directories building up dependencies
#-----------------------------------------------------------------------------

add_subdirectory(tangram)


# In addition to the include directories of the source, we need to
# include the build or directory to get the autogenerated
# tangram-config.h (The first of these is needed if Wonton is included
# as a submodule, the second is needed for the auto-generated config
# file if Tangram is included as a submodule, the third is to get the
# autogenerated config header if Tangram is being compiled separately
# and the last is for dependencies in installations)

target_include_directories(tangram INTERFACE
  $<BUILD_INTERFACE:${tangram_SOURCE_DIR}>
  $<BUILD_INTERFACE:${tangram_BINARY_DIR}>
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>)


# Tangram targets

install(TARGETS tangram
  EXPORT tangram_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include
  INCLUDES DESTINATION include
  )


#-----------------------------------------------------------------------------
# Add any applications built upon tangram
#-----------------------------------------------------------------------------

add_subdirectory(app)


#-----------------------------------------------------------------------------
# Prepare output for configuration files to be used by projects importing Tangram
#-----------------------------------------------------------------------------

# Write a configuration file from template replacing only variables enclosed
# by the @ sign.
configure_file(${PROJECT_SOURCE_DIR}/cmake/tangramConfig.cmake.in 
  tangramConfig.cmake @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/tangramConfig.cmake DESTINATION lib/cmake/tangram)


# write out a version file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(tangramConfigVersion.cmake
  VERSION "${TANGRAM_MAJOR_VERSION}.${TANGRAM_MINOR_VERSION}.${TANGRAM_PATCH_VERSION}"
  COMPATIBILITY SameMajorVersion)
install(FILES ${PROJECT_BINARY_DIR}/tangramConfigVersion.cmake
  DESTINATION lib/cmake/tangram)


# export targets

install(EXPORT tangram_LIBRARIES
  FILE tangramTargets.cmake
  NAMESPACE tangram::
  EXPORT_LINK_INTERFACE_LIBRARIES
  DESTINATION lib/cmake/tangram)



# Dynamically configured header files that contains defines like
# TANGRAM_ENABLE_MPI etc. if enabled

configure_file(${PROJECT_SOURCE_DIR}/config/tangram-config.h.in
  ${PROJECT_BINARY_DIR}/tangram-config.h @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/tangram-config.h
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include)
