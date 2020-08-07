#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

#-----------------------------------------------------------------------------~#
project(tangram_support)

add_library(tangram_support INTERFACE)

target_include_directories(tangram_support INTERFACE
  $<BUILD_INTERFACE:${tangram_SOURCE_DIR}>
  $<BUILD_INTERFACE:${tangram_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>)

set(tangram_support_HEADERS
  tangram.h
  MatPoly.h
  bfgs.h
  )

# Not yet allowed for INTERFACE libraries
# 
# set_target_properties(tangram_support PROPERTIES
#   PUBLIC_HEADER "${tangram_support_HEADERS}")
#
# Directly install files instead
#
install(FILES ${tangram_support_HEADERS} DESTINATION include/tangram/support)

target_link_libraries(tangram_support INTERFACE ${WONTON_LIBRARIES})

if (WONTON_ENABLE_MPI)
  target_compile_definitions(tangram_support INTERFACE WONTON_ENABLE_MPI)
endif ()

if (WONTON_ENABLE_THRUST)
  target_compile_definitions(tangram_support INTERFACE WONTON_ENABLE_THRUST)
endif ()

target_link_libraries(tangram INTERFACE tangram_support)

install(TARGETS tangram_support
  EXPORT tangram_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/tangram/support
  INCLUDES DESTINATION include/tangram/support
  )

# Unit tests

if (ENABLE_UNIT_TESTS)

  tangram_add_unittest(test_matpoly_2d
    SOURCES test/test_matpoly_2d.cc
    POLICY SERIAL
    LIBRARIES tangram_support
    )

  tangram_add_unittest(test_matpoly_3d
    SOURCES test/test_matpoly_3d.cc
    POLICY SERIAL
    LIBRARIES tangram_support
    )
endif(ENABLE_UNIT_TESTS)
