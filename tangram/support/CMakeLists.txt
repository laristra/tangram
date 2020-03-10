#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

#-----------------------------------------------------------------------------~#
project(tangram_support)

add_library(tangram_support INTERFACE)

target_include_directories(tangram_support INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIRECTORY}>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIRECTORY}>
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

  cinch_add_unit(test_MatPoly_2D
    SOURCES test/test_MatPoly_2D.cc
    POLICY SERIAL
    LIBRARIES tangram_support
    )

  cinch_add_unit(test_MatPoly_3D
    SOURCES test/test_MatPoly_3D.cc
    POLICY SERIAL
    LIBRARIES tangram_support
    )
endif(ENABLE_UNIT_TESTS)
