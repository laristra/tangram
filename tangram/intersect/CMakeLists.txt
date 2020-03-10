#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

#-----------------------------------------------------------------------------~#
project(tangram_intersect)


add_library(tangram_intersect INTERFACE)

target_include_directories(tangram_intersect INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIRECTORY}>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIRECTORY}>
  $<INSTALL_INTERFACE:include>)

set(tangram_intersect_HEADERS
  split_r2d.h
  split_r3d.h)

# Not yet allowed for INTERFACE libraries
#
# set_target_properties(tangram_intersect PROPERTIES
#   PUBLIC_HEADER "${tangram_intersect_HEADERS}")
#
# Directly install files instead
#
install(FILES ${tangram_intersect_HEADERS}
  DESTINATION include/tangram/intersect)


target_link_libraries(tangram_intersect INTERFACE tangram_support)


target_link_libraries(tangram INTERFACE tangram_intersect)

install(TARGETS tangram_intersect
  EXPORT tangram_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/tangram/intersect
  INCLUDES DESTINATION include/tangram/intersect
  )

# Unit tests

if (ENABLE_UNIT_TESTS) 

  cinch_add_unit(split_r3d
    SOURCES test/test_split_r3d.cc
    LIBRARIES tangram_intersect
    POLICY SERIAL
    )
  cinch_add_unit(split_r2d
    SOURCES test/test_split_r2d.cc
    LIBRARIES tangram_intersect
    POLICY SERIAL
    )
  cinch_add_unit(split_order
    SOURCES test/test_split_order.cc
    LIBRARIES tangram_intersect
    POLICY SERIAL
    )
endif(ENABLE_UNIT_TESTS)

