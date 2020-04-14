#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

#-----------------------------------------------------------------------------~#
project(tangram_reconstruct)


add_library(tangram_reconstruct INTERFACE)

target_include_directories(tangram_reconstruct INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIRECTORY}>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIRECTORY}>
  $<INSTALL_INTERFACE:include>)


if (XMOF2D_DIR)
  set(xmof2D_wrapper_HEADERS
    xmof2D_wrapper.h
    )
endif (XMOF2D_DIR)

set(tangram_reconstruct_HEADERS
  SLIC.h
  cutting_distance_solver.h
  nested_dissections.h
  VOF.h
  MOF.h
  )

# Not yet allowed for INTERFACE libraries
#
# set_target_properties(tangram_reconstruct PROPERTIES
#   PUBLIC_HEADER "${xmof2D_wrapper_HEADERS};${tangram_reconstruct_HEADERS}")
#
# Directly install files instead
#
install(FILES ${xmof2D_wrapper_HEADERS} ${tangram_reconstruct_HEADERS}
  DESTINATION include/tangram/reconstruct)


target_link_libraries(tangram_reconstruct INTERFACE tangram_support)


target_link_libraries(tangram INTERFACE tangram_reconstruct)

install(TARGETS tangram_reconstruct
  EXPORT tangram_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/tangram/reconstruct
  INCLUDES DESTINATION include/tangram/reconstruct
  )


# Unit tests

if (ENABLE_UNIT_TESTS)
  cinch_add_unit(cutting_distance_solver
    SOURCES test/test_cutting_distance_solver.cc
    LIBRARIES tangram_reconstruct
    POLICY SERIAL
    )
endif (ENABLE_UNIT_TESTS)
