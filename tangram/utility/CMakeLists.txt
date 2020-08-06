#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

#-----------------------------------------------------------------------------~#
project(tangram_utility)

add_library(tangram_utility INTERFACE)

target_include_directories(tangram_utility INTERFACE
  $<BUILD_INTERFACE:${tangram_SOURCE_DIR}>
  $<BUILD_INTERFACE:${tangram_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>)

set(tangram_utility_HEADERS
  read_material_data.h
  get_material_moments.h
  get_mat_sym_diff_vol.h
  advance_moments.h
  rk4.h
  )
set(rpgtools_HEADERS
  rpgtools/primitives.h
  rpgtools/cuts.h
  )
set(rpgtools_examples_HEADERS
  rpgtools/examples/matdata_rotor3d.h
  )

# Not yet allowed for INTERFACE libraries
#
# set_target_properties(tangram_utility PROPERTIES
#  PUBLIC_HEADER
#  "${tangram_utility_HEADERS} ${rpgtools_HEADERS} ${rpgtools_examples_HEADERS}"
#  )
#
# Instead directly install files
#
install(FILES ${tangram_utility_HEADERS} DESTINATION include/tangram/utility)
install(FILES ${rpgtools_HEADERS} DESTINATION include/tangram/utility/rpgtools)
install(FILES ${rpgtools_examples_HEADERS}
  DESTINATION include/tangram/utility/rpgtools/examples)


target_link_libraries(tangram_utility INTERFACE tangram_support)

target_link_libraries(tangram INTERFACE tangram_utility)

install(TARGETS tangram_utility
  EXPORT tangram_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/tangram/utility
  INCLUDES DESTINATION include/tangram/utility
  )

