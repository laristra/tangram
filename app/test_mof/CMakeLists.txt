#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

#-------------------------------------------------------------------------------
# test_mof_3d - uses 3D mesh from SimpleMesh or Jali
#-------------------------------------------------------------------------------
add_executable(test_mof_3d test_mof_3d.cc)
target_link_libraries(test_mof_3d tangram)

#-------------------------------------------------------------------------------
# test_single_cell_3d - uses a single-cell 3D mesh from SimpleMesh
#-------------------------------------------------------------------------------
add_executable(test_single_cell_3d test_single_cell_3d.cc)
target_link_libraries(test_single_cell_3d tangram)

#-------------------------------------------------------------------------------
# test_mof_2d - uses 2D mesh from SimpleMesh or Jali
#-------------------------------------------------------------------------------
add_executable(test_mof_2d test_mof_2d.cc)
target_link_libraries(test_mof_2d tangram)

#-------------------------------------------------------------------------------
# Testing
#-------------------------------------------------------------------------------
if (ENABLE_APP_TESTS)
  add_subdirectory(test)
endif (ENABLE_APP_TESTS)
