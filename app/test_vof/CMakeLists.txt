#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

#-------------------------------------------------------------------------------
# test_vof_3d - uses 3D mesh from SimpleMesh or Jali
#-------------------------------------------------------------------------------
add_executable(test_vof_3d test_vof_3d.cc)
target_link_libraries(test_vof_3d tangram)


#-------------------------------------------------------------------------------
# test_vof_2d - uses 2D mesh from SimpleMesh or Jali
#-------------------------------------------------------------------------------
add_executable(test_vof_2d test_vof_2d.cc)
target_link_libraries(test_vof_2d tangram)

#-------------------------------------------------------------------------------
# Testing
#-------------------------------------------------------------------------------
if (ENABLE_APP_TESTS)
  add_subdirectory(test)
endif (ENABLE_APP_TESTS)
