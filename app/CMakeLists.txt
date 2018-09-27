#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

add_subdirectory(apptest_cmp)
add_subdirectory(tangramapp)
add_subdirectory(test_mof)
add_subdirectory(test_vof)
add_subdirectory(demo-slic-app3d)
add_subdirectory(test_mof_3d)
if (ENABLE_MPI AND Jali_DIR)
  add_subdirectory(demo-vof-app)
  add_subdirectory(demo-mof-app)
endif (ENABLE_MPI AND Jali_DIR)
if(XMOF2D_FOUND)
  add_subdirectory(test_xmof2d)
  if (ENABLE_MPI AND Jali_DIR)
    add_subdirectory(demo-xmof-app)
    add_subdirectory(vfgen)
  endif (ENABLE_MPI AND Jali_DIR)
endif(XMOF2D_FOUND)
