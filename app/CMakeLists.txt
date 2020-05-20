#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

add_subdirectory(multidef)
add_subdirectory(apptest_cmp)
add_subdirectory(tangramapp)
add_subdirectory(test_mof)
add_subdirectory(test_vof)
add_subdirectory(test_lvira)
add_subdirectory(test_lvira+)
add_subdirectory(test_slic)
add_subdirectory(demo_rotor3d)
if (WONTON_ENABLE_MPI AND WONTON_ENABLE_Jali)
  add_subdirectory(demo-vof-app)
  add_subdirectory(demo-mof-app)
  add_subdirectory(demo-slic-app)
endif (WONTON_ENABLE_MPI AND WONTON_ENABLE_Jali)
if(XMOF2D_FOUND)
  add_subdirectory(test_xmof2d)
  if (WONTON_ENABLE_MPI AND WONTON_ENABLE_Jali)
    add_subdirectory(demo-xmof-app)
    add_subdirectory(vfgen)
  endif (WONTON_ENABLE_MPI AND WONTON_ENABLE_Jali)
endif(XMOF2D_FOUND)
