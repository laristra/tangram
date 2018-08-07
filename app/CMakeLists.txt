#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

add_subdirectory(tangramapp)
add_subdirectory(demoapp)
add_subdirectory(demoapp3d)
if (ENABLE_MPI AND Jali_DIR)
  add_subdirectory(demo-vof-app)
  add_subdirectory(demo-mof-app3d)
endif (ENABLE_MPI AND Jali_DIR)
if(XMOF2D_FOUND)
  add_subdirectory(xmof-linetest-simplemesh-app)
  if (ENABLE_MPI AND Jali_DIR)
    add_subdirectory(demo-xmof-app)
    add_subdirectory(vfgen)
  endif (ENABLE_MPI AND Jali_DIR)
endif(XMOF2D_FOUND)
