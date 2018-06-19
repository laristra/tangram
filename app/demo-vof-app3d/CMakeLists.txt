#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]


#-------------------------------------------------------------------------------
# demo-vof-app3d - uses 3D Simple_Mesh
#-------------------------------------------------------------------------------
add_executable(demo-vof-app3d demo-vof-app3d.cc)
target_link_libraries(demo-vof-app3d tangram ${EXTRA_LIBS})

#-------------------------------------------------------------------------------
# Other libraries
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Testing
#-------------------------------------------------------------------------------
# if (ENABLE_APP_TESTS)
#   add_subdirectory(test)
# endif (ENABLE_APP_TESTS)
