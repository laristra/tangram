#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

#-------------------------------------------------------------------------------
# demo-slic-app3d - uses 3D Simple_Mesh
#-------------------------------------------------------------------------------
add_executable(demo-slic-app demo-slic-app.cc)
target_link_libraries(demo-slic-app tangram)
