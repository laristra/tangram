#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

#-------------------------------------------------------------------------------
# demo-vof-app - uses 2D/3D mesh from Jali
#-------------------------------------------------------------------------------
add_executable(demo-vof-app demo-vof-app.cc)
target_link_libraries(demo-vof-app tangram)

