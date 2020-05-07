#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

#-------------------------------------------------------------------------------
# demo-mof-app - uses 2D/3D mesh from Jali
#-------------------------------------------------------------------------------
add_executable(demo-mof-app demo-mof-app.cc)
target_link_libraries(demo-mof-app tangram)

