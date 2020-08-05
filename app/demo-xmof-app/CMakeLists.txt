#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

#-------------------------------------------------------------------------------
# demo-xmof-app - uses 2D mesh from Jali
#-------------------------------------------------------------------------------
add_executable(demo-xmof-app demo-xmof-app.cc)
target_link_libraries(demo-xmof-app tangram)

