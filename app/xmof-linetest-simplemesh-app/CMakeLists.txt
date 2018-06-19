#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

#-------------------------------------------------------------------------------
# xmof-linetest-simplemesh-app - uses 2D mesh SimpleMesh
#-------------------------------------------------------------------------------
add_executable(xmof-linetest-simplemesh-app xmof-linetest-simplemesh-app.cc)
target_link_libraries(xmof-linetest-simplemesh-app 
                      ${XMOF2D_LIBRARIES} ${EXTRA_LIBS})
