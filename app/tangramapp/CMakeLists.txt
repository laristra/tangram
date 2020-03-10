#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

#-------------------------------------------------------------------------------
# tangramapp
#-------------------------------------------------------------------------------
add_executable(tangramapp tangramapp.cc)
target_link_libraries(tangramapp tangram)

