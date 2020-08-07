#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

#-------------------------------------------------------------------------------
# Rotor demo for SimpleMesh
#-------------------------------------------------------------------------------
add_executable(demo_rotor demo_rotor.cc)
target_link_libraries(demo_rotor tangram)

