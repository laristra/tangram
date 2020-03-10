#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

#-------------------------------------------------------------------------------
# genvf - generate volume fractions and centroids
#-------------------------------------------------------------------------------
add_executable(vfgen-sampling vfgen-sampling.cc)
target_link_libraries(vfgen-sampling tangram)
