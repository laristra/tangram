#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

#-------------------------------------------------------------------------------
# genvf - generate volume fractions and centroids
#-------------------------------------------------------------------------------
  if (WONTON_ENABLE_MPI AND WONTON_ENABLE_Jali)
    add_executable(vfgen-sampling-jali vfgen-sampling-jali.cc)
    target_link_libraries(vfgen-sampling-jali tangram)
  endif (WONTON_ENABLE_MPI AND WONTON_ENABLE_Jali)

  add_executable(vfgen-sampling-sm vfgen-sampling-sm.cc)
  target_link_libraries(vfgen-sampling-sm tangram)

