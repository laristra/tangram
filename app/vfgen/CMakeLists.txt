#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

#-------------------------------------------------------------------------------
# genvf - generate volume fractions and centroids
#-------------------------------------------------------------------------------
add_executable(vfgen-sampling vfgen-sampling.cc)
target_link_libraries(vfgen-sampling ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
  ${MPI_CXX_LIBRARIES} ${EXTRA_LIBS})

#TCMalloc option
set(ENABLE_TCMALLOC TRUE CACHE BOOL "Use TCMalloc")
if(ENABLE_TCMALLOC)
  set(TCMALLOC_LIB "${HOME}" CACHE STRING "Set the TCMalloc library")
  target_link_libraries(vfgen-sampling ${TCMALLOC_LIB})
endif(ENABLE_TCMALLOC)
