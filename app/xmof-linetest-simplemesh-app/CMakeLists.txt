#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

#-------------------------------------------------------------------------------
# xmof-linetest-simplemesh-app - uses 2D SimpleMesh
#-------------------------------------------------------------------------------
add_executable(xmof-linetest-simplemesh-app xmof-linetest-simplemesh-app.cc)
target_link_libraries(xmof-linetest-simplemesh-app ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}                      ${XMOF2D_LIBRARIES} ${MPI_CXX_LIBRARIES} ${EXTRA_LIBS})

if(ENABLE_MPI)
  if(MPI_CXX_COMPILE_FLAGS)
    set_target_properties(xmof-linetest-simplemesh-app PROPERTIES
      COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
  endif(MPI_CXX_COMPILE_FLAGS)

  if(MPI_CXX_LINK_FLAGS)
    set_target_properties(xmof-linetest-simplemesh-app PROPERTIES
      LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
  endif(MPI_CXX_LINK_FLAGS)
endif(ENABLE_MPI)

#TCMalloc option
set(ENABLE_TCMALLOC TRUE CACHE BOOL "Use TCMalloc")
if(ENABLE_TCMALLOC)
  set(TCMALLOC_LIB "${HOME}" CACHE STRING "Set the TCMalloc library")
  target_link_libraries(xmof-linetest-simplemesh-app ${TCMALLOC_LIB})
endif(ENABLE_TCMALLOC)                      
