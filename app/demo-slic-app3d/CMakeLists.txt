#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

#-------------------------------------------------------------------------------
# demo-slic-app3d - uses 3D Simple_Mesh
#-------------------------------------------------------------------------------
add_executable(demo-slic-app3d demo-slic-app3d.cc)
target_link_libraries(demo-slic-app3d wonton ${MPI_CXX_LIBRARIES} ${EXTRA_LIBS})

if(ENABLE_MPI)
  if(MPI_CXX_COMPILE_FLAGS)
    set_target_properties(demo-slic-app3d PROPERTIES
      COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
  endif(MPI_CXX_COMPILE_FLAGS)

  if(MPI_CXX_LINK_FLAGS)
    set_target_properties(demo-slic-app3d PROPERTIES
      LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
  endif(MPI_CXX_LINK_FLAGS)
endif(ENABLE_MPI)

#TCMalloc option
set(ENABLE_TCMALLOC TRUE CACHE BOOL "Use TCMalloc")
if(ENABLE_TCMALLOC)
  set(TCMALLOC_LIB "${HOME}" CACHE STRING "Set the TCMalloc library")
  target_link_libraries(demo-slic-app3d ${TCMALLOC_LIB})
endif(ENABLE_TCMALLOC)

#-------------------------------------------------------------------------------
# Other libraries
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Testing
#-------------------------------------------------------------------------------
# if (ENABLE_APP_TESTS)
#   add_subdirectory(test)
# endif (ENABLE_APP_TESTS)
