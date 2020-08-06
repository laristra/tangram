#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

#-------------------------------------------------------------------------------
# test_mof_3d - uses 3D mesh from SimpleMesh or Jali
#-------------------------------------------------------------------------------
add_executable(test_mof_3d test_mof_3d.cc)
target_link_libraries(test_mof_3d wonton ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES} 
                                  ${MPI_CXX_LIBRARIES} ${EXTRA_LIBS})

if(ENABLE_MPI)
  if(MPI_CXX_COMPILE_FLAGS)
    set_target_properties(test_mof_3d PROPERTIES
      COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
  endif(MPI_CXX_COMPILE_FLAGS)

  if(MPI_CXX_LINK_FLAGS)
    set_target_properties(test_mof_3d PROPERTIES
      LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
  endif(MPI_CXX_LINK_FLAGS)
endif(ENABLE_MPI)

#TCMalloc option
set(ENABLE_TCMALLOC TRUE CACHE BOOL "Use TCMalloc")
if(ENABLE_TCMALLOC)
  set(TCMALLOC_LIB "${HOME}" CACHE STRING "Set the TCMalloc library")
  target_link_libraries(test_mof_3d ${TCMALLOC_LIB})
endif(ENABLE_TCMALLOC)

#-------------------------------------------------------------------------------
# test_mof_2d - uses 2D mesh from SimpleMesh or Jali
#-------------------------------------------------------------------------------
add_executable(test_mof_2d test_mof_2d.cc)
target_link_libraries(test_mof_2d wonton ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES} 
                                  ${MPI_CXX_LIBRARIES} ${EXTRA_LIBS})

if(ENABLE_MPI)
  if(MPI_CXX_COMPILE_FLAGS)
    set_target_properties(test_mof_2d PROPERTIES
      COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
  endif(MPI_CXX_COMPILE_FLAGS)

  if(MPI_CXX_LINK_FLAGS)
    set_target_properties(test_mof_2d PROPERTIES
      LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
  endif(MPI_CXX_LINK_FLAGS)
endif(ENABLE_MPI)

#TCMalloc option
set(ENABLE_TCMALLOC TRUE CACHE BOOL "Use TCMalloc")
if(ENABLE_TCMALLOC)
  set(TCMALLOC_LIB "${HOME}" CACHE STRING "Set the TCMalloc library")
  target_link_libraries(test_mof_2d ${TCMALLOC_LIB})
endif(ENABLE_TCMALLOC)
#-------------------------------------------------------------------------------
# Other libraries
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Testing
#-------------------------------------------------------------------------------
if (ENABLE_APP_TESTS)
  add_subdirectory(test)
endif (ENABLE_APP_TESTS)
