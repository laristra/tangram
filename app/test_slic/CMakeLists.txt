#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

#-------------------------------------------------------------------------------
# test_slic_2d - uses 2D mesh from SimpleMesh or Jali
#-------------------------------------------------------------------------------
add_executable(test_slic_2d test_slic_2d.cc)
target_link_libraries(test_slic_2d tangram ${Jali_TPL_LIBRARIES})

if(Jali_DIR)
  target_link_libraries(test_slic_2d tangram wonton ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES})
endif(Jali_DIR)

if(ENABLE_MPI)
  target_link_libraries(test_slic_2d tangram ${MPI_CXX_LIBRARIES} ${EXTRA_LIBS})
  if(MPI_CXX_COMPILE_FLAGS)
    set_target_properties(test_slic_2d PROPERTIES
      COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
  endif(MPI_CXX_COMPILE_FLAGS)

  if(MPI_CXX_LINK_FLAGS)
    set_target_properties(test_slic_2d PROPERTIES
      LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
  endif(MPI_CXX_LINK_FLAGS)
endif(ENABLE_MPI)

#TCMalloc option
set(ENABLE_TCMALLOC TRUE CACHE BOOL "Use TCMalloc")
if(ENABLE_TCMALLOC)
  set(TCMALLOC_LIB "${HOME}" CACHE STRING "Set the TCMalloc library")
  target_link_libraries(test_slic_2d ${TCMALLOC_LIB})
endif(ENABLE_TCMALLOC)

#-------------------------------------------------------------------------------
# test_slic_3d - uses 3D mesh from SimpleMesh or Jali
#-------------------------------------------------------------------------------
add_executable(test_slic_3d test_slic_3d.cc)
target_link_libraries(test_slic_3d tangram wonton ${Jali_TPL_LIBRARIES})

if(Jali_DIR)
  target_link_libraries(test_slic_3d tangram ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES})
endif(Jali_DIR)

if(ENABLE_MPI)
  target_link_libraries(test_slic_3d tangram ${MPI_CXX_LIBRARIES} ${EXTRA_LIBS})
  if(MPI_CXX_COMPILE_FLAGS)
    set_target_properties(test_slic_3d PROPERTIES
      COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
  endif(MPI_CXX_COMPILE_FLAGS)

  if(MPI_CXX_LINK_FLAGS)
    set_target_properties(test_slic_3d PROPERTIES
      LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
  endif(MPI_CXX_LINK_FLAGS)
endif(ENABLE_MPI)

#TCMalloc option
if(ENABLE_TCMALLOC)
  target_link_libraries(test_slic_3d ${TCMALLOC_LIB})
endif(ENABLE_TCMALLOC)

#-------------------------------------------------------------------------------
# Testing
#-------------------------------------------------------------------------------
if (ENABLE_APP_TESTS)
  add_subdirectory(test)
endif (ENABLE_APP_TESTS)
