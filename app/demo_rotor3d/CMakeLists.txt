#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

#-------------------------------------------------------------------------------
# Rotor demos for SimpleMesh and Jali
#-------------------------------------------------------------------------------
if(ENABLE_MPI AND Jali_DIR)
  add_executable(demo_rotor_jali demo_rotor_jali.cc)
  target_link_libraries(demo_rotor_jali tangram wonton ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES} 
                        ${MPI_CXX_LIBRARIES} ${EXTRA_LIBS})

  if(MPI_CXX_COMPILE_FLAGS)
    set_target_properties(demo_rotor_jali PROPERTIES
      COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
  endif(MPI_CXX_COMPILE_FLAGS)

  if(MPI_CXX_LINK_FLAGS)
    set_target_properties(demo_rotor_jali PROPERTIES
      LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
  endif(MPI_CXX_LINK_FLAGS)
endif(ENABLE_MPI AND Jali_DIR)

add_executable(demo_rotor demo_rotor.cc)
target_link_libraries(demo_rotor tangram wonton ${MPI_CXX_LIBRARIES} ${EXTRA_LIBS})

if(ENABLE_MPI)
  if(MPI_CXX_COMPILE_FLAGS)
    set_target_properties(demo_rotor PROPERTIES
      COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
  endif(MPI_CXX_COMPILE_FLAGS)

  if(MPI_CXX_LINK_FLAGS)
    set_target_properties(demo_rotor PROPERTIES
      LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
  endif(MPI_CXX_LINK_FLAGS)
endif(ENABLE_MPI)

#TCMalloc option
set(ENABLE_TCMALLOC TRUE CACHE BOOL "Use TCMalloc")
if(ENABLE_TCMALLOC)
  set(TCMALLOC_LIB "${HOME}" CACHE STRING "Set the TCMalloc library")
  target_link_libraries(demo_rotor ${TCMALLOC_LIB})

  if(ENABLE_MPI AND Jali_DIR)
    target_link_libraries(demo_rotor_jali ${TCMALLOC_LIB})
  endif(ENABLE_MPI AND Jali_DIR)
endif(ENABLE_TCMALLOC)

#-------------------------------------------------------------------------------
# Other libraries
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Testing
#-------------------------------------------------------------------------------
#if (ENABLE_APP_TESTS)
#  add_subdirectory(test)
#endif (ENABLE_APP_TESTS)
