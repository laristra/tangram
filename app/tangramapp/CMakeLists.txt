#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

#-------------------------------------------------------------------------------
# tangramapp
#-------------------------------------------------------------------------------
add_executable(tangramapp tangramapp.cc)
target_link_libraries(tangramapp ${MPI_CXX_LIBRARIES} ${EXTRA_LIBS})

if(ENABLE_MPI)
  if(MPI_CXX_COMPILE_FLAGS)
    set_target_properties(tangramapp PROPERTIES
      COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
  endif(MPI_CXX_COMPILE_FLAGS)

  if(MPI_CXX_LINK_FLAGS)
    set_target_properties(tangramapp PROPERTIES
      LINK_FLAGS "${MPI_CXX_LINK_FLAGS}")
  endif(MPI_CXX_LINK_FLAGS)
endif(ENABLE_MPI)

#TCMalloc option
set(ENABLE_TCMALLOC TRUE CACHE BOOL "Use TCMalloc")
if(ENABLE_TCMALLOC)
  set(TCMALLOC_LIB "${HOME}" CACHE STRING "Set the TCMalloc library")
  target_link_libraries(tangramapp ${TCMALLOC_LIB})
endif(ENABLE_TCMALLOC)

#-------------------------------------------------------------------------------
# Other libraries
#-------------------------------------------------------------------------------
# taget_link_libraries(tangramapp ${WONTON_LIBS})
# include_directories(${WONTON_INCLUDE_DIRS})

#-------------------------------------------------------------------------------
# Testing
#-------------------------------------------------------------------------------
# if (ENABLE_APP_TESTS)
#   add_subdirectory(test)
# endif (ENABLE_APP_TESTS)
