#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

#-------------------------------------------------------------------------------
# demo-xmof-app - uses 2D mesh from Jali
#-------------------------------------------------------------------------------
add_executable(demo-xmof-app demo-xmof-app.cc)
target_link_libraries(demo-xmof-app ${XMOF2D_LIBRARIES} ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES} 
                                    ${MPI_CXX_LIBRARIES} ${EXTRA_LIBS})

#TCMalloc option
set(ENABLE_TCMALLOC TRUE CACHE BOOL "Use TCMalloc")
if(ENABLE_TCMALLOC)
  set(TCMALLOC_LIB "${HOME}" CACHE STRING "Set the TCMalloc library")
  target_link_libraries(demo-xmof-app ${TCMALLOC_LIB})
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
