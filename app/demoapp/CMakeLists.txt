#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

#-------------------------------------------------------------------------------
# demoapp - uses 3D Simple_Mesh
#-------------------------------------------------------------------------------
add_executable(demoapp demoapp.cc)
target_link_libraries(demoapp tangram ${EXTRA_LIBS})


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
