#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

#-------------------------------------------------------------------------------
# tangramapp
#-------------------------------------------------------------------------------
add_executable(tangramapp tangramapp.cc)
target_link_libraries(tangramapp ${EXTRA_LIBS})

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
