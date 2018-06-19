#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

#-------------------------------------------------------------------------------
# demoapp3 - uses Simple_Mesh from portage
#-------------------------------------------------------------------------------
add_executable(demoapp3d demoapp3d.cc)
target_link_libraries(demoapp3d tangram ${EXTRA_LIBS})

#-------------------------------------------------------------------------------
# Other libraries
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Testing
#-------------------------------------------------------------------------------
# if (ENABLE_APP_TESTS)
#   add_subdirectory(test)
# endif (ENABLE_APP_TESTS)
