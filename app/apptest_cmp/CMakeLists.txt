#[[
This file is part of the Ristra tangram project.
Please see the license file at the root of this repository, or at:
https://github.com/laristra/tangram/blob/master/LICENSE
]]


# Comparator for testing infrastructure
add_executable(apptest_cmp apptest_cmp.cc)

if (ENABLE_APP_TESTS)
  add_subdirectory(test)
endif (ENABLE_APP_TESTS)
