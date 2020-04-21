#[[
This file is part of the Ristra tangram project.
Please see the license file at the root of this repository, or at:
https://github.com/laristra/tangram/blob/master/LICENSE
]]
message(STATUS "Adding test_xmof2d test")

# this app should be run in serial
if(WONTON_ENABLE_MPI)
  set(RUN_COMMAND "mpirun -np 1")
else()
  set(RUN_COMMAND "")
endif(WONTON_ENABLE_MPI)

macro(ADD_XMOF2D_TEST TEST_NUM)
  file(COPY xmof2d_test${TEST_NUM}.sh
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY cell_sym_diff_gold${TEST_NUM}.txt
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  add_test(NAME xmof2d_test${TEST_NUM}
    COMMAND ./xmof2d_test${TEST_NUM}.sh)
  set_property(TEST xmof2d_test${TEST_NUM}
    PROPERTY ENVIRONMENT
    TESTAPPDIR=${CMAKE_CURRENT_BINARY_DIR}/..
    CMPAPPDIR=${CMAKE_BINARY_DIR}/app/apptest_cmp
    RUN_COMMAND=${RUN_COMMAND})
  set_property(TEST xmof2d_test${TEST_NUM}
    PROPERTY PROCESSORS 1)
endmacro(ADD_XMOF2D_TEST)

if (WONTON_ENABLE_Jali)
  ADD_XMOF2D_TEST(1)

  file(COPY voronoi124.exo DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
else(WONTON_ENABLE_Jali)
  ADD_XMOF2D_TEST(0)
endif(WONTON_ENABLE_Jali)