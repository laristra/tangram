#[[
This file is part of the Ristra tangram project.
Please see the license file at the root of this repository, or at:
https://github.com/laristra/tangram/blob/master/LICENSE
]]
message(STATUS "Adding test_slic_3d test")

# this app should be run in serial
if(ENABLE_MPI)
  set(RUN_COMMAND "mpirun -np 1")
else()
  set(RUN_COMMAND "")
endif(ENABLE_MPI)


# 3D TEST
macro(ADD_SLIC_3D_TEST TEST_NUM)
  file(COPY slic_3d_test${TEST_NUM}.sh
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY cell_sym_diff_gold${TEST_NUM}_3d.txt
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  add_test(NAME slic_3d_test${TEST_NUM}
    COMMAND ./slic_3d_test${TEST_NUM}.sh)
  set_property(TEST slic_3d_test${TEST_NUM}
    PROPERTY ENVIRONMENT
    TESTAPPDIR=${CMAKE_CURRENT_BINARY_DIR}/..
    CMPAPPDIR=${CMAKE_BINARY_DIR}/app/apptest_cmp
    RUN_COMMAND=${RUN_COMMAND})
  set_property(TEST slic_3d_test${TEST_NUM}
    PROPERTY PROCESSORS 1)
endmacro(ADD_SLIC_3D_TEST)

if (Jali_DIR)
  ADD_SLIC_3D_TEST(2)
  ADD_SLIC_3D_TEST(3)

  file(COPY cubic27.exo DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY voronoi95.exo DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
else(Jali_DIR)
  ADD_SLIC_3D_TEST(0)
  ADD_SLIC_3D_TEST(1)
endif(Jali_DIR)

# 2D TEST
macro(ADD_SLIC_2D_TEST TEST_NUM)
  file(COPY slic_2d_test${TEST_NUM}.sh
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY cell_sym_diff_gold${TEST_NUM}_2d.txt
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  add_test(NAME slic_2d_test${TEST_NUM}
    COMMAND ./slic_2d_test${TEST_NUM}.sh)
  set_property(TEST slic_2d_test${TEST_NUM}
    PROPERTY ENVIRONMENT
    TESTAPPDIR=${CMAKE_CURRENT_BINARY_DIR}/..
    CMPAPPDIR=${CMAKE_BINARY_DIR}/app/apptest_cmp
    RUN_COMMAND=${RUN_COMMAND})
  set_property(TEST slic_2d_test${TEST_NUM}
    PROPERTY PROCESSORS 1)
endmacro(ADD_SLIC_2D_TEST)

if (Jali_DIR)
  ADD_SLIC_2D_TEST(2)
  file(COPY voronoi124.exo DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
else(Jali_DIR)
  ADD_SLIC_2D_TEST(0)
  ADD_SLIC_2D_TEST(1)
endif(Jali_DIR)
