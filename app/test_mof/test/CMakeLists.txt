#[[
This file is part of the Ristra tangram project.
Please see the license file at the root of this repository, or at:
https://github.com/laristra/tangram/blob/master/LICENSE
]]
message(STATUS "Adding test_mof_3d test")

# this app should be run in serial
if(WONTON_ENABLE_MPI)
  set(RUN_COMMAND "mpirun -np 1")
else()
  set(RUN_COMMAND "")
endif(WONTON_ENABLE_MPI)

# 3D TEST
macro(ADD_MOF_3D_TEST TEST_NUM)
  file(COPY mof_3d_test${TEST_NUM}.sh
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY cell_sym_diff_gold${TEST_NUM}_3d.txt
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  add_test(NAME mof_3d_test${TEST_NUM}
    COMMAND ./mof_3d_test${TEST_NUM}.sh)
  set_property(TEST mof_3d_test${TEST_NUM}
    PROPERTY ENVIRONMENT
    TESTAPPDIR=${CMAKE_CURRENT_BINARY_DIR}/..
    CMPAPPDIR=${CMAKE_BINARY_DIR}/app/apptest_cmp
    RUN_COMMAND=${RUN_COMMAND})
  set_property(TEST mof_3d_test${TEST_NUM}
    PROPERTY PROCESSORS 1)
endmacro(ADD_MOF_3D_TEST)

if (WONTON_ENABLE_Jali)
  ADD_MOF_3D_TEST(2)
  ADD_MOF_3D_TEST(3)

  file(COPY cubic27.exo DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY voronoi95.exo DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
else(WONTON_ENABLE_Jali)
  ADD_MOF_3D_TEST(0)
  ADD_MOF_3D_TEST(1)
endif(WONTON_ENABLE_Jali)

# SINGLE-CELL 3D SINGULAR MATRIX TEST
file(COPY test_single_cell_3d.sh
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME test_single_cell_3d
         COMMAND ./test_single_cell_3d.sh)
set_property(TEST test_single_cell_3d
             PROPERTY ENVIRONMENT 
             TESTAPPDIR=${CMAKE_CURRENT_BINARY_DIR}/..
             RUN_COMMAND=${RUN_COMMAND})
set_property(TEST test_single_cell_3d
             PROPERTY PROCESSORS 1) 

# 2D TEST
macro(ADD_MOF_2D_TEST TEST_NUM)
  file(COPY mof_2d_test${TEST_NUM}.sh
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  file(COPY cell_sym_diff_gold${TEST_NUM}_2d.txt
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  add_test(NAME mof_2d_test${TEST_NUM}
    COMMAND ./mof_2d_test${TEST_NUM}.sh)
  set_property(TEST mof_2d_test${TEST_NUM}
    PROPERTY ENVIRONMENT
    TESTAPPDIR=${CMAKE_CURRENT_BINARY_DIR}/..
    CMPAPPDIR=${CMAKE_BINARY_DIR}/app/apptest_cmp
    RUN_COMMAND=${RUN_COMMAND})
  set_property(TEST mof_2d_test${TEST_NUM}
    PROPERTY PROCESSORS 1)
endmacro(ADD_MOF_2D_TEST)

if (WONTON_ENABLE_Jali)
  ADD_MOF_2D_TEST(2)
  file(COPY voronoi124.exo DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
  ADD_MOF_2D_TEST(1_rz)
else(WONTON_ENABLE_Jali)
  ADD_MOF_2D_TEST(0)
  ADD_MOF_2D_TEST(1)
  ADD_MOF_2D_TEST(1_rz)
endif(WONTON_ENABLE_Jali)

