#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

#-----------------------------------------------------------------------------~#


#-----------------------------------------------------------------------------~#

# Add header files
set(driver_HEADERS
  driver.h
  CellMatPoly.h
  write_to_gmv.h
  PARENT_SCOPE
)

# Add source files
set(driver_SOURCES
    PARENT_SCOPE
)

# Unit tests

if (ENABLE_UNIT_TESTS)
  if (Jali_DIR)

    cinch_add_unit(test_cell_matpoly_1d
      SOURCES test/test_cell_matpoly_1d.cc
      LIBRARIES ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
      POLICY MPI
      THREADS 1
      )
  
    cinch_add_unit(test_cell_matpoly_2d
      SOURCES test/test_cell_matpoly_2d.cc
      LIBRARIES ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
      POLICY MPI
      THREADS 1
      )
    
    cinch_add_unit(test_cell_matpoly_3d
      SOURCES test/test_cell_matpoly_3d.cc
      LIBRARIES ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
      POLICY MPI
      THREADS 1
      )

    cinch_add_unit(test_write_gmv_2d
      SOURCES test/test_write_gmv_2d.cc
      LIBRARIES ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
      POLICY MPI
      THREADS 1
      )
    
    cinch_add_unit(test_write_gmv_3d
      SOURCES test/test_write_gmv_3d.cc
      LIBRARIES ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
      POLICY MPI
      THREADS 1
      )
  endif(Jali_DIR)
endif (ENABLE_UNIT_TESTS)
