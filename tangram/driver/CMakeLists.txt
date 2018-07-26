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

    cinch_add_unit(test_CellMatPoly_1D
      SOURCES test/test_CellMatPoly_1D.cc
      LIBRARIES ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
      POLICY MPI
      THREADS 1
      )
  
    cinch_add_unit(test_CellMatPoly_2D
      SOURCES test/test_CellMatPoly_2D.cc
      LIBRARIES ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
      POLICY MPI
      THREADS 1
      )
    
    cinch_add_unit(test_CellMatPoly_3D
      SOURCES test/test_CellMatPoly_3D.cc
      LIBRARIES ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
      POLICY MPI
      THREADS 1
      )

    cinch_add_unit(test_write_gmv2D
      SOURCES test/test_write_gmv2D.cc
      LIBRARIES ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
      POLICY MPI
      THREADS 1
      )
    
    cinch_add_unit(test_write_gmv3D
      SOURCES test/test_write_gmv3D.cc
      LIBRARIES ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
      POLICY MPI
      THREADS 1
      )
  endif(Jali_DIR)
endif (ENABLE_UNIT_TESTS)
