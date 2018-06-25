# Add header files
set(test_vfgen_2D_fill_HEADERS
  PARENT_SCOPE
)

set(test_vfgen_2D_fill_SOURCES
  test_vfgen_2D_fill.cc
  PARENT_SCOPE
)

cinch_add_unit(test_vfgen_2D_fill
  SOURCES test_vfgen_2D_fill.cc 
  LIBRARIES ${Jali_LIBRARIES} ${Jali_TPL_LIBRARIES}
  POLICY MPI
  THREADS 1
  )
  
  
