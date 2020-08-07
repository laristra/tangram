#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

#-----------------------------------------------------------------------------~#
project(tangram_driver)

add_library(tangram_driver INTERFACE)

target_include_directories(tangram_driver INTERFACE
  $<BUILD_INTERFACE:${tangram_SOURCE_DIR}>
  $<BUILD_INTERFACE:${tangram_BINARY_DIR}>
  $<INSTALL_INTERFACE:include>)


set(tangram_driver_HEADERS
  driver.h
  CellMatPoly.h
  write_to_gmv.h
)

# Not yet allowed for INTERFACE libraries
#
# set_target_properties(tangram_driver PROPERTIES
#   PUBLIC_HEADER "${tangram_driver_HEADERS}")
#
# Directly install files instead
#
install(FILES ${tangram_driver_HEADERS} DESTINATION include/tangram/driver)


target_link_libraries(tangram_driver INTERFACE tangram_support)
target_link_libraries(tangram_driver INTERFACE tangram_reconstruct)
target_link_libraries(tangram_driver INTERFACE tangram_intersect)


target_link_libraries(tangram INTERFACE tangram_driver)

install(TARGETS tangram_driver
  EXPORT tangram_LIBRARIES
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION bin
  PUBLIC_HEADER DESTINATION include/tangram/driver
  INCLUDES DESTINATION include/tangram/driver
  )


# Unit tests

if (ENABLE_UNIT_TESTS)
  if (WONTON_ENABLE_Jali)

    tangram_add_unittest(test_cell_matpoly_1d
      SOURCES test/test_cell_matpoly_1d.cc
      LIBRARIES tangram_driver
      POLICY MPI
      THREADS 1
      )
  
    tangram_add_unittest(test_cell_matpoly_2d
      SOURCES test/test_cell_matpoly_2d.cc
      LIBRARIES tangram_driver
      POLICY MPI
      THREADS 1
      )
    
    tangram_add_unittest(test_cell_matpoly_3d
      SOURCES test/test_cell_matpoly_3d.cc
      LIBRARIES tangram_driver
      POLICY MPI
      THREADS 1
      )

    tangram_add_unittest(test_write_gmv_2d
      SOURCES test/test_write_gmv_2d.cc
      LIBRARIES tangram_driver
      POLICY MPI
      THREADS 1
      )
    
    tangram_add_unittest(test_write_gmv_3d
      SOURCES test/test_write_gmv_3d.cc
      LIBRARIES tangram_driver
      POLICY MPI
      THREADS 1
      )
  endif(WONTON_ENABLE_Jali)
endif (ENABLE_UNIT_TESTS)
