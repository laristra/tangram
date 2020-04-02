#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

#-----------------------------------------------------------------------------~#
project(tangram_driver)

add_library(tangram_driver INTERFACE)

target_include_directories(tangram_driver INTERFACE
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIRECTORY}>
  $<BUILD_INTERFACE:${PROJECT_BINARY_DIRECTORY}>
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

    cinch_add_unit(test_CellMatPoly_1D
      SOURCES test/test_CellMatPoly_1D.cc
      LIBRARIES tangram_driver
      POLICY MPI
      THREADS 1
      )
  
    cinch_add_unit(test_CellMatPoly_2D
      SOURCES test/test_CellMatPoly_2D.cc
      LIBRARIES tangram_driver
      POLICY MPI
      THREADS 1
      )
    
    cinch_add_unit(test_CellMatPoly_3D
      SOURCES test/test_CellMatPoly_3D.cc
      LIBRARIES tangram_driver
      POLICY MPI
      THREADS 1
      )

    cinch_add_unit(test_write_gmv2D
      SOURCES test/test_write_gmv2D.cc
      LIBRARIES tangram_driver
      POLICY MPI
      THREADS 1
      )
    
    cinch_add_unit(test_write_gmv3D
      SOURCES test/test_write_gmv3D.cc
      LIBRARIES tangram_driver
      POLICY MPI
      THREADS 1
      )
  endif(WONTON_ENABLE_Jali)
endif (ENABLE_UNIT_TESTS)
