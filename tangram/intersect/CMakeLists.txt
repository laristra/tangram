#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

#-----------------------------------------------------------------------------~#


#-----------------------------------------------------------------------------~#

# Add header files
set(intersect_HEADERS
    r2d.h
    r3d.h
    split_r2d.h
    split_r3d.h
    PARENT_SCOPE
)

# Add source files
set(intersect_SOURCES
    r2d.c
    r3d.c
    PARENT_SCOPE
)

# Unit tests

if (ENABLE_UNIT_TESTS) 

  cinch_add_unit(split_r3d
    SOURCES test/test_split_r3d.cc r3d.c
    LIBRARIES tangram
    POLICY SERIAL
    )
  cinch_add_unit(split_r2d
    SOURCES test/test_split_r2d.cc r2d.c
    LIBRARIES tangram
    POLICY SERIAL
    )
  cinch_add_unit(split_order
    SOURCES test/test_split_order.cc r2d.c r3d.c
    LIBRARIES tangram
    POLICY SERIAL
    )
endif(ENABLE_UNIT_TESTS)

