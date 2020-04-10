#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

#-----------------------------------------------------------------------------~#


#-----------------------------------------------------------------------------~#

# Add header files
set(support_HEADERS
    tangram.h
    MatPoly.h
    bfgs.h
    PARENT_SCOPE
)

# Add source files
set(support_SOURCES
    PARENT_SCOPE
)

# Unit tests

if (ENABLE_UNIT_TESTS)

  cinch_add_unit(test_matpoly_2d
    SOURCES test/test_matpoly_2d.cc
    POLICY SERIAL
    )

  cinch_add_unit(test_matpoly_3d
    SOURCES test/test_matpoly_3d.cc
    POLICY SERIAL
    )
endif(ENABLE_UNIT_TESTS)
