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
    rk4.h
    PARENT_SCOPE
)

# Add source files
set(support_SOURCES
    PARENT_SCOPE
)

# Unit tests

if (ENABLE_UNIT_TESTS)

  cinch_add_unit(test_MatPoly_2D
    SOURCES test/test_MatPoly_2D.cc
    POLICY SERIAL
    )

  cinch_add_unit(test_MatPoly_3D
    SOURCES test/test_MatPoly_3D.cc
    POLICY SERIAL
    )
endif(ENABLE_UNIT_TESTS)
