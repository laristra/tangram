#[[
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
]]

#-----------------------------------------------------------------------------~#


#-----------------------------------------------------------------------------~#

if (XMOF2D_DIR)
    set(xmof2D_wrapper_headers
        reconstruct/xmof2D_wrapper.h
        )
endif (XMOF2D_DIR)

# Add header files
set(reconstruct_HEADERS
    SLIC.h
    xmof2D_wrapper.h
    cutting_distance_solver.h
    nested_dissections.h
    VOF.h
    PARENT_SCOPE
)

# Add source files
set(reconstruct_SOURCES
    PARENT_SCOPE
)

# Unit tests

if (ENABLE_UNIT_TESTS)
  cinch_add_unit(cutting_distance_solver
    SOURCES test/test_cutting_distance_solver.cc
    LIBRARIES tangram
    POLICY SERIAL
    )
endif (ENABLE_UNIT_TESTS)
