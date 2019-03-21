#!/bin/bash
: <<'END'
This file is part of the Ristra tangram project.
Please see the license file at the root of this repository, or at:
https://github.com/laristra/tangram/blob/master/LICENSE
END


# Exit on error
set -e
# Echo each command
set -x

# 2D MOF interface reconstruction on a voronoi mesh
# with two planar material interfaces forming a T-junction.
# Mesh cells ARE decomposed into triangles.
# Uses Jali. 
${RUN_COMMAND} ${TESTAPPDIR}/test_mof_2d 1 voronoi124.exo

# Compare the values for the field
${CMPAPPDIR}/apptest_cmp cell_sym_diff_gold2_2d.txt cell_sym_diff_2d_voronoi124_decomposed.txt 1e-11
