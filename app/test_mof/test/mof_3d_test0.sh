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

# 3D MOF interface reconstruction on a regular 3x5x7 grid
# with two planar material interfaces forming a T-junction.
# Mesh cells are NOT decomposed into tetrahedrons.
# Uses SimpleMesh.
${RUN_COMMAND} ${TESTAPPDIR}/test_mof_3d 0 3 5 7

# Compare the values for the field
${CMPAPPDIR}/apptest_cmp cell_sym_diff_gold0_3d.txt cell_sym_diff_simple_mesh_3x5x7.txt 1e-09
