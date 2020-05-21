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

# 2D LVIRA interface reconstruction on a regular 3x5 grid
# with two planar material interfaces forming a T-junction.
# Mesh cells ARE decomposed into triangles..
# Uses SimpleMesh.
${RUN_COMMAND} ${TESTAPPDIR}/test_lvira_2d 1 3 5 

# Compare the values for the field
${CMPAPPDIR}/apptest_cmp cell_sym_diff_gold1_2d.txt cell_sym_diff_2d_simple_mesh_3x5_decomposed.txt 1e-9
