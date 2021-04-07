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

# 2D MOF interface reconstruction on a regular 10x10 RZ grid
# with two planar material interfaces forming a tilted T-junction.
# Uses SimpleMesh.
${RUN_COMMAND} ${TESTAPPDIR}/test_mof_2d_rz 1 10 10 1

# Compare the values for the field
${CMPAPPDIR}/apptest_cmp cell_sym_diff_gold1_rz_2d.txt cell_sym_diff_2d_simple_mesh_decomposed.txt 1e-11
