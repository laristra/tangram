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

# 3D MOF interface reconstruction on a regular 3x3 grid
# with two linear material interfaces forming a T-junction.
# Uses SimpleMesh.
${RUN_COMMAND} ${TESTAPPDIR}/test_xmof2d 3 3

# Compare the values for the field
${CMPAPPDIR}/apptest_cmp cell_sym_diff_gold0.txt cell_sym_diff_simple_mesh_3x3.txt 1e-14