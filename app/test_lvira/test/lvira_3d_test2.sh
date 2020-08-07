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

# 3D LVIRA interface reconstruction on a regular 3x3x3 grid
# with two planar material interfaces forming a T-junction.
# Mesh cells are NOT decomposed into tetrahedrons.
# Uses Jali.
${RUN_COMMAND} ${TESTAPPDIR}/test_lvira_3d 0 cubic27.exo

# Compare the values for the field
${CMPAPPDIR}/apptest_cmp cell_sym_diff_gold2_3d.txt cell_sym_diff_cubic27.txt 1e-06
