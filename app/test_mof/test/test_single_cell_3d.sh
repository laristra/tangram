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

# 3D MOF interface reconstruction on a single two-material cell.
# Reproduces a test case that at some point resulted
# in a numerically singular matrix B being inverted in BFGS 
# Uses SimpleMesh.
${RUN_COMMAND} ${TESTAPPDIR}/test_single_cell_3d
