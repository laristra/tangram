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

# The tests that are expected to fail are prefixed by the "!" Bash operator
# that flips the exit code. The command must be surrounded in a sub-shell with
# the (), otherwise it does not work with set -e. See
# http://stackoverflow.com/a/31549913/479532 for details.

${TESTAPPDIR}/apptest_cmp field1.txt field1.txt 1e-12
(! ${TESTAPPDIR}/apptest_cmp field1.txt field2.txt 1e-12)
${TESTAPPDIR}/apptest_cmp field1.txt field2.txt 2
(! ${TESTAPPDIR}/apptest_cmp field1.txt field3.txt 1e-12)
(! ${TESTAPPDIR}/apptest_cmp field1.txt field4.txt 1e-12)
