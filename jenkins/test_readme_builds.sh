#!/usr/bin/env bash
# This script is executed on Jenkins using
#
#     bash $WORKSPACE/jenkins/test_readme_builds.sh
#
# The exit code determines if the test succeeded or failed.

# Exit on error
set -e
# Echo each command
set -x

# Run build/test commands from README
python $WORKSPACE/jenkins/parseREADME.py $WORKSPACE/README.md $WORKSPACE
