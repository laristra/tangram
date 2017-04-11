#!/usr/bin/env bash
# This script is executed on Jenkins using
#
#     $WORKSPACE/jenkins/build_matrix_entry.sh <compiler> <build_type>
#
# The exit code determines if the test succeeded or failed.
# Note that the environment variable WORKSPACE must be set (Jenkins
# will do this automatically).

# Exit on error
set -e
# Echo each command
set -x

compiler=$1
build_type=$2

# set modules and install paths

openmpi_version=1.10.3

# compiler-specific settings
if [[ $compiler == "intel" ]]; then
  cxxmodule=intel/16.0.3
elif [[ $compiler == "gcc" ]]; then
  cxxmodule=gcc/5.3.0
fi
  
# build-type-specific settings
cmake_build_type=Release
extra_flags=
if [[ $build_type == "debug" ]]; then
  cmake_build_type=Debug
fi

export SHELL=/bin/sh

export MODULEPATH=""
. /opt/local/packages/Modules/default/init/sh
module load $cxxmodule
module load openmpi/${openmpi_version}
module load cmake # 3.0 or higher is required

echo $WORKSPACE
cd $WORKSPACE

mkdir build
cd build

cmake \
  -D CMAKE_C_COMPILER=`which mpicc` \
  -D CMAKE_CXX_COMPILER=`which mpiCC` \
  -D CMAKE_BUILD_TYPE=$cmake_build_type \
  -D ENABLE_UNIT_TESTS=True \
  -D ENABLE_APP_TESTS=True \
  -D ENABLE_JENKINS_OUTPUT=True \
  -D ENABLE_MPI=True \
  -D ENABLE_MPI_CXX_BINDINGS=True \
  $extra_flags \
  ..
make -j2
ctest --output-on-failure

