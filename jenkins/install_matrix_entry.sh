#!/usr/bin/env bash
# This script is executed on Jenkins using
#
#     $WORKSPACE/jenkins/install_matrix_entry.sh <compiler> <build_type>
#
# The exit code determines if the test succeeded or failed.
# Note that the environment variable WORKSPACE must be set (Jenkins
# will do this automatically).

# Exit on error
set -e
# Echo each command
set -x

# Set umask so installation directories will have rwx permissions for group
umask 007

compiler=$1
build_type=$2

# set modules and install paths
export NGC="/usr/local/codes/ngc"

# versions
xmof2d_version=0.9.5
wonton_version=$3
tangram_version=$4

# suffixes
compiler_suffix=""
mpi_suffix=""
thrust_suffix=""

# flags
mpi_flags=""
thrust_flags=""
jali_flags=""
flecsi_flags=""

# compiler-specific settings
if [[ $compiler =~ "intel" ]]; then

  compiler_version=18.0.1
  cxxmodule=intel/${compiler_version}
  compiler_suffix="-intel-${compiler_version}"

  openmpi_version=2.1.2
  mpi_module=openmpi/2.1.2
  mpi_suffix="-openmpi-${openmpi_version}"

elif [[ $compiler =~ "gcc" ]]; then

  openmpi_version=2.1.2
  if [[ $compiler == "gcc6" ]]; then
    compiler_version=6.4.0
  elif [[ $compiler == "gcc7" ]]; then
    compiler_version=7.3.0
  fi

  cxxmodule=gcc/${compiler_version}
  compiler_suffix="-gcc-${compiler_version}"

  mpi_module=openmpi/${openmpi_version}
  mpi_suffix="-openmpi-${openmpi_version}"
fi

# build-type-specific settings
mpi_flags="-D TANGRAM_ENABLE_MPI=True"
if [[ $build_type == "serial" ]]; then
  mpi_flags=""
  mpi_suffix=""
fi

cmake_build_type=Release
if [[ $build_type == "debug" ]]; then
  cmake_build_type=Debug
fi

if [[ $build_type == "thrust" ]]; then
  thrust_flags="-D TANGRAM_ENABLE_THRUST=True"
  thrust_suffix="-thrust"
fi

if [[ $compiler == "gcc6" && $build_type != "serial" ]]; then
    flecsi_flags="-D TANGRAM_ENABLE_FleCSI:BOOL=True"  # FleCSI found through Wonton
fi
if [[ $build_type != "serial" ]]; then
    jali_flags="-D TANGRAM_ENABLE_Jali:BOOL=True"  # Jali found through Wonton
fi

# install paths
xmof2d_install_dir="$NGC/private/xmof2d/${xmof2d_version}${compiler_suffix}"
wonton_install_dir="$NGC/private/wonton/${wonton_version}${compiler_suffix}${mpi_suffix}${thrust_suffix}"
tangram_install_dir="$NGC/private/tangram"/${tangram_version}${compiler_suffix}${mpi_suffix}${thrust_suffix}

xmof2d_flags="-D TANGRAM_ENABLE_XMOF2D=True -D XMOF2D_ROOT:FILEPATH=$xmof2d_install_dir/share/cmake"
wonton_flags="-D WONTON_ROOT:FILEPATH=$wonton_install_dir"

export SHELL=/bin/sh
export MODULEPATH=""
. /opt/local/packages/Modules/default/init/sh
module load $cxxmodule
if [[ -n "$mpi_flags" ]]; then
  module load ${mpi_module}
fi
module load cmake/3.14.0 # 3.13 or higher is required
module load git

echo "Jenkins workspace: $WORKSPACE"
cd $WORKSPACE

rm -rf build
mkdir -p build
cd build

cmake \
  -D CMAKE_BUILD_TYPE=$cmake_build_type \
  -D CMAKE_INSTALL_PREFIX=$tangram_install_dir \
  -D ENABLE_UNIT_TESTS=True \
  -D ENABLE_APP_TESTS=True \
  -D ENABLE_JENKINS_OUTPUT=True \
  $mpi_flags \
  $wonton_flags \
  $xmof2d_flags \
  $thrust_flags \
  $jali_flags \
  $flecsi_flags \
  ..
make -j2
ctest --output-on-failure
make install
