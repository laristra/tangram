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

# special case for README builds
if [[ $build_type == "readme" ]]; then
  python2 $WORKSPACE/jenkins/parseREADME.py $WORKSPACE/README.md $WORKSPACE
  exit
fi

# set modules and install paths

xmof2d_version=0.9.5

export NGC=/usr/local/codes/ngc
ngc_include_dir=$NGC/private/include

# compiler-specific settings
if [[ $compiler == "intel18" ]]; then

    intel_version=18.0.1
    cxxmodule=intel/${intel_version}
    compiler_suffix="-intel-${intel_version}"
    
    openmpi_version=2.1.2
    mpi_module=openmpi/${openmpi_version}
    mpi_suffix="-openmpi-${openmpi_version}"

elif [[ $compiler =~ "gcc" ]]; then

    openmpi_version=2.1.2
    if [[ $compiler == "gcc6" ]]; then
	gcc_version=6.4.0
    elif [[ $compiler == "gcc7" ]]; then
	gcc_version=7.4.0
    elif [[ $compiler == "gcc8" ]]; then
	gcc_version=8.2.0
	openmpi_version=3.1.3
    fi
    
    cxxmodule=gcc/${gcc_version}
    compiler_suffix="-gcc-${gcc_version}"
    
    mpi_module=openmpi/${openmpi_version}
    mpi_suffix="-openmpi-${openmpi_version}"
fi


# build-type-specific settings
mpi_flags="-D ENABLE_MPI=True"
if [[ $build_type == "serial" ]]; then
    mpi_flags=
    mpi_suffix=
fi

cmake_build_type=Release
if [[ $build_type == "debug" ]]; then
    cmake_build_type=Debug
fi

thrust_flags=
thrust_suffix=
if [[ $build_type == "thrust" ]]; then
    thrust_flags="-D ENABLE_THRUST=True"
    thrust_suffix="-thrust"
fi



xmof2d_install_dir=$NGC/private/xmof2d/${xmof2d_version}${compiler_suffix}
xmof2d_flags="-D ENABLE_XMOF2D=True -D XMOF2D_ROOT:FILEPATH=$xmof2d_install_dir/share/cmake"

wonton_install_dir=$NGC/private/wonton/new-cmake${compiler_suffix}${mpi_suffix}${thrust_suffix}
wonton_flags="-D WONTON_ROOT:FILEPATH=$wonton_install_dir"

export SHELL=/bin/sh

export MODULEPATH=""
. /opt/local/packages/Modules/default/init/sh
module load $cxxmodule
module load ${mpi_module}
module load cmake/3.14.0  # 3.13 or higher is required

echo $WORKSPACE
cd $WORKSPACE

mkdir build
cd build

cmake \
  -D CMAKE_BUILD_TYPE=$cmake_build_type \
  -D ENABLE_UNIT_TESTS=True \
  -D ENABLE_APP_TESTS=True \
  -D ENABLE_JENKINS_OUTPUT=True \
  $mpi_flags \
  $wonton_flags \
  $xmof2d_flags \
  $thrust_flags \
  ..
make -j2
ctest --output-on-failure
