#!/usr/bin/env bash
# This script is executed on Jenkins using
#
# $WORKSPACE/jenkins/build_matrix_entry_varan.sh BUILD_TYPE <VER> <WONTON_VER>
#
# BUILD_TYPE  -  pr, nightly, install
#
# if VER is abset, the HEAD of the master branch will be taken. If
# WONTON_VER is absent, the HEAD of the master branch of wonton will
# be taken. If BUILD_TYPE is 'install' and VER is specified, it will
# install it to /install_prefix/tangram/$VER-blah-blah; if VER is not
# specified, it will install to /install_prefix/wonton/dev-blah-blah
#
# Note that the following environment variables must be set (Jenkins
# will do this automatically).
#
# WORKSPACE   -  where the code is checked out
# CONFIG_TYPE -  base, debug, serial, readme, thrust, kokkos
# COMPILER    -  intel, gcc6, gcc7
# BRANCH_NAME -  master, kokkos
#
# The exit code determines if the test succeeded or failed.

# Exit on error
set -e
# Echo each command
set -x

# Set umask so installation directories will have rwx permissions for group
umask 007


BUILD_TYPE=$1
version=$2
if [[ $version == "" ]]; then
    version=dev
fi
wonton_version=$3
if [[ $wonton_version == "" ]]; then
    wonton_version=dev
fi

echo "inside build_matrix on PLATFORM=$PLATFORM with BUILD_TYPE=$BUILD_TYPE $CONFIG_TYPE=$CONFIG_TYPE COMPILER=$COMPILER"



# special case for README builds
if [[ $BUILD_TYPE != "install" && $CONFIG_TYPE == "readme" ]]; then

   # Put a couple of settings in place to generate test output even if
   # the README doesn't ask for it.
   export CTEST_OUTPUT_ON_FAILURE=1
   CACHE_OPTIONS="-D ENABLE_JENKINS_OUTPUT=True"
   sed "s/^ *cmake/& $CACHE_OPTIONS/g" $WORKSPACE/README.md >$WORKSPACE/README.md.1
   python2 $WORKSPACE/jenkins/parseREADME.py \
	   $WORKSPACE/README.md.1 \
	   $WORKSPACE \
	   varan
   exit
   
fi

# set modules and install paths

xmof2d_version=0.9.5

export NGC=/usr/local/codes/ngc
ngc_include_dir=$NGC/private/include

# compiler-specific settings
if [[ $COMPILER =~ "intel" ]]; then

    compiler_version=18.0.1
    cxxmodule=intel/${compiler_version}
    compiler_suffix="-intel-${compiler_version}"
    
    openmpi_version=2.1.2
    mpi_module=openmpi/${openmpi_version}
    mpi_suffix="-openmpi-${openmpi_version}"

elif [[ $COMPILER =~ "gcc" ]]; then

    openmpi_version=2.1.2
    if [[ $COMPILER == "gcc6" ]]; then
	compiler_version=6.4.0
    elif [[ $COMPILER == "gcc7" ]]; then
	compiler_version=7.3.0
    elif [[ $COMPILER == "gcc8" ]]; then
	compiler_version=8.2.0
	openmpi_version=3.1.3
    fi
    
    cxxmodule=gcc/${compiler_version}
    compiler_suffix="-gcc-${compiler_version}"
    
    mpi_module=openmpi/${openmpi_version}
    mpi_suffix="-openmpi-${openmpi_version}"
fi


# XMOF2D
xmof2d_install_dir=$NGC/private/xmof2d/${xmof2d_version}${compiler_suffix}
xmof2d_flags="-D TANGRAM_ENABLE_XMOF2D=True -D XMOF2D_ROOT:FILEPATH=$xmof2d_install_dir/share/cmake"

# Jali
jali_flags="-D TANGRAM_ENABLE_Jali::BOOL=True"

# FleCSI
flecsi_flags=
if [[ $compiler == "gcc6" && $build_type != "serial" ]]; then
    flecsi_flags="-D TANGRAM_ENABLE_FleCSI:BOOL=True"  # FleCSI found through Wonton
fi



# THRUST
thrust_flags=
thrust_suffix=
if [[ $CONFIG_TYPE == "thrust" ]]; then
    thrust_flags="-D TANGRAM_ENABLE_THRUST=True"
    thrust_suffix="-thrust"
fi

# MPI or not
mpi_flags="-D TANGRAM_ENABLE_MPI=True"
if [[ $CONFIG_TYPE == "serial" ]]; then
    mpi_flags="-D TANGRAM_ENABLE_MPI=False"
    mpi_suffix=
    jali_flags=
    flecsi_flags=
fi

# Debug or Optimized build
cmake_build_type=Release
debug_suffix=
if [[ $CONFIG_TYPE == "debug" ]]; then
    cmake_build_type=Debug
    debug_suffix="-debug"
fi

# WONTON
wonton_install_dir=$NGC/private/wonton/${wonton_version}${compiler_suffix}${mpi_suffix}${kokkos_suffix}${thrust_suffix}${debug_suffix}
wonton_flags="-D WONTON_ROOT:FILEPATH=$wonton_install_dir"

# Build up an install dir name
tangram_install_dir=$NGC/private/tangram/${version}${compiler_suffix}${mpi_suffix}${thrust_suffix}${kokkos_suffix}${debug_suffix}


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
  -D CMAKE_CXX_FLAGS="-Wall -Werror" \
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

if [[ $BUILD_TYPE == "install" ]]; then
    make install
fi
