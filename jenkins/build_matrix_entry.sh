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

jali_version=1.0.4
xmof2d_version=0.9.5
lapack_version=3.8.0
thrust_version=1.8.3

export NGC=/usr/local/codes/ngc
ngc_include_dir=$NGC/private/include
thrust_dir=$NGC/private/thrust/${thrust_version}

# compiler-specific settings
if [[ $compiler == "intel18" ]]; then
  intel_version=18.0.1
  cxxmodule=intel/${intel_version}
  # openmpi version that libs were built against
  openmpi_version=2.1.2
  # openmpi module for compiling and linking
  mpi_module=openmpi/2.1.2
  jali_install_dir=$NGC/private/jali/${jali_version}-intel-${intel_version}-openmpi-${openmpi_version}
  xmof2d_install_dir=$NGC/private/xmof2d/${xmof2d_version}-intel-${intel_version}
  lapacke_dir=$NGC/private/lapack/${lapack_version}-patched-intel-${intel_version}
elif [[ $compiler == "gcc6" ]]; then
  gcc_version=6.4.0
  cxxmodule=gcc/${gcc_version}
  # openmpi version that libs were built against
  openmpi_version=2.1.2
  # openmpi module for compiling and linking
  mpi_module=openmpi/2.1.2
  jali_install_dir=$NGC/private/jali/${jali_version}-gcc-${gcc_version}-openmpi-${openmpi_version}
  xmof2d_install_dir=$NGC/private/xmof2d/${xmof2d_version}-gcc-${gcc_version}
  lapacke_dir=$NGC/private/lapack/${lapack_version}-patched-gcc-${gcc_version}
elif [[ $compiler == "gcc7" ]]; then
  gcc_version=7.3.0
  cxxmodule=gcc/${gcc_version}
  # openmpi version that libs were built against
  openmpi_version=2.1.2
  # openmpi module for compiling and linking
  mpi_module=openmpi/2.1.2
  jali_install_dir=$NGC/private/jali/${jali_version}-gcc-${gcc_version}-openmpi-${openmpi_version}
  xmof2d_install_dir=$NGC/private/xmof2d/${xmof2d_version}-gcc-${gcc_version}
  lapacke_dir=$NGC/private/lapack/${lapack_version}-patched-gcc-${gcc_version}
elif [[ $compiler == "gcc8" ]]; then
  gcc_version=8.2.0
  cxxmodule=gcc/${gcc_version}
  # openmpi version that libs were built against
  openmpi_version=3.1.3
  # openmpi module for compiling and linking
  mpi_module=openmpi/3.1.3
  jali_install_dir=$NGC/private/jali/${jali_version}-gcc-${gcc_version}-openmpi-${openmpi_version}
  xmof2d_install_dir=$NGC/private/xmof2d/${xmof2d_version}-gcc-${gcc_version} 
  lapacke_dir=$NGC/private/lapack/${lapack_version}-patched-gcc-${gcc_version} 
fi

jali_flags="-D Jali_DIR:FILEPATH=$jali_install_dir/lib"
xmof2d_flags="-D XMOF2D_DIR:FILEPATH=$xmof2d_install_dir/share/cmake"
mpi_flags="-D ENABLE_MPI=True"
lapacke_flags="-D LAPACKE_DIR:FILEPATH=$lapacke_dir"

# build-type-specific settings
cmake_build_type=Release
mpi_flags="-D ENABLE_MPI=True"
thrust_flags=
if [[ $build_type == "debug" ]]; then
    cmake_build_type=Debug
elif [[ $build_type == "serial" ]]; then
    mpi_flags=
    jali_flags=    # jali is not available in serial
elif [[ $build_type == "thrust" ]]; then
    thrust_flags="_D ENABLE_THRUST=True -DTHRUST_DIR:FILEPATH=${thrust_dir}"
fi

export SHELL=/bin/sh

export MODULEPATH=""
. /opt/local/packages/Modules/default/init/sh
module load $cxxmodule
module load ${mpi_module}
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
  -D NGC_INCLUDE_DIR:FILEPATH=$ngc_include_dir \
  -D BOOST_ROOT:FILEPATH=$boost_dir \
  $mpi_flags \
  $jali_flags \
  $xmof2d_flags \
  $lapacke_flags \
  $thrust_flags \
  ..
make -j2
ctest --output-on-failure
