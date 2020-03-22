#!/usr/bin/env bash
# This script is executed on Jenkins using
#
#     bash $WORKSPACE/jenkins/build_pr.sh
#
# The exit code determines if the test succeeded or failed.

# Exit on error
set -e
# Echo each command
set -x

#Setup proxies
export http_proxy="http://proxyout.lanl.gov:8080"
export https_proxy="http://proxyout.lanl.gov:8080"
export ftp_proxy="http://proxyout.lanl.gov:8080"
export HTTP_PROXY="http://proxyout.lanl.gov:8080"
export HTTPS_PROXY="http://proxyout.lanl.gov:8080"
export FTP_PROXY="http://proxyout.lanl.gov:8080"

# Put a couple of settings in place to generate test output even if
# the README doesn't ask for it.
export CTEST_OUTPUT_ON_FAILURE=1
CACHE_OPTIONS="-D ENABLE_JENKINS_OUTPUT=True"
sed "s/^ *cmake/& $CACHE_OPTIONS/g" $WORKSPACE/README.md >$WORKSPACE/README.md.1

# Run build/test commands from README
python $WORKSPACE/jenkins/parseREADME.py $WORKSPACE/README.md.1 $WORKSPACE

# General NGC include directory
NGC_DIR=/usr/local/codes/ngc/private
NGC_INCLUDE_DIR=${NGC_DIR}/include

git config user.email ""
git config user.name "Jenkins"

export SHELL=/bin/sh

export MODULEPATH=""
. /opt/local/packages/Modules/default/init/sh
intel_version=18.0.1
module load intel/${intel_version}
openmpi_version=2.1.2
module load openmpi/${openmpi_version}
module load cmake/3.14.0  # Need 3.13 or higher

echo $WORKSPACE
cd $WORKSPACE
git clean -dfx

# Build XMOF2D

git clone https://github.com/laristra/XMOF2D.git xmof2d-repo
cd xmof2d-repo
XMOF2D_INSTALL_PREFIX=`pwd`/xmof2d-inst
mkdir build
cd build

cmake \
cmake \
  -D CMAKE_BUILD_TYPE=Release \
  -D CMAKE_C_COMPILER=`which mpicc` \
  -D CMAKE_CXX_COMPILER=`which mpiCC` \
  -D CMAKE_CXX_FLAGS='-std=c++11' \
  -D INSTALL_DIR:FILEPATH=$XMOF2D_INSTALL_PREFIX \
  ..
make -j2
ctest -j2 --output-on-failure
make install


# Wonton

WONTON_INSTALL_PREFIX=$NGC_DIR/wonton/new-cmake-intel-${intel_version}-openmpi-${openmpi-version}

# Build Tangram with Thrust

cd $WORKSPACE
mkdir build
cd build

cmake \
  -D CMAKE_BUILD_TYPE=Release \
  -D ENABLE_UNIT_TESTS=True \
  -D ENABLE_MPI=True \
  -D ENABLE_JENKINS_OUTPUT=True \
  -D ENABLE_XMOF2D=True
  -D XMOF2D_ROOT:FILEPATH=$XMOF2D_INSTALL_PREFIX/share/cmake \
  -D WONTON_ROOT:FILEPATH=$TPL_INSTALL_PREFIX \
  -D ENABLE_THRUST=True \
  ..
make -j2
ctest --output-on-failure
for test_name in intersectClipper test_matfuncs; do
  valgrind --gen-suppressions=all --suppressions=../jenkins/portage_valgrind.supp --error-exitcode=1 --leak-check=full test/$test_name
done


# Build Tangram without Thrust

cd $WORKSPACE
mkdir build-nothrust
cd build-nothrust

cmake \
  -D CMAKE_BUILD_TYPE=Release \
  -D ENABLE_UNIT_TESTS=True \
  -D ENABLE_MPI=True \
  -D XMOF2D_ROOT:FILEPATH=$XMOF2D_INSTALL_PREFIX/share/cmake \
  -D ENABLE_THRUST=False \
  ..
make -j2
ctest --output-on-failure
