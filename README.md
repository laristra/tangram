* Welcome to tangram!

To obtain the code, simply clone it recursively from our BitBucket Server:

```git clone --recursive ssh://git@xcp-stash.lanl.gov:7999/laristra/tangram.git```

Example of the configure script:

```c++
#!/bin/bash

BUILD_TYPE=Release
JALI_INSTALL_PREFIX=/path/to/jali/installation
TPLS_INSTALL_PREFIX=/path/to/jali/tpls/installation
XMOF2D_INSTALL_PREFIX=/path/to/XMOF2D/installation
THRUST_PATH=/path/to/thrust
BOOST_PATH=/path/to/boost
TCMALLOC_PATH=/path/to/google/performance/tools/installation

CC=`which mpicc`
CXX=`which mpiCC`

rm -rf build
mkdir build
cd build

cmake \
    -D CMAKE_C_COMPILER=${CC} \
    -D CMAKE_CXX_COMPILER=${CXX} \
    -D CMAKE_BUILD_TYPE=${BUILD_TYPE} \
    -D ENABLE_UNIT_TESTS=True \
    -D ENABLE_APP_TESTS=True \
    -D ENABLE_MPI=True \
    -D ENABLE_MPI_CXX_BINDINGS=True \
    -D XMOF2D_DIR:FILEPATH=${XMOF2D_INSTALL_PREFIX}/lib \
    -D BOOST_ROOT:FILEPATH=${BOOST_PATH} \
    -D ENABLE_THRUST:BOOL=True \
    -D THRUST_DIR:PATH=${THRUST_PATH} \
    -D Jali_DIR:FILEPATH=${JALI_INSTALL_PREFIX}/lib \
    -D ENABLE_TCMALLOC=True \
    -D TCMALLOC_LIB:FILEPATH={TCMALLOC_PATH}/libtcmalloc.so \
    ..
make -j
```
# Example builds

Below we list copy & paste instructions for several local machines; we
have a script that parses this README file to execute the examples
below to ensure they build.

## Varan

Execute the following from the tangram root directory:

```c++
# machine=varan
export MODULEPATH=""
. /opt/local/packages/Modules/default/init/sh
module load intel/17.0.1 openmpi/1.10.5 cmake/3.8.2
TPL_INSTALL_PREFIX=/usr/local/codes/ngc/private/jali-tpl/1.0.9-intel-17.0.1-openmpi-1.10.5
NGC_INCLUDE_DIR=/usr/local/codes/ngc/private/include

TWORKSPACE=`pwd`
git clone --recursive ssh://git@xcp-stash.lanl.gov:7999/laristra/tangram.git
cd $TWORKSPACE/
mkdir build

cd $TWORKSPACE/build
git clone ssh://git@xcp-stash.lanl.gov:7999/laristra/jali.git jali-repo
cd jali-repo
JALI_INSTALL_PREFIX=`pwd`/jali-inst
mkdir build
cd build
cmake  \
 -C $TPL_INSTALL_PREFIX/share/cmake/Jali-tpl-config.cmake  \
 -D CMAKE_BUILD_TYPE=Release   -D CMAKE_CXX_FLAGS='-std=c++11'  \
 -D CMAKE_INSTALL_PREFIX:FILEPATH=$JALI_INSTALL_PREFIX  \
 -D HDF5_NO_SYSTEM_PATHS:BOOL=TRUE  \
 -D BOOST_ROOT:FILEPATH=$TPL_INSTALL_PREFIX  \
 -D ENABLE_MSTK_Mesh:BOOL=TRUE  \
 -D ENABLE_STK_Mesh:BOOL=FALSE  \
 -D ENABLE_MOAB_Mesh:BOOL=FALSE   ..
make -j2
ctest -j2 --output-on-failure
make install

cd $TWORKSPACE/build
git clone ssh://git@xcp-stash.lanl.gov:7999/~rgertl/xmof2d.git xmof2d-repo
cd xmof2d-repo
XMOF2D_INSTALL_PREFIX=`pwd`/xmof2d-inst
mkdir build
cd build
cmake  \
 -D CMAKE_BUILD_TYPE=Release  \
 -D CMAKE_C_COMPILER=`which mpicc`  \
 -D CMAKE_CXX_COMPILER=`which mpiCC`  \
 -D CMAKE_CXX_FLAGS='-std=c++11'  \
 -D INSTALL_DIR:FILEPATH=$XMOF2D_INSTALL_PREFIX  \
 -D INSTALL_ADD_VERSION=yes  \
 -D XMOF2D_VERSION_MAJOR=0  \
 -D XMOF2D_VERSION_MINOR=9  \
 -D INSTALL_PREFIX_ARCHOS=no  \
 ..
make -j2
ctest -j2 --output-on-failure
make install

cd $TWORKSPACE/build
cmake \
  -D CMAKE_C_COMPILER=`which mpicc` \
  -D CMAKE_CXX_COMPILER=`which mpiCC` \
  -D CMAKE_BUILD_TYPE=Release \
  -D ENABLE_UNIT_TESTS=True \
  -D ENABLE_MPI=True \
  -D ENABLE_MPI_CXX_BINDINGS=True \
  -D ENABLE_JENKINS_OUTPUT=True \
  -D Jali_DIR:FILEPATH=$JALI_INSTALL_PREFIX/lib \
  -D XMOF2D_DIR:FILEPATH=$XMOF2D_INSTALL_PREFIX/lib \
  -D BOOST_ROOT:FILEPATH=$TPL_INSTALL_PREFIX \
  -D NGC_INCLUDE_DIR:FILEPATH=$NGC_INCLUDE_DIR \
  -D ENABLE_THRUST=True \
  ..
make -j2
ctest --output-on-failure
```

---
