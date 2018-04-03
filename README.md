* Welcome to tangram!

To obtain the code, simply clone it recursively from our BitBucket Server:

```
git clone --recursive ssh://git@xcp-stash.lanl.gov:7999/laristra/tangram.git
```

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
    -D XMOF2D_DIR:FILEPATH=${XMOF2D_INSTALL_PREFIX}/share/cmake \
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
module load intel/18.0.1 openmpi/2.1.2 cmake/3.10.2
TPL_INSTALL_PREFIX=/usr/local/codes/ngc/private/jali-tpl/1.0.9-intel-18.0.1-openmpi-2.1.2
NGC_INCLUDE_DIR=/usr/local/codes/ngc/private/include
XMOF2D_INSTALL_PREFIX=/usr/local/codes/ngc/private/xmof2d/6023dea445c-intel-18.0.1
JALI_INSTALL_PREFIX=/usr/local/codes/ngc/private/jali/0.9.8-intel-18.0.1-openmpi-2.1.2
mkdir build
cd build
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
  ..
make -j2
ctest --output-on-failure
```

---
