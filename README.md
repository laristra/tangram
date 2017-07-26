* Welcome to tangram!

To obtain the code, simply clone it recursively from our BitBucket Server:

```git clone --recursive ssh://git@xcp-stash.lanl.gov:7999/laristra/tangram.git```

Example of the configure script:

```c++
#!/bin/bash

BUILD_TYPE=Release
BASE_PATH=
JALI_INSTALL_PREFIX=
TPLS_INSTALL_PREFIX=
XMOF2D_INSTALL_PREFIX=${BASE_PATH}/XMOF2D/install-${BUILD_TYPE}
THRUST_PATH=
BOOST_PATH=
TCMALLOC_PATH=

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
    -D XMOF2D_DIR:FILEPATH=$XMOF2D_INSTALL_PREFIX/lib \
    -D BOOST_ROOT:FILEPATH=${BOOST_PATH} \
    -D ENABLE_THRUST:BOOL=True \
    -D THRUST_DIR:PATH=$THRUST_PATH \
    -D Jali_DIR:FILEPATH=$JALI_INSTALL_PREFIX/lib \
    -D ENABLE_TCMALLOC=True \
    -D TCMALLOC_LIB:FILEPATH=TCMALLOC_PATH/libtcmalloc.so \
    ..
make -j
```
