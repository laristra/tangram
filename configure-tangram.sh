#!/bin/bash

BUILD_TYPE=Release
BASE_PATH=/Users/kikinzon/research
JALI_INSTALL_PREFIX=${BASE_PATH}/jali/inst-jali
TPL_INSTALL_PREFIX=${BASE_PATH}/jali/inst-tpl
XMOF2D_INSTALL_PREFIX=${BASE_PATH}/XMOF2D/install-${BUILD_TYPE}
THRUST_PATH=${BASE_PATH}/thrust
BOOST_PATH=${BASE_PATH}/boost

CC=/usr/bin/mpicc
CXX=/usr/bin/mpic++
#CC=/usr/local/bin/mpicc
#CXX=/usr/local/bin/mpic++
#CC=/usr/local/opt/llvm/bin/clang
#CXX=${CC}++

rm -rf build
mkdir build
cd build

cmake \
    -D CMAKE_C_COMPILER=${CC} \
    -D CMAKE_CXX_COMPILER=${CXX} \
    -D CMAKE_BUILD_TYPE=${BUILD_TYPE} \
    -D ENABLE_UNIT_TESTS=False \
    -D ENABLE_APP_TESTS=True \
    -D ENABLE_MPI=True \
    -D ENABLE_MPI_CXX_BINDINGS=True \
    -D XMOF2D_DIR:FILEPATH=$XMOF2D_INSTALL_PREFIX/lib \
    -D BOOST_ROOT:FILEPATH=${BOOST_PATH} \
    -D ENABLE_THRUST:BOOL=False \
    -D THRUST_DIR:PATH=$THRUST_PATH \
    -D Jali_DIR:FILEPATH=$JALI_INSTALL_PREFIX/lib \
    .. 
make
