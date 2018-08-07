#!/bin/bash

BUILD_TYPE=Debug
BASE_PATH=/Users/kikinzon/research
TANGRAM_DIR=${BASE_PATH}/tangram-fork
JALI_INSTALL_PREFIX=${BASE_PATH}/jali/inst-jali
TPL_INSTALL_PREFIX=${BASE_PATH}/jali/inst-tpl
XMOF2D_INSTALL_PREFIX=${BASE_PATH}/github-xmof2d/install-${BUILD_TYPE}
THRUST_PATH=${BASE_PATH}/thrust
BOOST_PATH=${BASE_PATH}/boost
LAPACKE_PATH=${BASE_PATH}/lapack/install-Release

BUILD_DIR=build-${BUILD_TYPE}
INSTALL_DIR=install-${BUILD_TYPE}

CC=/usr/local/bin/mpicc
CXX=/usr/local/bin/mpic++

#rm -rf ${BUILD_DIR}
#mkdir ${BUILD_DIR}
cd ${BUILD_DIR}

cmake \
    -D CMAKE_C_COMPILER=${CC} \
    -D CMAKE_CXX_COMPILER=${CXX} \
    -D CMAKE_BUILD_TYPE=${BUILD_TYPE} \
    -D CMAKE_INSTALL_PREFIX=${TANGRAM_DIR}/${INSTALL_DIR} \
    -D ENABLE_UNIT_TESTS=True \
    -D ENABLE_APP_TESTS=True \
    -D ENABLE_MPI=True \
    -D Jali_DIR:FILEPATH=$JALI_INSTALL_PREFIX/lib \
    -D XMOF2D_DIR:FILEPATH=$XMOF2D_INSTALL_PREFIX/share/cmake \
    -D BOOST_ROOT:FILEPATH=${BOOST_PATH} \
    -D ENABLE_THRUST:BOOL=True \
    -D THRUST_DIR:PATH=$THRUST_PATH \
    -D ENABLE_TCMALLOC=True \
    -D TCMALLOC_LIB:FILEPATH=${BASE_PATH}/tcmalloc/lib/libtcmalloc.a \
    -D LAPACKE_DIR:FILEPATH=${LAPACKE_PATH} \
    .. 
make -j4 install
#ctest -j4 --output-on-failure
# -D Jali_DIR:FILEPATH=$JALI_INSTALL_PREFIX/lib \
