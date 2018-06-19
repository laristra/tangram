[![Build Status](https://travis-ci.org/laristra/tangram.svg?branch=master)](https://travis-ci.org/laristra/tangram)
[![codecov.io](https://codecov.io/github/laristra/tangram/coverage.svg?branch=master)](https://codecov.io/github/laristra/tangram/tangram?branch=master)
[![Quality Gate](https://sonarqube.com/api/badges/gate?key=tangram%3A%2Fmaster)](https://sonarqube.com/dashboard?id=tangram%3A%2Fmaster)

# tangram

The tangram library provides a framework for interface reconstruction
in computational physics applications. Interface reconstruction is
facilitated through the use of user-supplied _wrappers_ around
meshes with their materials data. Interface reconstruction algorithms 
can be customized (e.g. VOF, MOF) and, through the wrappers, take 
advantage of hybrid parallelism (MPI+X).

## Getting Started

To obtain a copy of tangram and its submodules from GitHub, clone
recursively:

```sh
git clone --recursive https://github.com/laristra/tangram
```

If you are familiar with Docker, take a look at
our
[Dockerfile](https://github.com/laristra/tangram/blob/master/docker/Dockerfile) for
a working build environment.  In particular, the Dockerfile builds off
of
the [tangram-buildenv](https://github.com/laristra/tangram-buildenv)
Dockerfile, and uses
our
[travis.yml](https://github.com/laristra/tangram/blob/master/.travis.yml) file
with Travis CI.

### Prerequisites

Tangram uses standard C++11 features, so a fairly modern compiler is
needed.  We regularly test with Intel 18+ or GCC 6.4+.  Utilizing the
full capabilities of portage will require an MPI implementation; we
regularly test with OpenMPI 2.1.2+ The build system _requires_ CMake
version 3.0+.

The following libraries are also _required_ (see examples below):

- **__Either__** Boost (1.68.0+) **__or__** Thrust (1.6.0+):
  We wrap some features of either one of these packages.  If you would
  like to run with OpenMP or TBB threads, then you _must_ use Thrust.

Tangram provides wrappers for a few third-party mesh types.  Building
support for these is _optional_:

- [Jali](http://github.com/lanl/jali):

  We regularly test with verison 0.9.8.  You will need to set the
  `Jali_Dir` CMake variable if you wish to build support for Jali and
  its tests (see examples below).

### Installing

In the simplest case where you have the appropriate versions mentioned
above and Boost is in the usual locations that CMake
searches, then the build step is:

```sh
tangram $ mkdir build
tangram $ cd build
tangram/build $ cmake -DENABLE_APP_TESTS=True ..
tangram/build $ make
```

This compiles the serial code and about a dozen application tests.  To
run the tests, simply execute

```sh
tangram/build $ make test
```

If you wish to install the code into the `CMAKE_INSTALL_PREFIX` then
simply execute
```sh
tangram/build $ make install
```

# License

This project is licensed under a modified 3-clause BSD license - see
the [LICENSE](https://github.com/laristra/tangram/blob/master/LICENSE)
file for details.

# Release

This software has been approved for open source release and has been
assigned **LA-CC-17-133**.

# Example of the configure script:

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
XMOF2D_INSTALL_PREFIX=/usr/local/codes/ngc/private/xmof2d/529f2dcdbe4-intel-18.0.1
JALI_INSTALL_PREFIX=/usr/local/codes/ngc/private/jali/0.9.8-intel-18.0.1-openmpi-2.1.2
mkdir build
cd build
cmake \
  -D CMAKE_C_COMPILER=`which mpicc` \
  -D CMAKE_CXX_COMPILER=`which mpiCC` \
  -D CMAKE_BUILD_TYPE=Release \
  -D ENABLE_UNIT_TESTS=True \
  -D ENABLE_MPI=True \
  -D ENABLE_JENKINS_OUTPUT=True \
  -D Jali_DIR:FILEPATH=$JALI_INSTALL_PREFIX/lib \
  -D XMOF2D_DIR:FILEPATH=$XMOF2D_INSTALL_PREFIX/share/cmake \
  -D BOOST_ROOT:FILEPATH=$TPL_INSTALL_PREFIX \
  -D NGC_INCLUDE_DIR:FILEPATH=$NGC_INCLUDE_DIR \
  ..
make -j2
ctest --output-on-failure
```

---
