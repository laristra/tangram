[![Build Status](https://travis-ci.org/laristra/tangram.svg?branch=master)](https://travis-ci.org/laristra/tangram)
[![codecov.io](https://codecov.io/github/laristra/tangram/coverage.svg?branch=master)](https://codecov.io/github/laristra/tangram/tangram?branch=master)
[![Quality Gate](https://sonarqube.com/api/badges/gate?key=tangram%3A%2Fmaster)](https://sonarqube.com/dashboard?id=tangram%3A%2Fmaster)

# Tangram

The Tangram library provides a framework for interface reconstruction
in computational physics applications. Interface reconstruction is
facilitated through the use of user-supplied _wrappers_ around
meshes with their materials data. Interface reconstruction algorithms 
can be customized (e.g. VOF, MOF) and, through the wrappers, take 
advantage of hybrid parallelism (MPI+X).

## Getting Started

To obtain a copy of Tangram and its submodules from GitHub, clone
recursively:

```sh
git clone --recursive https://github.com/laristra/tangram
```

If you are familiar with Docker, take a look at
our
[Dockerfile](https://github.com/laristra/tangram/blob/master/docker/Dockerfile) for
a working build environment.  In particular, the Dockerfile builds off
of
the [portage-buildenv](https://github.com/laristra/portage-buildenv)
Dockerfile, and uses
our
[travis.yml](https://github.com/laristra/tangram/blob/master/.travis.yml) file
with Travis CI.

### Prerequisites

Tangram uses standard C++11 features, so a fairly modern compiler is
needed.  We regularly test with Intel 18.0.1, GCC 6.4.0, and GCC 7.3.0.  
Utilizing the full capabilities of Tangram will require an MPI implementation; 
we regularly test with OpenMPI 2.1.2. The build system _requires_ CMake
version 3.13+.

The following libraries are also _required_ (see examples below):

- **__Either__** Boost (1.68.0+) **__or__** Thrust (1.8.0+):
  We wrap some features of either one of these packages.  If you would
  like to run with OpenMP, then you _must_ use Thrust.

Tangram requires the [Wonton](https://github.com/laristra/wonton)
library. The path where Wonton is installed must be specified either
through the **WONTON_ROOT** variable or included in the
**CMAKE_PREFIX_PATH**. Tangram will pick up the third party libraries
that Wonton is built with (Jali, FleCSI and LAPACKE) based on the
installation. Some options requested for Tangram must match those of
the Wonton installation. For example, you cannot ask for Tangram to be
built with MPI (or Thrust) and link to an installation of Wonton built
without MPI (or Thrust).

Tangram can optionally be built with the
[XMOF2D](https://github.com/laristra/xmof2d) library. Enable XMOF2D
support using the option "TANGRAM_ENABLE_XMOF2D=True" and point to the
location where *XMOF2DConfig.cmake* is installed either through the
**XMOF2D_ROOT** variable or include it in the **CMAKE_PREFIX_PATH**.

### Installing

In the simplest case where you have the appropriate versions mentioned
above and LAPACKE and Boost are in the usual locations that CMake
searches, then the build step is:

```sh
tangram $ mkdir build
tangram $ cd build
tangram/build $ cmake -DWONTON_ROOT:FILEPATH=/path/to/wonton ..
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
XMOF2D_INSTALL_PREFIX=/path/to/XMOF2D/installation

rm -rf build
mkdir build
cd build

cmake \
    -D CMAKE_BUILD_TYPE=${BUILD_TYPE} \
	-D CMAKE_INSTALL_PREFIX=/path/where/to/install/tangram \
    -D ENABLE_UNIT_TESTS=True \
    -D ENABLE_APP_TESTS=True \
    -D TANGRAM_ENABLE_MPI=True \
	-D WONTON_ROOT=/path/where/wonton/is/installed \
	-D TANGRAM_ENABLE_Jali=True \
	-D TANGRAM_ENABLE_XMOF2D=True \
    -D XMOF2D_DIR:FILEPATH=${XMOF2D_INSTALL_PREFIX}/share/cmake \
    -D TANGRAM_ENABLE_THRUST:BOOL=True \
    ..
make -j
```
# Example builds

Below we list copy & paste instructions for several local machines; we
have a script that parses this README file to execute the examples
below to ensure they build.

## Varan

Execute the following from the Tangram root directory:

```c++
# machine=varan
export MODULEPATH=""
. /opt/local/packages/Modules/default/init/sh
module load intel/18.0.1 openmpi/2.1.2 cmake/3.14.0
XMOF2D_INSTALL_PREFIX=/usr/local/codes/ngc/private/xmof2d/0.9.5-intel-18.0.1
WONTON_INSTALL_PREFIX=/usr/local/codes/ngc/private/wonton/1.2.2-intel-18.0.1-openmpi-2.1.2
mkdir build
cd build
cmake \
  -D CMAKE_C_COMPILER=`which mpicc` \
  -D CMAKE_CXX_COMPILER=`which mpiCC` \
  -D CMAKE_BUILD_TYPE=Release \
  -D ENABLE_UNIT_TESTS=True \
  -D TANGRAM_ENABLE_MPI=True \
  -D ENABLE_JENKINS_OUTPUT=True \
  -D WONTON_ROOT=$WONTON_INSTALL_PREFIX \
  -D TANGRAM_ENABLE_Jali=True \
  -D TANGRAM_ENABLE_XMOF2D=True \
  -D XMOF2D_ROOT:FILEPATH=$XMOF2D_INSTALL_PREFIX/share/cmake \
  ..
make -j2
ctest --output-on-failure
```

---
