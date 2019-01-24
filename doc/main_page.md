# Welcome to Tangram!   {#mainpage}

Tangram is a framework that computational physics applications can use
to build a highly customized, hybrid parallel (MPI+X) interface
reconstruction library for constructing in-cell material polytopes
based on given volume fractions and, optionally, centroids data.

We aim to provide:
- A modern, modular design - pick and choose your preferred
  intersection and interface reconstruction methods
- Support for general polytopal meshes and higher-order material moments
- High flexibility for application customization
- Algorithms that take advantage of both distributed and on-node parallelism
- An _Open Source Community_ for these tools!
- Use of client application's native mesh data structures

See the [Concepts](@ref concepts) page for a high-level discussion of
the components of Tangram.

See the [Example Use](@ref example) page for a simple example of
hooking Tangram up to a mesh and interface reconstructor.

---

# Details and Requirements

At a minimum, Tangram requires:
- A C++-11 compatible compiler; regular testing is performed with GCC
  6.3+ and Intel 17+.
- [CMake](https://cmake.org) 3.8+
- [LAPACKE](https://https://github.com/Reference-LAPACK/lapack/tree/master/LAPACKE) 3.8.0+
- [Boost](https://www.boost.org) 1.68.0+ **or** [Thrust](https://thrust.github.io/) 1.8.1+

In addition to the minimum set of libraries, Tangram has been known to
work with version 1.0.0 of the [Jali](https://github.com/lanl/jali).

Distributed parallelism of Tangram is currently supported through MPI;
regular testing is performed with OpenMPI 1.10.3+ . Most application
tests are currently only built if MPI is
used.  MPI is enabled in Tangram by setting the CMake variable
`ENABLE_MPI=True`.

On-node parallelism is exposed through
the [Thrust](https://thrust.github.io) library.  Enabling Thrust
within Tangram requires setting at least two CMake variables:
`ENABLE_THRUST=True` and `THRUST_DIR=<path_to_thrust_directory>`.
Additionally, one can specify the Thrust backend to utilize, with the
default being the OpenMP backend
`THRUST_BACKEND=THRUST_DEVICE_SYSTEM_OMP`.  Tangram also supports the
`THRUST_DEVICE_SYSTEM_TBB` backend.  Regular testing happens with
Thrust 1.8. **If you turn
on Thrust for multi-threading-enabled executables, the team strongly
recommends linking to the TCMalloc library available in [Google
Performance Tools](https://github.com/gperftools/gperftools) 
to see the expected scaling.**

## Obtaining Tangram

The latest release of [Tangram](https://github.com/laristra/Tangram)
can be found on GitHub.  Tangram makes use of git submodules, so it must be
cloned recursively:

```sh
git clone --recursive https://github.com/laristra/Tangram
```

## Building

Tangram uses the CMake build system.  In the simplest case where you
want to build a serial version of the code, and CMake knows where to
find your Boost and LAPACKE installations, one can do

```sh
Tangram/ $ mkdir build
Tangram/ $ cd build
Tangram/build/ $ cmake ..
Tangram/build/ $ make
```

This will build a serial version of the code into a library (without
any tests).  A more complete build with MPI, Thrust and TCMalloc (for on-node
parallelism), unit and application test support, documentation
support, and support for both [Jali](https://github.com/lanl/jali) and
[XMOF2D](https://github.com/laristra/XMOF2D) libraries would look
like:

~~~sh
Tangram/ $ mkdir build
Tangram/ $ cd build
Tangram/build/ $ cmake -DENABLE_UNIT_TESTS=True \
                       -DENABLE_MPI=True \
                       -DENABLE_THRUST=True 
                       -DTHRUST_DIR=/path/to/thrust/include/directory \
                       -DENABLE_TCMALLOC=True \
                       -DTCMALLOC_LIB=path/to/TCMalloc/library \
                       -DJali_DIR=path/to/Jali/lib \
                       -DXMOF2D_DIR=path/to/XMOF2D/lib \
                       -DENABLE_DOXYGEN=True \
                       -DLAPACKE_DIR=/path/to/LAPACKE
                       ..
Tangram/build/ $ make           # builds the library and tests
Tangram/build/ $ make test      # runs the tests
Tangram/build/ $ make doxygen   # builds this HTML and a PDF form of the documentation
Tangram/build/ $ make install   # installs the Tangram library and headers into CMAKE_INSTALL_PREFIX
~~~

## Useful CMake Flags
Below is a non-exhaustive list of useful CMake flags for building
Tangram.

| CMake flag:type | Description | Default |
|:----------|:------------|:--------|
| `CMAKE_BUILD_TYPE:STRING`| `Debug` or optimized `Release` build | `Debug` |
| `CMAKE_INSTALL_PREFIX:PATH` | Location for the Tangram library and headers to be installed | `/usr/local` |
| `CMAKE_PREFIX_PATH:PATH` | Locations where CMake can look for packages | "" |
| `ENABLE_APP_TESTS:BOOL` | Turn on compilation and test harness of application tests | `False` |
| `ENABLE_DOXYGEN:BOOL` | Create a target to build this documentation | `False` |
| `ENABLE_MPI:BOOL` | Build with support for MPI | `False` |
| `ENABLE_TCMALLOC:BOOL` | Build with support for TCMalloc | `False` |
| `ENABLE_THRUST:BOOL` | Turn on Thrust support for on-node parallelism | `False` |
| `ENABLE_UNIT_TESTS:BOOL` | Turn on compilation and test harness of unit tests | `False` |
| `LAPACKE_DIR:PATH` | Hint location for CMake to find LAPACKE include and library files | NO_DEFAULT |
| `Jali_DIR:PATH` | Hint location for CMake to find Jali.  This version of Tangram works with version 0.9.8 of Jali | NO_DEFAULT |
| `XMOF2D_DIR:PATH` | Hint location for CMake to find XMOF2D | NO_DEFAULT |
| `TCMALLOC_LIB:PATH` | The TCMalloc library to use | NO_DEFAULT |
| `THRUST_DIR:PATH` | Directory of the Thrust install | NO_DEFAULT |
| `BOOST_ROOT:PATH` | Directory of the Boost install | NO_DEFAULT |
| `THRUST_BACKEND:STRING` | Backend to use for Thrust | `"THRUST_DEVICE_SYSTEM_OMP"` |
