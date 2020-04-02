# Welcome to Tangram!   {#mainpage}

Tangram is a framework that computational physics applications can use
to build a highly customized, hybrid parallel (MPI+X) interface
reconstruction library for constructing in-cell material polytopes
based on given volume fractions and, optionally, material centroid data.

Tangram aims to provide:
- A modern, modular design - pick and choose your preferred
  intersection and interface reconstruction methods
- Support for general polytopal meshes and higher-order material moments
- High flexibility for application customization
- Algorithms that take advantage of both distributed and on-node parallelism
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
- [CMake](https://cmake.org) 3.13+
- [Wonton](https://github.com/laristra/wonton)

Distributed parallelism of Tangram is currently supported through MPI;
regular testing is performed with OpenMPI 2.1.0+ . Most application
tests are currently only built if MPI is used.
MPI is enabled in Tangram by setting the CMake variable
`ENABLE_MPI=True` **and** using a version of Wonton built with MPI.

On-node parallelism is exposed through
the [Thrust](https://thrust.github.io) library.  Enabling Thrust
within Tangram requires setting at least two CMake variables:
`ENABLE_THRUST=True` **and** using a version of Wonton built with Thrust.

## Obtaining Tangram

The latest release of [Tangram](https://github.com/laristra/tangram)
can be found on GitHub.  Tangram makes use of git submodules, so it must be
cloned recursively:

```sh
git clone --recursive https://github.com/laristra/tangram
```

## Building

Tangram uses the CMake build system.  In the simplest case where you
want to build a serial version of the code, and CMake knows where to
find your Boost and LAPACKE installations, one can do

```sh
Tangram/ $ mkdir build
Tangram/ $ cd build
Tangram/build/ $ cmake -DWONTON_ROOT:FILEPATH=/path/to/wonton ..
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
                       -DENABLE_THRUST=True \
                       -DWONTON_ROOT=path/to/Jali/lib \
					   -DENABLE_XMOF2D=True \
                       -DXMOF2D_ROOT=path/to/XMOF2D/lib \
                       -DENABLE_DOXYGEN=True \
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
| `ENABLE_THRUST:BOOL` | Turn on Thrust support for on-node parallelism | `False` |
| `ENABLE_UNIT_TESTS:BOOL` | Turn on compilation and test harness of unit tests | `False` |
| `WONTON_ROOT:PATH` | Path to Wonton installation under which wontonConfig.cmake may be found | NO_DEFAULT |
| `XMOF2D_ROOT:PATH` | Hint location for CMake to find XMOF2DConfig.cmake | NO_DEFAULT |
