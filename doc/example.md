# Simple Mesh Example    {#example}

Tangram provides a very crude mesh framework aptly
called `Simple_Mesh` through its support package, Wonton.  
The goal of this framework
is to show how one can wrap their favorite mesh infrastructure for
use with Tangram - _they should not be used in any production sense._
The details on the `Wonton::Simple_Mesh` framework, its wrapper 
`Wonton::Simple_Mesh_Wrapper`, and the helper class
`Wonton::AuxMeshTopology`, which builds the additional mesh entities and
connectivities from a basic mesh framework wrapper, can be found in the
Wonton documentation.

Tangram also provides implementations of the 2D and 3D VOF (Youngs') and MOF
interface reconstruction algorithms, as well as a wrapper for XMOF2D,
which is a library implementing the 2D X-MOF method. These serve to
demonstrate how to wrap or implement a specific interface reconstruction
algorithm.

# The Intersector

The built-in interface reconstruction method allows one to pick different 
intersection routines.  The default option is 
[r3d](https://github.com/devonmpowell/r3d), for which two functors,
ClipR3D and SplitR3D, are implemented.  The former calculates the moments
of a collection of material polyhedra below a plane, and the latter
return two collections of resulting material polyhedra below and above
the plane.  If a custom intersection algorithm is to be used with the
built-in interface reconstruction methods, the respective Clip and Split
functors should be written for it. Similarly, ClipR2D and SplitR2D functors
are implemented for 2D problems by using the r2d intersector, which is a part 
of the r3d package.

# Wrappers

As mentioned on the [Concepts](@ref concepts) page, Tangram interfaces
with interface reconstruction algorithms through their _wrappers_.  
The reason for this is that Tangram and/or reconstructors may need 
some information that may not be readily available within the original 
mesh frameworks. Moreover, Tangram uses a specific material data format 
and stores the results of interface reconstruction in `CellMatPoly`
objects, which should be constructed within a wrapper.

Tangram's support package, Wonton, provides a helper class, 
Wonton::AuxMeshTopology, that assists in
extending a basic mesh's topological entities needed for some remap
capabilities.  One does not _need_ to use the AuxMeshTopology class,
especially if one's mesh already efficiently supports the advanced
mesh topologies.

We also provide several reconstructor wrappers, which demonstrate how
to implement the interface that is compatible with the standard driver.


## Tangram::XMOF2D_Wrapper

This wraps the XMOF2D reconstructor, which is the implementation of the
2D X-MOF method. XMOF2D uses its own data structures to represent the mesh
and material data, as well as the results of interface reconstruction. 
XMOF2D_Wrapper demonstrates how to interface such custom reconstructor
with the standard Tangram driver.

# Applications and Tests

There are several demo applications and tests that allow to perform interface reconstruction
using different provided reconstructors.  For example, the `app/demo-mof-app/demo-mof-app.cc`
program shows how to wrap the mesh object, set the material data, 
and utilize Tangram::Driver with Jali and built-in MOF to perform interface 
reconstruction in 2D or 3D on an unstructured mesh given in Exodus II format.
The Tangram::Driver is templated on mesh and reconstructor wrapper type, 
as well as the intersection routine, and can be used with other frameworks.  
It need not be used at all, but is a nice starting point for writing one's 
own interface reconstruction application.
For every implemented reconstructor, we also include app tests. If only
SimpleMesh mesh framework is available, they are run for structured grids
(including the case of partitioning cells into simplices), and when Tangram is
built with Jali, reconstructors are also tested on Voronoi meshes. In those tests,
the domain contains three materials with planar/linear interfaces between then,
which form so-called T-junction. Reconstructors are tested in 2D, or 3D, or both,
depending on dimensionality they support.