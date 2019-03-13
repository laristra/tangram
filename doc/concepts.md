# Tangram Concepts      {#concepts}

In some multi-material, multi-physics problems the computational mesh may not
exactly follow the material boundaries, which results in some of mesh cells 
containing several different materials. This situation is typical for problems
where the Arbitrary Lagrangian-Eulerian (ALE) method is applied. When exact 
conservation is critical, interface reconstruction methods are usually the method 
of choice to compute the material interfaces inside the multi-material cells. These
methods take volume fractions of materials in a cell as the input data and partition
the cell into non-overlapping single-material regions that completely cover the cell.
For every material, the volume of corresponding regions inside the cell matches the
volume given by the provided volume fraction to a user-specified tolerance.

Tangram is an interface reconstruction framework that provides a uniform interface
to implementations of interface reconstruction methods. It features several built-in
reconstructors and supports the addition of other reconstructors through the use of
wrappers. Tangram supports on-node and distributed parallelism and is therefore suitable 
for enabling interface reconstruction in high-performance simulations.

Tangram is templated on both meshes and reconstructors and interacts with them through
the wrappers. The key component of Tangram is the _driver_, which takes material volume 
fractions and optionally centroids, and outputs material polytopes that represent the 
subdivision of a cell. The application only needs to interact with the driver by
specifying the mesh and the reconstructor to be used, providing the material data, and
receiving reconstruction results. Custom drivers can be written and a standard driver
is provided.

Interface reconstruction within Tangram is based on provided material
moments: the required data is material volume fractions and, for some 
methods (e.g. MOF), material centroids. It is assumed that in a single
cell materials are not repeated, i.e. for every contained material its
moments are given only once. Tangram requires material data to be provided
for all mesh cells, including single-material cells.

The results of reconstruction in Tangram are stored in `CellMatPoly`
objects, where a single `CellMatPoly` contains material polytopes for
the corresponding cell.  The standard driver performs interface reconstruction
for all multi-material cells, but custom drivers may be written to only operate
on specific sets of cells or to create `CellMatPoly` objects even for single 
material cells. While `CellMatPoly` stores all the material
polytopes in a cell, a specific material polytope can be added or 
extracted as a `MatPoly` object.

In the simplest case, the reconstructor uses material volume fractions
to create `MatPoly` objects for every material in a cell, then adds
them to the corresponding `CellMatPoly` object.  Note that it is possible
for a single material to be represented by several `MatPoly` objects:
this is typical for non-convex cells that are partitioned into simplices
on a pre-processing step.
The current implementation of `CellMatPoly` attempts to identify identical
nodes and faces of stored material polytopes and to establish connectivity
between them. Depending on the interface reconstruction method, the resulting
topological data might be incorrect.

Tangram works with your underlying mesh through a wrapper that provides 
an interface to the queries needed to perform interface reconstruction.  
For an example of the requirements of the wrappers, see the 
[Example Use](@ref example) page.  When we refer to _mesh_ in terms
of the operations, we really mean _mesh wrappers_. More details about the
requirements of mesh wrappers are given in the documentation
of the support package, [Wonton](https://github.com/laristra/wonton).

## MatPoly

Light-weight data structure for storing material polytopes in boundary-face
representation.  Provides methods for faceting faces, decomposing star-convex 
polytopes into simplices, and computing moments
that are consistent with Wonton::AuxMeshTopology. Pre-computed moments can also
be assigned and stored, in which case they will not be recomputed when queried.
Note that faceting procedure involves the following steps:

- Using the geometric center of the face to partition it into triangles;
- Finding the centroid of the resulting union of triangles;
- Using that centroid to obtain the final triangulation.

Consequently, we require faces have to be star-convex with respect to both 
the geometric center and the centroid.
Similarly, finding the centroid of the polytope involves:

- Faceting faces that are not simplices;
- Finding the centroid of the resulting polytope;
- Using that centroid to partition the polytope into simplices. 

## CellMatPoly

The current version stores all material polytopes in a cell as a local
mesh, and is able to answer topological queries.  When a material polytope
is added to `CellMatPoly`, its connectivity is established and the local 
topology is updated.  Note that for material polytopes to be identified
as adjacent, their shared faces have to match exactly, i.e. the number and the coordinates of 
corresponding face vertices should be exactly the same. If provided, `CellMatPoly`
can also store parentage data for local mesh entities, i.e. the indices of
entities in the computational mesh which the local mesh entities were derived from.

## Reconstructor

Tangram works with reconstructors through a wrapper.  All such wrappers should
expose methods for setting the material moments, for specifying the sets of cells
for which `CellMatPoly` objects are to be generated, and the operator to
perform interface reconstruction and create the corresponding `CellMatPoly` objects.
Tangram provides several built-in reconstructors:

- Tangram::SLIC - implements the Simple Linear Interface Calculation algorithm in 2D and 3D [1], 
recommended for testing only.
- Tangram::VOF - implements the Volume-of-Fluid (Youngs') method in 2D and 3D [2,3];
requires _only_ material volume fractions, but is first order accurate and needs 
material data from the neighboring cells.
- Tangram::MOF - implements the Moment-of-Fluid method in 2D and 3D [4,5]; requires 
material volume fractions _and_ centroids, but is second order accurate. Due to being a
completely local method, requires no information about the neighboring cells and is 
therefore more amenable to parallel implementations.

The standard intersection algorithm used in Tangram is implemented in the
[r3d](https://github.com/devonmpowell/r3d) library, but custom intersectors 
can also be used by implementing operators for clipping and splitting of polytopes. 
Standalone reconstructors can also be used by writing a wrapper, which is demonstrated 
with the wrapper for the XMOF2D library. Note that XMOF2D is originally a serial 
interface reconstructor, but by wrapping it in Tangram we can use it with on-node and
distributed parallelism, thus also enabling its use for distributed multi-material remap
in [Portage](https://github.com/laristra/portage). Similar results can be achieved with
other serial reconstructors.

## Intersector

Built-in reconstructors require an implementation of two functors, Split and Clip.
The included functors utilize the r3d intersector and are named respectively:

- Tangram::SplitR2D - for splitting a convex or non-convex polygon with a line.
- Tangram::SplitR3D - for splitting a convex or non-convex polyhedron with a plane.
- Tangram::ClipR2D - for computing the moments of the part of a convex or non-convex 
  polygon below a line.
- Tangram::ClipR3D - for computing the moments of the part of a convex or non-convex 
  polyhedron below a plane.

These functors are template parameters of the built-in reconstructors and can be replaced
by custom implementations. Note that Split and Clip functors are generally independent of each other
and can be based on different intersectors as long as the respective moments are consistent.
In general, the Clip functor is invoked much more often and is more critical for performance.

----

## Driver

Tangram comes with a standard _driver_ that is templated on computational mesh, 
reconstructor and intersection routine, thus allowing high degree of customization.
The standard driver performs interface reconstruction on all multi-material cells
in batches based on the number of materials in a cell: if shared memory
parallelism is utilized, this generally improves the load-balancing.

The driver is used within our application tests in the
`apps` directory.  The applications choose a particular mesh and 
select a specific reconstructor.  Users are encouraged to write their
own specialized drivers, but the standard driver should serve as a starting
point.

Every driver should provide to the reconstructor the material data for all mesh 
cells on the mesh partition, which is done using the method `set_volume_fractions`
exposed by the reconstructor, and select the subset of cells for which reconstruction
is to be performed using the method `set_cell_indices_to_operate_on`. The latter
can be used to limit the cells for which `CellMatPoly` objects are created or to perform
reconstruction in batches when mesh cells are grouped based on certain cretiria (e.g.
the number of materials). The `set_volume_fractions` requires the data on the number of
materials in every cell, their material indices, volume fractions, and, optionally, 
centroids. This data should be provided for all mesh cells, whether owned or ghosts,
single or multi-material, in the form of flat vectors of the respective type. The
standard driver uses its internal logic to issue batches of cells to the reconstructor
and provides a method for the application to setup the material data in the format
required by the reconstructor. Custom drivers can be written to expose either less or
more of this machinery to their target applications.

----

## BFGS

Tangram includes an implementation of the Broyden–Fletcher–Goldfarb–Shanno algorithm
as described in ["Numerical Optimization"](https://www.springer.com/us/book/9780387303031) 
by J. Nocedal and S. Wright. This optimization algorithm is used in the implementation of 
the MOF interface reconstruction method both in 2D and 3D. 
Note that XMOF2D uses the golden-section search instead. Two additional variants of 
the BFGS algorithm (MWWP and D-BFGS) are also included and can be used if 
convergence issues are observed.

----

## Utilities

Tangram comes with several support functions that can assist in testing and verification:
they can be found in the `utilities` subdirectory. In particular, `get_material_moments` can
be used to generate material data corresponding to a sequence of planar interfaces, and
`rpgtools` would help with generation of material data for complex material distributions. 
Material data generation tools use the r3d intersector and provide both the material moments
and reference polytopes, which allows for measuring the symmetric difference between the
reconstructed and reference material polytopes using the `get_mat_sym_diff_vol` function. 

[1] Noh W.F., Woodward P. "SLIC (Simple Line Interface Calculation)". In: van de Vooren
  A.I., Zandbergen P.J. (eds) Proceedings of the Fifth International Conference on
  Numerical Methods in Fluid Dynamics June 28 – July 2, 1976 Twente University, 
  Enschede. Lecture Notes in Physics, vol. 59, pages 330-340, Springer, Berlin, 
  Heidelberg, 1976.

[2] Youngs D.L., "Time-Dependent Multi-Material Flow with Large Fluid Distortion".
  In Morton K.W. and Baines M.J., (eds) Numerical Methods for Fluid Dynamics, pages
  273-285. Academic Press, 1982.
  
[3] Youngs D.L., "An Interface Tracking Method for a 3D Eulerian Hydrodynamics Code".
  Technical Report 44/92/35, AWRE, 1984.

[4] Dyadechko V., Shashkov M., "Reconstruction of multi-material interfaces from moment 
  data". Journal of Computational Physics, vol. 227, issue 11, pages 5361-5384, 2008.

[5] Ahn H.T., Shashkov M., "Multi-material interface reconstruction on generalized
  polyhedral meshes". Journal of Computational Physics, vol. 226, issue 2, pages 
  2096-2132, 2007.


