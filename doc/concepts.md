# Tangram Concepts      {#concepts}

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
coincident entities in the computational mesh.

## Reconstructor

Tangram works with reconstructors through a wrapper.  All such wrappers should
expose methods for setting the material moments, for specifying the sets of cells
for which `CellMatPoly` objects are to be generated, and the operator to
perform interface reconstruction and create the corresponding `CellMatPoly` objects.
Tangram provides several built-in reconstructors:

- Tangram::SLIC - implements the Simple Linear Interface Calculation algorithm in 2D and 3D, 
recommended for testing only.
- Tangram::VOF - implements the Volume-of-Fluid (Youngs') method in 2D and 3D.
- Tangram::MOF - implements the Moment-of-Fluid method in 2D and 3D.

The standard intersection algorithm used in Tangram is implemented in the
[r3d](https://github.com/devonmpowell/r3d) library, but custom intersectors 
can also be used by implementing operators for clipping and splitting of polytopes. 
Standalone reconstructors can also be used by writing a wrapper, 
which is demonstrated with the wrapper for the XMOF2D library.

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