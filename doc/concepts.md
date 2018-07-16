# tangram Concepts      {#concepts}

Interface reconstruction within tangram is based on provided material
moments: the required data is material volume fractions and, for some 
methods (e.g. MOF), material centroids. It is assumed that in a single
cell materials are not repeated, i.e. for every contained material its
moments are given only once.

The results of reconstruction in tangram are stored in `CellMatPoly`
objects, where a single `CellMatPoly` contains material polytopes for
the corresponding cell.  Standard driver performs interface reconstruction
for all multi-material cell, but it can be customized to only operate
on specific sets of cells or to create `CellMatPoly` objects even for
single-material cells.  While `CellMatPoly` stores all the material
polytopes in a cell, a specific material polytope can be added or 
extracted as a `MatPoly` object.

In the simplest case, the reconstructor uses material volume fractions
to create `MatPoly` objects for every material in a cell, then adds
them to the corresponding `CellMatPoly` object.  Note that it is possible
for a single material to be represented by several `MatPoly` objects.
The current implementation of `CellMatPoly` attemps to identify identical
nodes and faces of stores material polytopes and to establish connectivity
between them.

Tangram works with your underlying mesh through a wrapper that provides 
an interface to the queries needed to perform interface reconstruction.  
For an example of the requirements of the wrappers, see the 
[Example Use](@ref example) page.  When we refer to _mesh_ in terms
of the operations, we really mean _mesh wrappers_.

## MatPoly

Light-weight data structure for storing material polytopes in boundary-face
representation.  Provides methods for faceting faces and computing moments
that are consistent with Tangram::AuxMeshTopology.

## CellMatPoly

The current version stores all material polytopes in a cell as a local
mesh, and is able to answer topological queries.  When a material polytope
is added to `CellMatPoly`, its connectivity is established and local 
topology is updated.  Note that for material polytopes to be identified
as adjacent, their faces have to match exactly.  If provided, `CellMatPoly`
can also store parentage data for local mesh entities, i.e. the indices of
coincident entities of the computational mesh.

## Reconstructor

Tangram works with reconstructors through a wrapper.  All such wrappers should
expose methods for setting the material moments, for specifying the sets of cells
for which `CellMatPoly` objects are to be generated, and the operator to
perform interface reconstruction and create the corresponding `CellMatPoly` objects.
Tangram provides several built-in reconstructors, such as SLIC and VOF, which are
templated on the intersection routine.  The standard intersection algorithm for
tangram is r3d, but custom intersectors can also be used.
Standalone reconstructors can also be used by writing a wrapper, which is 
demonstrated with the wrapper for the XMOF2D library.

----

## Driver

Tangram comes with a standard _driver_ that is templated on computational mesh, 
reconstructor and intersection routine, thus allowing high degree of customization.
The standard driver performs interface reconstruction on all multi-material cells
in batches formed based on the number of materials in a cell: if shared memory
parallelism is utilised, this generally improves the load-balancing.

The driver is used within our application tests within the
`apps` directory.  The applications choose a particular mesh and 
select a specfic reconstructor.  Users are encouraged to write their
own specialized drivers, but the above should serve as a starting
point.