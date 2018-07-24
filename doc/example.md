# Simple Mesh Example    {#example}

Tangram provides very crude built-in mesh framework aptly
called `Simple_Mesh`.  The goal of this framework
is to show how one can wrap their favorite mesh infrastructure for
use with tangram - _they should not be used in any production sense._
Tangram also provides implementations of a 3D SLIC and VOF (Youngs') 
interface reconstruction algorithms, as well as a wrapper for XMOF2D,
which is a library implementing the 2D X-MOF method. These serve to
demonstrate how to wrap or implement a specific interface recontruction
algorithm.

----

# The Mesh

## Tangram::Simple_Mesh

This mesh framework is a non-distributed (i.e. has no ghost
information), 3D, regular Cartesian mesh framework.  A Simple_Mesh is
constructed by specifying the extents of the box and the number of
cells in each direction.  The constructor then builds connectivity
information between cell, node, and face indices.  The ordering of
things like the nodes making up a cell, or the faces making up a cell
are specified consistently, but the choice of ordering does not
matter.  There are a few helper functions like
Portage::Simple_Mesh::cell_get_nodes() that will retrieve the indices
of some connected mesh entities given another mesh entity's index.

# The Intersector

The built-in interface reconstruction method allows one to pick different 
intersection routines.  The default option is 
[r3d](https://github.com/devonmpowell/r3d), for which two functors,
ClipR3D and SplitR3D, are implemented.  The former calculates the moments
of a collection of material polyhedra clipped by a plane, and the latter
return two collections of resulting material polyhedra below and above
the plane.  If a custom intersection algorithm is to be used with the
built-in interface reconstruction methods, the respective Clip and Split
functors should be written for it.

# Wrappers

As mentioned on the [Concepts](@ref concepts) page, tangram interfaces
with interface reconstruction algorithms through their _wrappers_.  
The reason for this is that tangram and/or reconstructors may need 
some information that may not be readily available within the original 
mesh frameworks. Moreover, tangram uses a specific material data format 
and stores the results of interface reconstruction in `CellMatPoly`
objects, which should be constructed within a wrapper.

We provide a helper class, Tangram::AuxMeshTopology, that assists in
extending a basic mesh's topological entities needed for some remap
capabilities.  One does not _need_ to use the AuxMeshTopology class,
especially if one's mesh already efficiently supports the advanced
mesh topologies.

We also provide several reconstructor wrappers, which demonstrate how
to implement the interface that is compatible with the standard driver.

## Tangram::AuxMeshTopology

This helper class will build the additional mesh entities and
connectivities from a basic mesh framework wrapper.  The basic mesh
framework wrapper must support cells, nodes, and faces as well as
connectivity and number queries about these entities.

This class is not a complete class design.  In particular, it is
designed to be used within
the
[Curiously Recurring Template Pattern](https://en.m.mwikipedia.org/wiki/Curiously_recurring_template_pattern)(CRTP)
design pattern to achieve static polymorphism.  Under the CRTP in this
case, the basic mesh framework wrapper looks something like

~~~{.cc}
class Basic_Mesh_Wrapper : public AuxMeshTopology<Basic_Mesh_Wrapper> {...};
~~~

In this way, the `Basic_Mesh_Wrapper` can use its own methods, or
defer to AuxMeshTopology to perform more advanced queries.

In addition to the advanced mesh entities (sides, wedges, and
corners), AuxMeshTopology also resolves some advanced connectivity
information.  An example is
Portage::AuxMeshTopology::node_get_cell_adj_nodes(), which, given a node index
in the mesh, returns a vector of all the nodes that are attached to
all cells attached to the given node.  AuxMeshTopology additionally
creates iterators over all the various types of mesh entities,
provides geometric information (e.g. volumes, centroids, etc.) of the
mesh entities, and a method to determine if an entitiy is on the
domain boundary (needed for limiting in higher-order remaps).

_If one does **not** use AuxMeshTopology to help extend their mesh
wrapper's functionality, one should ensure that their mesh wrapper at
least has the same public functions as AuxMeshTopology._ For this
reason, it is advised to utilize AuxMeshTopology where possible, but
to defer to one's own mesh framework wrapper when more efficient
methods are available.

## Tangram::Simple_Mesh_Wrapper

This wraps the Tangram::Simple_Mesh framework.  Simple_Mesh, as its
name suggests, is quite simple and does not know about advanced mesh
entities nor connectivities.  It lets AuxMeshTopology do the vast
majority of the heavy lifting by automatically creating the advanced
mesh features.

Where possible, Simple_Mesh_Wrapper provides quick and efficient
answers to queries that AuxMeshTopology would otherwise solve in a
general sense.  Two trivial examples are:

1. Tangram::Simple_Mesh_Wrapper::cell_get_type(), which determines the
   Portage::Entity_type (e.g. PARALLEL_OWNED, PARALLEL_GHOST, etc.).
   We know Tangram::Simple_Mesh does not know anything about ghost
   information, so we simple return
   Tangram::Entity_type::PARALLEL_OWNED.
2. Tangram::Simple_Mesh_Wrapper::cell_get_element_type(), which
   determines the geometric shape of a given cell from one of the
   Tangram::Element_type's.  A 3D simple mesh is only a structured
   Cartesian mesh, so we always return Tangram::Element_type::HEX.

There are a few other examples within Tangram::Simple_Mesh_Wrapper
where the AuxMeshTopology methods are overwritten to take advantage of
information that is cached within the Tangram::Simple_Mesh.  This is a
prime example of how the CRTP and AuxMeshTopology are intended to be
used.  

## Tangram::XMOF2D_Wrapper

This wraps the XMOF2D reconstructor, which is the implementation of the
2D X-MOF method. XMOF2D uses its own data structures to represent the mesh
and material data, as well as the results of interface reconstruction. 
XMOF2D_Wrapper demonstrates how to interface such custom reconstructor
with the standard tangram driver.

# Applications and Tests

There are several applications that allow to perform interface reconstruction
using the provided reconstructor and mesh wrappers.  In particular, the
`app/xmof-linetest-simplemesh-app/xmof-linetest-simplemesh-app.cc` program 
shows how to wrap mesh and reconstructor objects, set the material data, 
and utilize Tangram::Driver with SimpleMesh and XMOF2D to perform interface 
reconstruction.  This test program expects a single linear material interface
and computes the Hausdorff distance between the reference and reconstructed 
interfaces.  Material data obtained by sampling is provided for several mesh
resolutions. The Tangram::Driver is templated on mesh and reconstructor wrapper type, 
as well as the intersection routine, and can be used with other frameworks.  
It need not be used at all, but is a nice starting point for writing one's 
own interface reconstruction application.