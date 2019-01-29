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
[Wonton](https://github.com/laristra/wonton) documentation.

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
returns two collections of resulting material polyhedra below and above
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
XMOF2D_Wrapper demonstrates how to interface a custom reconstructor
with the standard Tangram driver.

# Putting it all together

Let us illustrate the usage of Tangram with a few lines of code.
First, we create a cubic 10x10x10 mesh in a unitary domain:

```c++
nx = 10; ny = 10; nz = 10;
Wonton::Simple_Mesh mesh(0.0, 0.0, 0.0,
                         1.0, 1.0, 1.0,
                         nx, ny, nz);
```

Then, we wrap it with the Tangram-compatible `Wonton::Simple_Mesh_Wrapper`:

```c++
Wonton::Simple_Mesh_Wrapper mesh_wrapper(mesh);
```

Now, we choose the tolerances for the interface reconstruction method. All the
supplied reconstructors use an iterative procedure to find the position of a linear
material interface with the known orientation that matches the given volume fraction up
to a certain tolerance. Here, we provide that tolerance as the negligible discrepancy
in the objective function, `fun_eps`, and also specify the maximum allowed number of iterations,
`max_num_iter` and the change in argument between two steps of the procedure
that should result in an early termination, `arg_eps`. 

For the `SLIC` and `VOF` reconstructors no other tolerances are required. The `MOF`
reconstructor solves an optimization problem to find the orientation of the linear interface.
Its argument is the normal to the interface, and its objective function evaluates the distance
between the provided material centroid and the centroid of the material polytope corresponding
to that choice of the normal. Therefore, for this iterative method `fun_eps` will correspond
to the tolerance on distance between centroids, and `arg_eps` specifies the norm of the change 
in normal that results in early termination.

It can be seen that the first iterative method deals with volumes or zero-order moments of
material polytopes, while the second deals with centroids or the first-order moments. The way 
tolerances are specified in Tangram reflects that fact:

```c++
std::vector<Tangram::IterativeMethodTolerances_t> ims_tols(2);
ims_tols[0] = {.max_num_iter = 1000, .arg_eps = 1.0e-15, .fun_eps = 1.0e-14};
ims_tols[1] = {.max_num_iter = 100, .arg_eps = 1.0e-12, .fun_eps = 1.0e-14};
```

With tolerances initialized, we can now create a reconstructor object. Here, we choose
the 3D MOF reconstructor and create a driver object, which application use to perform
The interface reconstruction:

```c++
bool all_cells_are_convex = true;
Tangram::Driver<Tangram::MOF, 3, Wonton::Simple_Mesh_Wrapper,
                Tangram::SplitR3D, Tangram::ClipR3D>
  mof_driver(mesh_wrapper, ims_tols, all_cells_are_convex);
```

Here, we specified the chosen reconstructor, the problem dimension, the wrapper type, and 
the standard intersector functors as the template parameters. Because mesh cells are cubes
with planar faces, there is no need to decompose them into tetrahedra, so we set the flag
indicating that all mesh cells are convex.

The last required piece of information is the material data. Let us assume that the vectors
below are initialized by the application:

```c++
std::vector<int> cell_num_mats; //Number of materials in each mesh cell
std::vector<int> cell_mat_ids;  //Indices of contained material for all mesh cells, flattened
std::vector<double> cell_mat_volfracs; //Volume fractions of contained material for all mesh cells, flattened
std::vector<Wonton::Point3> cell_mat_centroids; //Centroids of contained material for all mesh cells, flattened
```

After the material data vectors are populated, we pass this information to the driver/reconstructor:

```c++
mof_driver.set_volume_fractions(cell_num_mats, cell_mat_ids,
                                cell_mat_volfracs, cell_mat_centroids);
```

Finally, we can perform interface reconstruction. By default, it will be performed for all the
Multi-material cells in the mesh:

```c++
mof_driver.reconstruct();
```

The results of the interface reconstruction are stored in ``CellMatPoly` objects:

```c++
std::vector<std::shared_ptr< Tangram::CellMatPoly<3> >> cellmatpoly_list =
  mof_driver.cell_matpoly_ptrs();
```

Note that by default ``CellMatPoly` objects are not created for single-material cells, so
the corresponding entries of the `cellmatpoly_list` vector will contain `nullptr`'s.

# Applications and Tests

There are several demo applications and tests that perform interface reconstruction
using different provided reconstructors.  For example, the `app/demo-mof-app/demo-mof-app.cc`
program shows how to wrap the mesh object, set the material data, 
and utilize Tangram::Driver with Jali and built-in MOF to perform interface 
reconstruction in 2D or 3D on an unstructured mesh given in Exodus II format.
The Tangram::Driver is templated on mesh and reconstructor wrapper type, 
as well as the intersection routine, and can be used with other frameworks.  
It need not be used at all, but is a nice starting point for writing one's 
own interface reconstruction application.
For every implemented reconstructor, we also include app tests. If only the
SimpleMesh mesh framework is available, tests are run for structured grids
(including the case of partitioning cells into simplices), and when Tangram is
built with Jali, reconstructors are also tested on Voronoi meshes. In those tests,
the domain contains three materials with planar/linear interfaces between then,
which form a so-called T-junction. Reconstructors are tested in 2D, or 3D, or both,
depending on dimensionality they support.