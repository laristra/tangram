/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#include <stdlib.h>

#include <memory>
#include <vector>
#include <fstream>
#include <ostream>
#include <iterator>
#include <algorithm>
#include <string>
#include <sstream>
#include <set>
#include <type_traits>

#include "mpi.h"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "tangram/support/tangram.h"
#include "tangram/wrappers/mesh/jali/jali_mesh_wrapper.h"

#include "tangram/intersect/split_r2d.h"
#include "tangram/intersect/split_r3d.h"
#include "tangram/driver/driver.h"
#include "tangram/reconstruct/VOF.h"
#include "tangram/driver/write_to_gmv.h"
#include "tangram/utility/read_material_data.h"

/* Demo app for an unstructured 2D/3D mesh
   and a given material data.
   Uses Jali.
   Reads mesh and material data from file,
   performs interface reconstruction, and outputs material 
   polygons to a gmv file */

template<size_t dim>
using SplitRnD = std::conditional_t<dim==2, Tangram::SplitR2D, Tangram::SplitR3D>;

template<size_t dim>
using ClipRnD = std::conditional_t<dim==2, Tangram::ClipR2D, Tangram::ClipR3D>;

template<size_t dim>
void run (std::shared_ptr<Jali::Mesh> inputMesh,
          std::string in_data_fname,
          std::string out_gmv_fname)
{
  assert(inputMesh != nullptr);
  Tangram::Jali_Mesh_Wrapper mesh_wrapper(*inputMesh, true, false, false);

  int ncells = mesh_wrapper.num_owned_cells();

  std::vector<int> cell_num_mats;
  std::vector<int> cell_mat_ids;
  std::vector<double> cell_mat_volfracs;
  std::vector<Tangram::Point<dim>> cell_mat_centroids;
  read_material_data<Tangram::Jali_Mesh_Wrapper, dim>(mesh_wrapper, in_data_fname, 
    cell_num_mats, cell_mat_ids, cell_mat_volfracs, cell_mat_centroids);
  
  // Volume fraction and angles tolerance
  std::vector< Tangram::IterativeMethodTolerances_t> ims_tols(1) ;
  ims_tols[0]= {.max_num_iter = 1000, .arg_eps = 1.0e-15, .fun_eps = 1.0e-15};


  Tangram::Driver<Tangram::VOF, dim, Tangram::Jali_Mesh_Wrapper, 
                  SplitRnD<dim>, ClipRnD<dim>>
  vof_driver(mesh_wrapper, ims_tols, true);
   

  vof_driver.set_volume_fractions(cell_num_mats, cell_mat_ids, 
                                  cell_mat_volfracs, cell_mat_centroids);
  vof_driver.reconstruct();    
    
  std::vector<std::shared_ptr<Tangram::CellMatPoly<dim>>> 
  cellmatpoly_list = vof_driver.cell_matpoly_ptrs();

  write_to_gmv(cellmatpoly_list, out_gmv_fname);
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  int world_size = 1;
  MPI_Comm_size(comm, &world_size);
  if (world_size > 1)
    throw std::runtime_error("This app is designed to run in serial!");

  if ((argc < 4) || (argc > 5)) {
      std::ostringstream os;
      os << std::endl <<
      "Correct usage: demo-vof-app <dim=2|3> <mat_data_filename> " << 
      "<base_mesh_file> <out_gmv_filename>" << std::endl;
      throw std::runtime_error(os.str());
  }
  
  int dim = atoi(argv[1]);
  std::string in_data_fname = argv[2];
  std::string out_gmv_fname;
  if (argc > 4) out_gmv_fname = argv[4];
  else out_gmv_fname = in_data_fname + ".gmv";

  Jali::MeshFactory mesh_factory(comm);
  mesh_factory.framework(Jali::MSTK);
 
  if (dim == 2)
  {
    mesh_factory.included_entities({Jali::Entity_kind::EDGE});
    std::shared_ptr<Jali::Mesh> mesh = mesh_factory(argv[3]);

    run<2>(mesh, in_data_fname, out_gmv_fname);
  }
  else if (dim == 3)
  {
    mesh_factory.included_entities({Jali::Entity_kind::EDGE, Jali::Entity_kind::FACE}); 
    std::shared_ptr<Jali::Mesh> mesh = mesh_factory(argv[3]);
    
    run<3>(mesh, in_data_fname, out_gmv_fname);
  }

  MPI_Finalize();

  return 0;
}
