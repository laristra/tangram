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

#include "mpi.h"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "tangram/support/tangram.h"
#include "tangram/wrappers/mesh/jali/jali_mesh_wrapper.h"

#include "tangram/intersect/split_r3d.h"
#include "tangram/driver/driver.h"
#include "tangram/reconstruct/MOF.h"
#include "tangram/driver/write_to_gmv.h"
#include "tangram/utility/read_material_data.h"

/* Demo app for an unstructured 3D mesh
   and a given material data.
   Uses Jali.
   Reads mesh and material data from file,
   performs interface reconstruction, and outputs material 
   polygons to a gmv file */


#include <set>

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  int world_size = 1;
  MPI_Comm_size(comm, &world_size);
  if (world_size > 1)
    throw std::runtime_error("This app is designed to run in serial!");

  if ((argc < 3) || (argc > 4)) {
      std::ostringstream os;
      os << std::endl <<
      "Correct usage: demo-mof-app <mat_data_filename> " << 
      "<base_mesh_file> <out_gmv_filename>" << std::endl;
      throw std::runtime_error(os.str());
  }

  std::string in_data_fname = argv[1];
  std::string out_gmv_fname;
  if (argc > 3) out_gmv_fname = argv[3];
  else out_gmv_fname = in_data_fname + ".gmv";

  Jali::MeshFactory mesh_factory(comm);
  mesh_factory.framework(Jali::MSTK);
  mesh_factory.included_entities({Jali::Entity_kind::EDGE, Jali::Entity_kind::FACE});
  std::shared_ptr<Jali::Mesh> mesh = mesh_factory(argv[2]);
  
  assert(mesh != nullptr);
  Tangram::Jali_Mesh_Wrapper mesh_wrapper(*mesh, true, false, false);

  int ncells = mesh_wrapper.num_owned_cells();

  std::vector<int> cell_num_mats;
  std::vector<int> cell_mat_ids;
  std::vector<double> cell_mat_volfracs;
  std::vector<Tangram::Point3> cell_mat_centroids;
  read_material_data<Tangram::Jali_Mesh_Wrapper, 3>(mesh_wrapper, in_data_fname, 
    cell_num_mats, cell_mat_ids, cell_mat_volfracs, cell_mat_centroids);

  // Volume and angles tolerance
  std::vector<Tangram::IterativeMethodTolerances_t> ims_tols(2);
  ims_tols[0] = {.max_num_iter = 1000, .arg_eps = 1.0e-15, .fun_eps = 1.0e-15};
  ims_tols[1] = {.max_num_iter = 100, .arg_eps = 1.0e-13, .fun_eps = 1.0e-13};

  // Build the driver
  Tangram::Driver<Tangram::MOF, 3, Tangram::Jali_Mesh_Wrapper, 
                  Tangram::SplitR3D, Tangram::ClipR3D> 
    mof_driver(mesh_wrapper, ims_tols, true);

  mof_driver.set_volume_fractions(cell_num_mats, cell_mat_ids, 
                                  cell_mat_volfracs, cell_mat_centroids);
  mof_driver.reconstruct();    

  std::vector<std::shared_ptr<Tangram::CellMatPoly<3>>> cellmatpoly_list = 
    mof_driver.cell_matpoly_ptrs();

  write_to_gmv(cellmatpoly_list, out_gmv_fname);

  MPI_Finalize();

  return 0;
}
