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
#include "tangram/reconstruct/MOF.h"
#include "tangram/driver/write_to_gmv.h"
#include "tangram/utility/read_material_data.h"

/* Demo app for an unstructured 2D/3D mesh
   and a given material data.
   Uses Jali.
   Reads mesh and material data from file,
   performs interface reconstruction, and outputs material 
   polygons to a gmv file */
/*
template<size_t dim>
using SplitRnD = std::conditional_t<dim==2, Tangram::SplitR2D, Tangram::SplitR3D>;

template<size_t dim>
using ClipRnD = std::conditional_t<dim==2, Tangram::ClipR2D, Tangram::ClipR3D>;
*/

template<size_t dim>
std::vector<std::shared_ptr<Tangram::CellMatPoly<dim>>> 
run_driver(Tangram::Jali_Mesh_Wrapper &mesh_wrapper,
           std::vector<Tangram::IterativeMethodTolerances_t> &ims_tols,
           bool &isconvex,
           std::vector<int> &cell_num_mats,
           std::vector<int> &cell_mat_ids,
           std::vector<double> &cell_mat_volfracs, 
           std::vector<Tangram::Point<dim>>& cell_mat_centroids)
{};

template<>
std::vector<std::shared_ptr<Tangram::CellMatPoly<2>>> 
run_driver<2>(Tangram::Jali_Mesh_Wrapper &mesh_wrapper,
              std::vector<Tangram::IterativeMethodTolerances_t> &ims_tols,
              bool &isconvex,
              std::vector<int> &cell_num_mats,
              std::vector<int> &cell_mat_ids,
              std::vector<double> &cell_mat_volfracs, 
              std::vector<Tangram::Point<2>> &cell_mat_centroids)
{

  Tangram::Driver<Tangram::MOF, 2, Tangram::Jali_Mesh_Wrapper, 
                  Tangram::SplitR2D, Tangram::ClipR2D>
  mof_driver(mesh_wrapper, ims_tols, isconvex);

  mof_driver.set_volume_fractions(cell_num_mats, cell_mat_ids, 
                                  cell_mat_volfracs, cell_mat_centroids);

  mof_driver.reconstruct();    

  std::vector<std::shared_ptr<Tangram::CellMatPoly<2>>>
  cellmatpoly_list = mof_driver.cell_matpoly_ptrs();
  
  return cellmatpoly_list; 
};

template<>
std::vector<std::shared_ptr<Tangram::CellMatPoly<3>>> 
run_driver<3>(Tangram::Jali_Mesh_Wrapper &mesh_wrapper,
              std::vector<Tangram::IterativeMethodTolerances_t> &ims_tols,
              bool &isconvex,
              std::vector<int> &cell_num_mats,
              std::vector<int> &cell_mat_ids,
              std::vector<double> &cell_mat_volfracs, 
              std::vector<Tangram::Point<3>> &cell_mat_centroids)
{

  Tangram::Driver<Tangram::MOF, 3, Tangram::Jali_Mesh_Wrapper, 
                  Tangram::SplitR3D, Tangram::ClipR3D>
  mof_driver(mesh_wrapper, ims_tols, isconvex);

  mof_driver.set_volume_fractions(cell_num_mats, cell_mat_ids, 
                                  cell_mat_volfracs, cell_mat_centroids);

  mof_driver.reconstruct();    

  std::vector<std::shared_ptr<Tangram::CellMatPoly<3>>>
  cellmatpoly_list = mof_driver.cell_matpoly_ptrs();

  return cellmatpoly_list; 
};

template<size_t dim>
void run (std::shared_ptr<Jali::Mesh> inputMesh,
          bool &isconvex,
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
  std::vector< Tangram::IterativeMethodTolerances_t> ims_tols(2) ;
  ims_tols[0].max_num_iter = 1000;
  ims_tols[0].arg_eps = 1.0e-15;
  ims_tols[0].fun_eps = 1.0e-15;
  ims_tols[1].max_num_iter = 100;
  ims_tols[1].arg_eps = 1.0e-13;
  ims_tols[1].fun_eps = 1.0e-13;

  std::vector<std::shared_ptr<Tangram::CellMatPoly<dim>>>
  cellmatpoly_list = run_driver<dim>(mesh_wrapper, ims_tols, isconvex, cell_num_mats,
                     cell_mat_ids, cell_mat_volfracs, cell_mat_centroids);

 //Create MatPoly's for single-material cells
  std::vector<int> offsets(ncells, 0);
  for (int icell = 0; icell < ncells - 1; icell++)
    offsets[icell + 1] = offsets[icell] + cell_num_mats[icell];
  
  for (int icell = 0; icell < ncells; icell++) {
    if (cell_num_mats[icell] == 1) {
      assert(cellmatpoly_list[icell] == nullptr);
      std::shared_ptr< Tangram::CellMatPoly<dim> >
        cmp_ptr(new Tangram::CellMatPoly<dim>(icell));
      Tangram::MatPoly<dim> cell_matpoly;
      cell_get_matpoly(mesh_wrapper, icell, &cell_matpoly);
      cell_matpoly.set_mat_id(cell_mat_ids[offsets[icell]]);
      cmp_ptr->add_matpoly(cell_matpoly);
      cellmatpoly_list[icell] = cmp_ptr;
    }
  }
    
  write_to_gmv(cellmatpoly_list, out_gmv_fname);
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  int world_size = 1;
  MPI_Comm_size(comm, &world_size);
  if (world_size > 1)
    throw std::runtime_error("This app is designed to run in serial!");

  if ((argc < 5) || (argc > 6)) {
      std::ostringstream os;
      os << std::endl <<
      "Correct usage: demo-mof-app <dim=2|3> <isconvex=1|0> <mat_data_filename> " << 
      "<base_mesh_file> <out_gmv_filename>" << std::endl;
      throw std::runtime_error(os.str());
  }
  
  int dim = atoi(argv[1]);
  bool isconvex = atoi(argv[2]);
  std::string in_data_fname = argv[3];
  std::string out_gmv_fname;
  if (argc > 5) out_gmv_fname = argv[5];
  else out_gmv_fname = in_data_fname + ".gmv";

  Jali::MeshFactory mesh_factory(comm);
  mesh_factory.framework(Jali::MSTK);
 
  if (dim == 2)
  {
    mesh_factory.included_entities({Jali::Entity_kind::EDGE});
    std::shared_ptr<Jali::Mesh> mesh = mesh_factory(argv[4]);

    run<2>(mesh, isconvex, in_data_fname, out_gmv_fname);
  }
  else if (dim == 3)
  {
    mesh_factory.included_entities({Jali::Entity_kind::EDGE, Jali::Entity_kind::FACE}); 
    std::shared_ptr<Jali::Mesh> mesh = mesh_factory(argv[4]);
    
    run<3>(mesh, isconvex, in_data_fname, out_gmv_fname);
  }

  MPI_Finalize();

  return 0;
}