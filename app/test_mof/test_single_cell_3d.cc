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

#include "tangram/support/tangram.h"
#include "wonton/mesh/simple/simple_mesh.h"
#include "wonton/mesh/simple/simple_mesh_wrapper.h"

// tangram includes
#include "tangram/driver/driver.h"
#include "tangram/reconstruct/MOF.h"
#include "tangram/driver/write_to_gmv.h"
#include "tangram/utility/get_mat_sym_diff_vol.h"

// wonton includes
#include "wonton/support/Vector.h"


/* Test app for a 3D single-cell two-material mesh.
   Uses SimpleMesh and MOF.
   Reproduces a test case that at some point resulted
   in a numerically singular matrix B being inverted in BFGS */

int main(int argc, char** argv) {
  //Single-cell mesh
  Wonton::Simple_Mesh mesh(0.25, 0.5, 0.25, 
		                       0.30000000000000004, 0.55000000000000004, 0.30000000000000004,
                           1, 1, 1);
  Wonton::Simple_Mesh_Wrapper mesh_wrapper(mesh);

  int ncells = mesh_wrapper.num_owned_cells();

  // Volume and angle tolerances
  double dst_tol = sqrt(3)*std::numeric_limits<double>::epsilon();
  double vol_tol = std::numeric_limits<double>::epsilon();
  std::vector< Tangram::IterativeMethodTolerances_t> ims_tols(2);
  ims_tols[0] = {1000, dst_tol, vol_tol};
  ims_tols[1] = {100, sqrt(2)*std::numeric_limits<double>::epsilon(), dst_tol};

  std::vector<int> mesh_materials = {0, 1};
  int nmesh_materials = static_cast<int>(mesh_materials.size());

  std::vector<int> cell_num_mats(ncells);
  std::vector<int> cell_mat_ids(ncells + 1);
  std::vector<double> cell_mat_volfracs(ncells + 1);
  std::vector<Tangram::Point3> cell_mat_centroids(ncells + 1);

  cell_num_mats[0] = 2;
  cell_mat_ids[0] = 0; cell_mat_ids[1] = 1; 
  cell_mat_volfracs[0] = 0.0027932085949562884; cell_mat_volfracs[1] = 0.99720679140504276;
  cell_mat_centroids[0] = Tangram::Point3(0.25230411145417619, 0.50197006108709041, 0.25717965932730663);
  cell_mat_centroids[1] = Tangram::Point3(0.27505350642975118, 0.52492912983766349, 0.27505535044851798);

  std::vector<int> offsets(ncells, 0);
  for (int icell = 0; icell < ncells - 1; icell++)
    offsets[icell + 1] = offsets[icell] + cell_num_mats[icell];

  int nmmcells = 0;
  std::vector<double> mmcells_material_volumes(nmesh_materials, 0.0);
  for (int icell = 0; icell < ncells; icell++)
    if (cell_num_mats[icell] > 1) {
      double cell_volume = mesh_wrapper.cell_volume(icell);
      for (int icmat = 0; icmat < cell_num_mats[icell]; icmat++) {
        int cur_mat_id = cell_mat_ids[offsets[icell] + icmat];
        int mesh_matid = std::distance(mesh_materials.begin(),
          std::find(mesh_materials.begin(),
                    mesh_materials.end(), cur_mat_id));
        mmcells_material_volumes[mesh_matid] +=
          cell_volume*cell_mat_volfracs[offsets[icell] + icmat];
      }
      nmmcells++;
    }

  std::cout << "Mesh has " << ncells << " cells and " << nmmcells << " multi-material cells" << std::endl;
  for (int imat = 0; imat < nmesh_materials; imat++)
    std::cout << "Material " << mesh_materials[imat] << " has volume of " << mmcells_material_volumes[imat] << std::endl;

  Tangram::Driver<Tangram::MOF, 3, Wonton::Simple_Mesh_Wrapper,
                  Tangram::SplitR3D, Tangram::ClipR3D>
    mof_driver(mesh_wrapper, ims_tols, true);

  mof_driver.set_volume_fractions(cell_num_mats, cell_mat_ids,
                                  cell_mat_volfracs, cell_mat_centroids);
  mof_driver.reconstruct();

#ifndef NDEBUG
  std::vector<std::shared_ptr<Tangram::CellMatPoly<3>>> cellmatpoly_list =
    mof_driver.cell_matpoly_ptrs();

  write_to_gmv(mesh_wrapper, nmesh_materials, cell_num_mats, cell_mat_ids,
               cellmatpoly_list, "single_cell.gmv");  
#endif

  std::cout << "Reconstruction completed!\n";

  return 0;
}
