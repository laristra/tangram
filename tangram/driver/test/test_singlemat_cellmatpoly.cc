/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#include "wonton/mesh/simple/simple_mesh.h"
#include "wonton/mesh/simple/simple_mesh_wrapper.h"

#include "tangram/support/tangram.h"
#include "tangram/driver/driver.h"
#include "tangram/reconstruct/VOF.h"
#include "tangram/intersect/split_rNd.h"

#include "gtest/gtest.h"

// Test the CellMatPoly creation for cells 
// with material volumes below tolerance

TEST(CellMatPoly, SingleMat) {
  // Make a 2x2 SimpleMesh
  Wonton::Simple_Mesh mesh(0.0, 0.0,
                           1.0, 1.0,
                           2, 2);
  Wonton::Simple_Mesh_Wrapper mesh_wrapper(mesh);

  // Distance tolerance
  double dst_tol = sqrt(2)*std::numeric_limits<double>::epsilon();
  // Volume Tolerance
  double vol_tol = 2.0e-16;
  // Volume fraction tolerance
  double vfrac_tol = 4*vol_tol;

  std::vector< Tangram::IterativeMethodTolerances_t> ims_tols(1);
  ims_tols[0] = {1000, dst_tol, vol_tol}; 

  std::vector<int> cell_num_mats = {1, 2, 3, 4};
  std::vector<int> cell_mat_ids = {0, 
                                   0, 1, 
                                   0, 1, 2, 
                                   0, 1, 2, 3};
  std::vector<double> cell_mat_volfracs = 
    {1.0,
     0.95*vfrac_tol, 1.0 - 0.05*vfrac_tol,
     1.5*vfrac_tol, 0.5*vfrac_tol, 1.0 - 2*vfrac_tol,
     0.5*vfrac_tol, 0.75*vfrac_tol, 0.25*vfrac_tol, 1.0 - 1.5*vfrac_tol};

  // Create a driver for VOF
  Tangram::Driver<Tangram::VOF, 2, Wonton::Simple_Mesh_Wrapper,
                  Tangram::SplitRnD<2>, Tangram::ClipRnD<2>>
    vof_driver(mesh_wrapper, ims_tols, false);

  // Perform reconstruction for the given data
  vof_driver.set_volume_fractions(cell_num_mats, cell_mat_ids,
                                  cell_mat_volfracs);
  vof_driver.reconstruct();

  std::vector<std::shared_ptr<Tangram::CellMatPoly<2>>> cellmatpoly_list =
    vof_driver.cell_matpoly_ptrs();

  // No CellMatPoly should be created for the single-material cell
  ASSERT_EQ(cellmatpoly_list[0], nullptr);
  // CellMatPoly's should be created for all cells with multiple materials
  // as per input data
  for (int icell = 0; icell < 3; icell++)
    ASSERT_NE(cellmatpoly_list[icell + 1], nullptr);

  // Verify materials in cell 1
  const std::vector<int>& cell1_mat_ids = cellmatpoly_list[1]->cell_matids();
  ASSERT_EQ(cell1_mat_ids.size(), 1);
  ASSERT_EQ(cell1_mat_ids[0], 1);
  double cell_volume = mesh_wrapper.cell_volume(1);
  ASSERT_NEAR(cellmatpoly_list[1]->material_moments(1)[0], cell_volume, vol_tol);

  // Verify materials in cell 2
  const std::vector<int>& cell2_mat_ids = cellmatpoly_list[2]->cell_matids();
  ASSERT_EQ(cell2_mat_ids.size(), 2);
  ASSERT_EQ(cell2_mat_ids[0], 0);
  ASSERT_EQ(cell2_mat_ids[1], 2);
  cell_volume = mesh_wrapper.cell_volume(2);
  ASSERT_NEAR(cellmatpoly_list[2]->material_moments(0)[0], 
              cell_volume*1.5*vfrac_tol/(1.0 - 0.5*vfrac_tol), vol_tol);
  ASSERT_NEAR(cellmatpoly_list[2]->material_moments(2)[0], 
              cell_volume*(1.0 - 2*vfrac_tol)/(1.0 - 0.5*vfrac_tol), vol_tol);

  // Verify materials in cell 3
  const std::vector<int>& cell3_mat_ids = cellmatpoly_list[3]->cell_matids();
  ASSERT_EQ(cell3_mat_ids.size(), 1);
  ASSERT_EQ(cell3_mat_ids[0], 3);
  cell_volume = mesh_wrapper.cell_volume(3);
  ASSERT_NEAR(cellmatpoly_list[3]->material_moments(3)[0], cell_volume, vol_tol);
}
