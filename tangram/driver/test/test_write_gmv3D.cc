/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#include "gtest/gtest.h"

#include "mpi.h"

// tangram includes
#include "tangram/driver/write_to_gmv.h"
#include "tangram/driver/CellMatPoly.h"
#include "tangram/support/tangram.h"

// wonton includes
#include "wonton/mesh/jali/jali_mesh_wrapper.h"

// Jali includes
#include "Mesh.hh"
#include "MeshFactory.hh"

/// Test the CellMatPoly structure for one-dimensional cells

TEST(WriteCellMatPoly, Mesh3D) {

  Jali::MeshFactory mf(MPI_COMM_WORLD);

  if (Jali::framework_available(Jali::MSTK))
    mf.framework(Jali::MSTK);
  mf.included_entities({Jali::Entity_kind::EDGE, Jali::Entity_kind::FACE});

  // Make an 8-cell three-dimensional mesh
  std::shared_ptr<Jali::Mesh> mesh = mf(0.0, 0.0, 0.0, 2.0, 2.0, 2.0, 2, 2, 2);
  Wonton::Jali_Mesh_Wrapper mesh_wrapper(*mesh);

  // Distance tolerance
  double dst_tol = sqrt(3)*std::numeric_limits<double>::epsilon();

  // Make a 2-material CellMatPoly object for cell 0

  Tangram::CellMatPoly<3> cellmatpoly;

  // Cell info
  int cellid = 0;
  cellmatpoly.set_cell(cellid);

  std::vector<int> cnodes(8);
  mesh->cell_get_nodes(cellid, &cnodes);

  // Two material polyhedra (3D extrusion of 2D figure shown here) 
  // So left polyhedron will be a distorted hex and the right polyhedron
  // will be a sideways triangular prism
  //
  //   o----------o
  //   |         /|
  //   |        / |
  //   |       /  |
  //   |      /   |
  //   o-----o----o
  //


  // Specify the first polygon with parent info for points and faces

  std::vector<Tangram::Point<3>> points0(8);
  std::vector<Tangram::Entity_kind> vparentkind0(8);
  std::vector<int> vparentid0(8);

  points0[0] = Tangram::Point<3>(0.0, 0.0, 0.0);
  vparentkind0[0] = Tangram::Entity_kind::NODE;
  vparentid0[0] = 0;                            // parent is node 0
  points0[1] = Tangram::Point<3>(0.5, 0.0, 0.0);
  vparentkind0[1] = Tangram::Entity_kind::EDGE;
  vparentid0[1] = 0;                            // parent is edge 0
  points0[2] = Tangram::Point<3>(1.0, 0.0, 1.0);
  vparentkind0[2] = Tangram::Entity_kind::NODE;
  vparentid0[2] = 10;                            // parent is node 2
  points0[3] = Tangram::Point<3>(0.0, 0.0, 1.0);
  vparentkind0[3] = Tangram::Entity_kind::NODE;
  vparentid0[3] = 9;                            // parent is node 3
  points0[4] = Tangram::Point<3>(0.0, 1.0, 0.0);
  vparentkind0[4] = Tangram::Entity_kind::NODE;
  vparentid0[4] = 3;                            // parent is node 0
  points0[5] = Tangram::Point<3>(0.5, 1.0, 0.0);
  vparentkind0[5] = Tangram::Entity_kind::EDGE;
  vparentid0[5] = 2;                            // parent is edge 2 ????
  points0[6] = Tangram::Point<3>(1.0, 1.0, 1.0);
  vparentkind0[6] = Tangram::Entity_kind::NODE;
  vparentid0[6] = 13;                            // parent is node 2
  points0[7] = Tangram::Point<3>(0.0, 1.0, 1.0);
  vparentkind0[7] = Tangram::Entity_kind::NODE;
  vparentid0[7] = 12;                            // parent is node 3


  // input vertices of faces for matpoly 1 (local ordering -
  // will be identical for first poly, different for subsequent ones)

  std::vector<std::vector<int>> fverts0_in;
  fverts0_in.push_back({0, 4, 5, 1});
  fverts0_in.push_back({3, 2, 6, 7});
  fverts0_in.push_back({0, 3, 7, 4});
  fverts0_in.push_back({0, 1, 2, 3});
  fverts0_in.push_back({1, 5, 6, 2});
  fverts0_in.push_back({5, 4, 7, 6});


  std::vector<int> nfv0;
  for (int i = 0; i < 6; i++)
    nfv0.push_back(fverts0_in[i].size());
  
  std::vector<int> fverts0_flat;
  for (int i = 0; i < 6; i++)
    for (int j = 0; j < nfv0[i]; j++)
       fverts0_flat.push_back(fverts0_in[i][j]);

  cellmatpoly.add_matpoly(0, 8, &(points0[0]), dst_tol, &(vparentkind0[0]),
                          &(vparentid0[0]), 6, &(nfv0[0]), &(fverts0_flat[0]),
                          nullptr, nullptr);


  // Specify the second polygon with parent info only for points

  std::vector<Tangram::Point<3>> points1(6);
  std::vector<Tangram::Entity_kind> vparentkind1(6);
  std::vector<int> vparentid1(6);

  points1[0] = Tangram::Point<3>(0.5, 0.0, 0.0);
  vparentkind1[0] = Tangram::Entity_kind::EDGE;
  vparentid1[0] = 0;                            // parent is edge 0
  points1[1] = Tangram::Point<3>(1.0, 0.0, 0.0);
  vparentkind1[1] = Tangram::Entity_kind::NODE;
  vparentid1[1] = 1;                            // parent is node 1
  points1[2] = Tangram::Point<3>(1.0, 0.0, 1.0);
  vparentkind1[2] = Tangram::Entity_kind::NODE;
  vparentid1[2] = 10;                           // parent is node 10
  points1[3] = Tangram::Point<3>(0.5, 1.0, 0.0);
  vparentkind1[3] = Tangram::Entity_kind::EDGE;
  vparentid1[3] = 2;                            // parent is edge 2 ?????
  points1[4] = Tangram::Point<3>(1.0, 1.0, 0.0);
  vparentkind1[4] = Tangram::Entity_kind::NODE;
  vparentid1[4] = 4;                            // parent is node 2
  points1[5] = Tangram::Point<3>(1.0, 1.0, 1.0);
  vparentkind1[5] = Tangram::Entity_kind::NODE;
  vparentid1[5] = 13;                            // parent is node 6

  // input vertices of faces for matpoly 1 (local ordering -
  // will be identical for first poly, different for subsequent ones)

  std::vector<std::vector<int>> fverts1_in;
  fverts1_in.push_back({0, 1, 2});
  fverts1_in.push_back({3, 5, 4});
  fverts1_in.push_back({0, 3, 4, 1});
  fverts1_in.push_back({1, 4, 5, 2});
  fverts1_in.push_back({0, 2, 5, 3});
  std::vector<int> nfv1;
  for (int i = 0; i < 5; i++)
    nfv1.push_back(fverts1_in[i].size());
  
  std::vector<int> fverts1_flat;
  for (int i = 0; i < 5; i++)
    for (int j = 0; j < nfv1[i]; j++)
       fverts1_flat.push_back(fverts1_in[i][j]);


  cellmatpoly.add_matpoly(1, 6, &(points1[0]), dst_tol, &(vparentkind1[0]),
                          &(vparentid1[0]), 5, &(nfv1[0]), &(fverts1_flat[0]),
                          nullptr, nullptr);


  // Most cells have 1 material
  std::vector<int> cell_num_mats(8, 1);
  cell_num_mats[0] = 2;  // cell 0 has two materials

  std::vector<int> cell_mat_ids({0, 1, 1, 0, 1, 1, 1, 1, 1});

  // Most cells don't have a CellMatPoly structure
  std::vector<std::shared_ptr<Tangram::CellMatPoly<3>>> cellmatpoly_list(8);
  cellmatpoly_list[0] = std::make_shared<Tangram::CellMatPoly<3>>(cellmatpoly);

  Tangram::write_to_gmv<Wonton::Jali_Mesh_Wrapper, 3>(mesh_wrapper, 2,
                                                       cell_num_mats,
                                                       cell_mat_ids,
                                                       cellmatpoly_list,
                                                       "interface_3D.gmv");
}
