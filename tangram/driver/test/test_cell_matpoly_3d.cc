/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#include "tangram/driver/CellMatPoly.h"

#include "gtest/gtest.h"
#include "mpi.h"

#include "Mesh.hh"
#include "MeshFactory.hh"

#include "tangram/support/tangram.h"

/// Test the CellMatPoly structure for three-dimensional cells

TEST(CellMatPoly, Mesh3D) {

  Jali::MeshFactory mf(MPI_COMM_WORLD);

  if (Jali::framework_available(Jali::MSTK))
    mf.framework(Jali::MSTK);
  mf.included_entities({Jali::Entity_kind::EDGE, Jali::Entity_kind::FACE});

  // Make a 1-cell three-dimensional mesh
  std::shared_ptr<Jali::Mesh> mesh = mf(0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1, 1, 1);

  // Distance tolerance
  double dst_tol = sqrt(3)*std::numeric_limits<double>::epsilon();

  // Make a 2-material CellMatPoly object for the cell

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
  vparentid0[2] = 2;                            // parent is node 2
  points0[3] = Tangram::Point<3>(0.0, 0.0, 1.0);
  vparentkind0[3] = Tangram::Entity_kind::NODE;
  vparentid0[3] = 3;                            // parent is node 3
  points0[4] = Tangram::Point<3>(0.0, 1.0, 0.0);
  vparentkind0[4] = Tangram::Entity_kind::NODE;
  vparentid0[4] = 0;                            // parent is node 0
  points0[5] = Tangram::Point<3>(0.5, 1.0, 0.0);
  vparentkind0[5] = Tangram::Entity_kind::EDGE;
  vparentid0[5] = 2;                            // parent is edge 0
  points0[6] = Tangram::Point<3>(1.0, 1.0, 1.0);
  vparentkind0[6] = Tangram::Entity_kind::NODE;
  vparentid0[6] = 6;                            // parent is node 2
  points0[7] = Tangram::Point<3>(0.0, 1.0, 1.0);
  vparentkind0[7] = Tangram::Entity_kind::NODE;
  vparentid0[7] = 3;                            // parent is node 3


  std::vector<Tangram::Entity_kind> fparentkind0(8);
  std::vector<int> fparentid0(8);

  fparentkind0[0] = Tangram::Entity_kind::FACE;
  fparentid0[0] = 0;                            // parent is face 0
  fparentkind0[1] = Tangram::Entity_kind::FACE;
  fparentid0[1] = 1;                            // parent is face 1
  fparentkind0[2] = Tangram::Entity_kind::FACE;
  fparentid0[2] = 2;                            // parent is face 2
  fparentkind0[3] = Tangram::Entity_kind::FACE;
  fparentid0[3] = 3;                            // parent is face 3
  fparentkind0[4] = Tangram::Entity_kind::CELL;
  fparentid0[4] = 0;                            // parent is cell 0
  fparentkind0[5] = Tangram::Entity_kind::FACE;
  fparentid0[5] = 5;                            // parent is face 5

  // expected vertices of faces for matpoly 0

  std::vector<std::vector<int>> fverts0_expected;
  fverts0_expected.push_back({0, 4, 5, 1});
  fverts0_expected.push_back({3, 2, 6, 7});
  fverts0_expected.push_back({0, 3, 7, 4});
  fverts0_expected.push_back({0, 1, 2, 3});
  fverts0_expected.push_back({1, 5, 6, 2});
  fverts0_expected.push_back({5, 4, 7, 6});

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

  // expected matpolys connected to each face
  std::vector<std::array<int, 2>> fmatpolys0;
  fmatpolys0.push_back({0, -1});
  fmatpolys0.push_back({0, -1});  
  fmatpolys0.push_back({0, -1});
  fmatpolys0.push_back({0, -1});
  fmatpolys0.push_back({0, 1});
  fmatpolys0.push_back({0, -1});

  cellmatpoly.add_matpoly(0, -1, 8, &(points0[0]), dst_tol, &(vparentkind0[0]),
                          &(vparentid0[0]), 6, &(nfv0[0]), &(fverts0_flat[0]),
                          &(fparentkind0[0]), &(fparentid0[0]));


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
  vparentid1[2] = 2;                            // parent is node 5
  points1[3] = Tangram::Point<3>(0.5, 1.0, 0.0);
  vparentkind1[3] = Tangram::Entity_kind::EDGE;
  vparentid1[3] = 2;                            // parent is edge 2
  points1[4] = Tangram::Point<3>(1.0, 1.0, 0.0);
  vparentkind1[4] = Tangram::Entity_kind::NODE;
  vparentid1[4] = 2;                            // parent is node 2
  points1[5] = Tangram::Point<3>(1.0, 1.0, 1.0);
  vparentkind1[5] = Tangram::Entity_kind::NODE;
  vparentid1[5] = 6;                            // parent is node 6

  // expected vertices of faces for matpoly 1

  std::vector<std::vector<int>> fverts1_expected;
  fverts1_expected.push_back({1, 8, 2});
  fverts1_expected.push_back({5, 6, 9});
  fverts1_expected.push_back({1, 5, 9, 8});
  fverts1_expected.push_back({8, 9, 6, 2});
  fverts1_expected.push_back({1, 5, 6, 2});

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


  // expected matpolys connected to each face
  std::vector<std::array<int, 2>> fmatpolys1;
  fmatpolys1.push_back({1, -1});
  fmatpolys1.push_back({1, -1});
  fmatpolys1.push_back({1, -1});
  fmatpolys1.push_back({1, -1});
  fmatpolys1.push_back({0, 1});  // 2nd face is internal; cncted to matpolys 0,1

  cellmatpoly.add_matpoly(1, -1, 6, &(points1[0]), dst_tol, &(vparentkind1[0]),
                          &(vparentid1[0]), 5, &(nfv1[0]), &(fverts1_flat[0]),
                          nullptr, nullptr);


  // Verify the matpoly info in cell 0

  ASSERT_EQ(cellid, cellmatpoly.cell());

  // There should be two material polygons, with material IDs 0 and 1

  ASSERT_EQ(2, cellmatpoly.num_matpolys());
  ASSERT_EQ(0, cellmatpoly.matpoly_matid(0));
  ASSERT_EQ(1, cellmatpoly.matpoly_matid(1));

  ASSERT_EQ(10, cellmatpoly.num_matvertices());
  ASSERT_EQ(10, cellmatpoly.num_matfaces());
  

  // Verify info for matpoly 0
  {
    std::vector<int> mverts_expected = {0, 4, 5, 1, 3, 2, 6, 7};
    std::vector<int> mverts_out = cellmatpoly.matpoly_vertices(0);
    ASSERT_EQ(mverts_expected.size(), mverts_out.size());
    for (int i = 0; i < 8; i++)
      ASSERT_EQ(mverts_expected[i], mverts_out[i]);
   
    for (int i = 0; i < 8; i++) {
      int v = mverts_out[i];
      bool found = false;
      for (int j = 0; j < 8; j++) {
        if (approxEq(points0[j], cellmatpoly.matvertex_point(v), 1.0e-16)) {
          found = true;
          ASSERT_EQ(vparentkind0[j], cellmatpoly.matvertex_parent_kind(v));
          ASSERT_EQ(vparentid0[j], cellmatpoly.matvertex_parent_id(v));
          break;
        }
      }
      ASSERT_TRUE(found);
    }

 //   std::vector<Tangram::Point<3>> mpoints = cellmatpoly.matpoly_points(0);
 //   ASSERT_EQ(mpoints.size(), 8);
 //   for (int i = 0; i < 8; i++) {
 //     ASSERT_TRUE(approxEq(points0[i], mpoints[i], 1.0e-8));
 //   }  
    
    std::vector<int> const& mfaces_out = cellmatpoly.matpoly_faces(0);
    ASSERT_EQ(unsigned(6), mfaces_out.size());
    ASSERT_EQ(0, mfaces_out[0]);
    ASSERT_EQ(1, mfaces_out[1]);
    ASSERT_EQ(2, mfaces_out[2]);
    ASSERT_EQ(3, mfaces_out[3]);
    ASSERT_EQ(4, mfaces_out[4]);
    ASSERT_EQ(5, mfaces_out[5]);
    
    for (int i = 0; i < 6; i++) {
      int f = mfaces_out[i];
      std::vector<int> const& mfverts_out = cellmatpoly.matface_vertices(f);
      int const nb_mfverts_out = mfverts_out.size();
      ASSERT_EQ(nfv0[i], nb_mfverts_out);

      for (int j = 0; j < nb_mfverts_out; j++)
        ASSERT_EQ(fverts0_expected[i][j], mfverts_out[j]);
      
      int fmatpolys_out[2];
      cellmatpoly.matface_matpolys(f, &(fmatpolys_out[0]), &(fmatpolys_out[1]));
      ASSERT_EQ(fmatpolys0[i][0], fmatpolys_out[0]);
      ASSERT_EQ(fmatpolys0[i][1], fmatpolys_out[1]);
      
      if (fmatpolys_out[0] != -1 && fmatpolys_out[1] != -1)
        ASSERT_TRUE(cellmatpoly.matface_is_interface(f));
      else
        ASSERT_TRUE(!cellmatpoly.matface_is_interface(f));

      ASSERT_EQ(fparentkind0[i], cellmatpoly.matface_parent_kind(f));
      ASSERT_EQ(fparentid0[i], cellmatpoly.matface_parent_id(f));
    }
    
    // Verify volume
    
    ASSERT_NEAR(0.75*mesh->cell_volume(0), cellmatpoly.matpoly_volume(0),
                1.0e-08);
    
    // Verify centroid
    
    Tangram::Point<3> expcen0(0.38888888, 0.5, 0.55555555);    
    Tangram::Point<3> mcen0 = cellmatpoly.matpoly_centroid(0);
    ASSERT_TRUE(approxEq(expcen0, mcen0, 1.0e-08));
  }

  // Verify info for matpoly 1

  {
    std::vector<int> mverts_expected = {1, 8, 2, 5, 6, 9};
    std::vector<int> mverts_out = cellmatpoly.matpoly_vertices(1);
    ASSERT_EQ(mverts_expected.size(), mverts_out.size());
    for (int i = 0; i < 6; i++)
      ASSERT_EQ(mverts_expected[i], mverts_out[i]);
   
    for (int i = 0; i < 6; i++) {
      int v = mverts_out[i];
      bool found = false;
      for (int j = 0; j < 6; j++) {
        if (approxEq(points1[j], cellmatpoly.matvertex_point(v), 1.0e-16)) {
          found = true;
          ASSERT_EQ(vparentkind1[j], cellmatpoly.matvertex_parent_kind(v));
          ASSERT_EQ(vparentid1[j], cellmatpoly.matvertex_parent_id(v));
          break;
        }
      }
      ASSERT_TRUE(found);
    }

//    std::vector<Tangram::Point<3>> mpoints = cellmatpoly.matpoly_points(1);
//    ASSERT_EQ(mpoints.size(), 6);
//    for (int i = 0; i < 6; i++) {
//      ASSERT_TRUE(approxEq(points1[i], mpoints[i], 1.0e-08));
//    }  

    std::vector<int> const& mfaces_out = cellmatpoly.matpoly_faces(1);
    ASSERT_EQ(unsigned(5), mfaces_out.size());
    ASSERT_EQ(6, mfaces_out[0]);
    ASSERT_EQ(7, mfaces_out[1]);
    ASSERT_EQ(8, mfaces_out[2]);
    ASSERT_EQ(9, mfaces_out[3]);
    ASSERT_EQ(4, mfaces_out[4]);

    for (int i = 0; i < 5; i++) {
      int f = mfaces_out[i];
      std::vector<int> const& mfverts_out = cellmatpoly.matface_vertices(f);
      int const nb_mfverts_out = mfverts_out.size();
      ASSERT_EQ(nfv1[i], nb_mfverts_out);

      for (int j = 0; j < nb_mfverts_out; j++)
        ASSERT_EQ(fverts1_expected[i][j], mfverts_out[j]);
      
      int fmatpolys_out[2];
      cellmatpoly.matface_matpolys(f, &(fmatpolys_out[0]), &(fmatpolys_out[1]));
      ASSERT_EQ(fmatpolys1[i][0], fmatpolys_out[0]);
      ASSERT_EQ(fmatpolys1[i][1], fmatpolys_out[1]);
      
      if (fmatpolys_out[0] != -1 && fmatpolys_out[1] != -1)
        ASSERT_TRUE(cellmatpoly.matface_is_interface(mfaces_out[i]));
      else {
        ASSERT_TRUE(!cellmatpoly.matface_is_interface(i));
      
	ASSERT_EQ(Tangram::Entity_kind::UNKNOWN_KIND,
		  cellmatpoly.matface_parent_kind(f));
	ASSERT_EQ(-1, cellmatpoly.matface_parent_id(f));
      }
    }

    ASSERT_NEAR(0.25*mesh->cell_volume(0), cellmatpoly.matpoly_volume(1),
                1.0e-08);
    
    Tangram::Point<3> expcen1(0.8333333333, 0.5, 0.3333333333);
    Tangram::Point<3> mcen1 = cellmatpoly.matpoly_centroid(1);
    ASSERT_TRUE(approxEq(expcen1, mcen1, 1.0e-08));
  }
  
  // Extract matpoly 1 as a MatPoly object
  {
    Tangram::MatPoly<3> MatPoly1 = cellmatpoly.get_ith_matpoly(1);
    
    //Verify material ID
    ASSERT_EQ(1, MatPoly1.mat_id());
    
    //Verify coordinates
    std::vector<int> exp_vrt_id = {0, 1, 2, 3, 5, 4};
    const std::vector<Tangram::Point3>& matpoly_points = MatPoly1.points();
    ASSERT_EQ(6, MatPoly1.num_vertices());
    for (int ivrt = 0; ivrt < 6; ivrt++)
      ASSERT_TRUE(approxEq(points1[exp_vrt_id[ivrt]], matpoly_points[ivrt], 1.0e-15));
    
    //Verify faces
    std::vector<std::vector<int>> exp_face_vrts = {
      {0, 1, 2}, {3, 4, 5}, {0, 3, 5, 1}, {1, 5, 4, 2}, {2, 4, 3, 0} };
    ASSERT_EQ(5, MatPoly1.num_faces());
    for (int iface = 0; iface < 5; iface++) {
      const std::vector<int>& face_vertices = MatPoly1.face_vertices(iface);
      ASSERT_EQ(exp_face_vrts[iface].size(), face_vertices.size());
      int const nb_elems = exp_face_vrts[iface].size();
      for (int ivrt = 0; ivrt < nb_elems; ivrt++)
        ASSERT_EQ(exp_face_vrts[iface][ivrt], face_vertices[ivrt]);
    }
  }
}
