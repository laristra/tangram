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
#include "tangram/support/MatPoly.h"

/// Test the CellMatPoly structure for two-dimensional cells

TEST(CellMatPoly, Mesh2D) {

  Jali::MeshFactory mf(MPI_COMM_WORLD);

  if (Jali::framework_available(Jali::MSTK))
    mf.framework(Jali::MSTK);
  mf.included_entities({Jali::Entity_kind::EDGE, Jali::Entity_kind::FACE});

  // Make a 1-cell two-dimensional mesh
  std::shared_ptr<Jali::Mesh> mesh = mf(0.0, 0.0, 1.0, 1.0, 1, 1);

  // Distance tolerance
  double dst_tol = sqrt(2)*std::numeric_limits<double>::epsilon();
  
  // Make a 2-material CellMatPoly object for the cell

  Tangram::CellMatPoly<2> cellmatpoly;

  // Cell info
  int cellid = 0;
  cellmatpoly.set_cell(cellid);

  std::vector<int> cnodes(4);
  mesh->cell_get_nodes(cellid, &cnodes);

  // Two material polygons 
  //
  //   o----------o
  //   |         /|
  //   |        / |
  //   |       /  |
  //   |      /   |
  //   o-----o----o


  // Specify the first polygon with parent info for points and faces

  std::vector<Tangram::Point<2>> points0(4);
  std::vector<Tangram::Entity_kind> vparentkind0(4);
  std::vector<int> vparentid0(4);

  points0[0] = Tangram::Point<2>(0.0, 0.0);
  vparentkind0[0] = Tangram::Entity_kind::NODE;
  vparentid0[0] = 0;                            // parent is node 0
  points0[1] = Tangram::Point<2>(0.5, 0.0);
  vparentkind0[1] = Tangram::Entity_kind::FACE;
  vparentid0[1] = 0;                            // parent is face 0
  points0[2] = Tangram::Point<2>(1.0, 1.0);
  vparentkind0[2] = Tangram::Entity_kind::NODE;
  vparentid0[2] = 2;                            // parent is node 2
  points0[3] = Tangram::Point<2>(0.0, 1.0);
  vparentkind0[3] = Tangram::Entity_kind::NODE;
  vparentid0[3] = 3;                            // parent is node 3

  std::vector<Tangram::Entity_kind> fparentkind0(4);
  std::vector<int> fparentid0(4);

  fparentkind0[0] = Tangram::Entity_kind::FACE;
  fparentid0[0] = 0;                            // parent is face 0
  fparentkind0[1] = Tangram::Entity_kind::CELL;
  fparentid0[1] = 0;                            // parent is cell 0
  fparentkind0[2] = Tangram::Entity_kind::FACE;
  fparentid0[2] = 2;                            // parent is face 2
  fparentkind0[3] = Tangram::Entity_kind::FACE;
  fparentid0[3] = 3;                            // parent is face 3

  // expected vertices of faces for matpoly 0

  std::vector<std::vector<int>> fverts0;
  fverts0.push_back({0, 1});
  fverts0.push_back({1, 2});
  fverts0.push_back({2, 3});
  fverts0.push_back({3, 0});

  // expected matpolys connected to each face
  std::vector<std::array<int, 2>> fmatpolys0;
  fmatpolys0.push_back({0, -1});
  fmatpolys0.push_back({0, 1});  // 2nd face is internal; cncted to matpolys 0,1
  fmatpolys0.push_back({0, -1});
  fmatpolys0.push_back({0, -1});

  cellmatpoly.add_matpoly(0, -1, 4, &(points0[0]), dst_tol, &(vparentkind0[0]),
                          &(vparentid0[0]), &(fparentkind0[0]),
                          &(fparentid0[0]));


  // Specify the second polygon with parent info only for points.
  // Trigger a particular check (and possibility of failure) by
  // specifying the internal face points first

  std::vector<Tangram::Point<2>> points1(3);
  std::vector<Tangram::Entity_kind> vparentkind1(3);
  std::vector<int> vparentid1(3);

  points1[0] = Tangram::Point<2>(1.0, 1.0);
  vparentkind1[0] = Tangram::Entity_kind::NODE;
  vparentid1[0] = 2;                            // parent is node 2
  points1[1] = Tangram::Point<2>(0.5, 0.0);
  vparentkind1[1] = Tangram::Entity_kind::FACE;
  vparentid1[1] = 0;                            // parent is face 0
  points1[2] = Tangram::Point<2>(1.0, 0.0);
  vparentkind1[2] = Tangram::Entity_kind::NODE;
  vparentid1[2] = 1;                            // parent is node 1

  // expected vertices of faces for matpoly 1

  std::vector<std::vector<int>> fverts1;
  fverts1.push_back({1, 2});  // this face got defined earlier as 1, 2
  fverts1.push_back({1, 4});
  fverts1.push_back({4, 2});

  // expected matpolys connected to each face
  std::vector<std::array<int, 2>> fmatpolys1;
  fmatpolys1.push_back({0, 1});  // 1st face is internal; cncted to matpolys 0,1
  fmatpolys1.push_back({1, -1});
  fmatpolys1.push_back({1, -1});

  cellmatpoly.add_matpoly(1, -1, 3, &(points1[0]), dst_tol, &(vparentkind1[0]),
                          &(vparentid1[0]), nullptr, nullptr);


  // Verify the matpoly info in cell 0

  ASSERT_EQ(cellid, cellmatpoly.cell());

  // There should be two material polygons, with material IDs 0 and 1

  ASSERT_EQ(2, cellmatpoly.num_matpolys());
  ASSERT_EQ(0, cellmatpoly.matpoly_matid(0));
  ASSERT_EQ(1, cellmatpoly.matpoly_matid(1));


  // Verify info for matpoly 0
  {
    std::vector<int> mverts_out = cellmatpoly.matpoly_vertices(0);
    ASSERT_EQ(4, mverts_out.size());
    ASSERT_EQ(0, mverts_out[0]);
    ASSERT_EQ(1, mverts_out[1]);
    ASSERT_EQ(2, mverts_out[2]);
    ASSERT_EQ(3, mverts_out[3]);
    
    for (int i = 0; i < 4; i++) {
      int v = mverts_out[i];
      ASSERT_EQ(vparentkind0[i], cellmatpoly.matvertex_parent_kind(v));
      ASSERT_EQ(vparentid0[i], cellmatpoly.matvertex_parent_id(v));
    }

    std::vector<Tangram::Point<2>> mpoints = cellmatpoly.matpoly_points(0);
    ASSERT_EQ(4, mpoints.size());
    for (int i = 0; i < 4; i++) {
      ASSERT_TRUE(approxEq(points0[i], mpoints[i], 1.0e-8));
    }  
    
    std::vector<int> const& mfaces_out = cellmatpoly.matpoly_faces(0);
    ASSERT_EQ(4, mfaces_out.size());
    ASSERT_EQ(0, mfaces_out[0]);
    ASSERT_EQ(1, mfaces_out[1]);
    ASSERT_EQ(2, mfaces_out[2]);
    ASSERT_EQ(3, mfaces_out[3]);
    
    for (int i = 0; i < 4; i++) {
      int f = mfaces_out[i];
      std::vector<int> const& mfverts_out = cellmatpoly.matface_vertices(f);
      for (int j = 0; j < mfverts_out.size(); j++)
        ASSERT_EQ(fverts0[i][j], mfverts_out[j]);
      
      int fmatpolys_out[2];
      cellmatpoly.matface_matpolys(f, &(fmatpolys_out[0]), &(fmatpolys_out[1]));
      ASSERT_EQ(fmatpolys0[i][0], fmatpolys_out[0]);
      ASSERT_EQ(fmatpolys0[i][1], fmatpolys_out[1]);
      
      if (fmatpolys_out[0] != -1 && fmatpolys_out[1] != -1)
        ASSERT_TRUE(cellmatpoly.matface_is_interface(mfaces_out[i]));
      else
        ASSERT_TRUE(!cellmatpoly.matface_is_interface(i));
      
      ASSERT_EQ(fparentkind0[i], cellmatpoly.matface_parent_kind(f));
      ASSERT_EQ(fparentid0[i], cellmatpoly.matface_parent_id(f));
    }
    
    // Verify volume
    
    ASSERT_NEAR(0.75*mesh->cell_volume(0), cellmatpoly.matpoly_volume(0),
                1.0e-08);
    
    // Verify centroid
    
    Tangram::Point<2> expcen0(0.388888889, 0.55555556);
    Tangram::Point<2> mcen0 = cellmatpoly.matpoly_centroid(0);
    ASSERT_TRUE(approxEq(expcen0, mcen0, 1.0e-08));
  }

  // Verify info for matpoly 1

  {
    std::vector<int> mverts_out = cellmatpoly.matpoly_vertices(1);
    ASSERT_EQ(3, mverts_out.size());
    ASSERT_EQ(2, mverts_out[0]);
    ASSERT_EQ(1, mverts_out[1]);
    ASSERT_EQ(4, mverts_out[2]);
    
    for (int i = 0; i < 3; i++) {
      int v = mverts_out[i];
      ASSERT_EQ(vparentkind1[i], cellmatpoly.matvertex_parent_kind(v));
      ASSERT_EQ(vparentid1[i], cellmatpoly.matvertex_parent_id(v));
    }

    std::vector<Tangram::Point<2>> mpoints = cellmatpoly.matpoly_points(1);
    ASSERT_EQ(3, mpoints.size());
    for (int i = 0; i < 3; i++) {
      ASSERT_TRUE(approxEq(points1[i], mpoints[i], 1.0e-08));
    }  

    std::vector<int> const& mfaces_out = cellmatpoly.matpoly_faces(1);
    ASSERT_EQ(3, mfaces_out.size());
    ASSERT_EQ(1, mfaces_out[0]);
    ASSERT_EQ(4, mfaces_out[1]);
    ASSERT_EQ(5, mfaces_out[2]);

    for (int i = 0; i < 3; i++) {
      int f = mfaces_out[i];
      std::vector<int> const& mfverts_out = cellmatpoly.matface_vertices(f);
      for (int j = 0; j < mfverts_out.size(); j++)
        ASSERT_EQ(fverts1[i][j], mfverts_out[j]);
      
      int fmatpolys_out[2];
      cellmatpoly.matface_matpolys(f, &(fmatpolys_out[0]), &(fmatpolys_out[1]));
      ASSERT_EQ(fmatpolys1[i][0], fmatpolys_out[0]);
      ASSERT_EQ(fmatpolys1[i][1], fmatpolys_out[1]);
      
      if (fmatpolys_out[0] != -1 && fmatpolys_out[1] != -1)
        ASSERT_TRUE(cellmatpoly.matface_is_interface(f));
      else
        ASSERT_TRUE(!cellmatpoly.matface_is_interface(f));
      
      if (f != 1)
        ASSERT_EQ(Tangram::Entity_kind::UNKNOWN_KIND,
                  cellmatpoly.matface_parent_kind(f));
      else
        ASSERT_EQ(Tangram::Entity_kind::CELL,
                  cellmatpoly.matface_parent_kind(f));


      if (f != 1)
        ASSERT_EQ(-1, cellmatpoly.matface_parent_id(f));
      else
        ASSERT_EQ(0, cellmatpoly.matface_parent_id(f));
    }

    ASSERT_NEAR(0.25*mesh->cell_volume(0), cellmatpoly.matpoly_volume(1),
                1.0e-08);
    
    Tangram::Point<2> expcen1(2.5/3.0, 1.0/3.0);
    Tangram::Point<2> mcen1 = cellmatpoly.matpoly_centroid(1);
    ASSERT_TRUE(approxEq(expcen1, mcen1, 1.0e-08));
  }
  
  // Extract matpoly 1 as a MatPoly object
  {
    Tangram::MatPoly<2> MatPoly1 = cellmatpoly.get_ith_matpoly(1);
    
    //Verify material ID
    ASSERT_EQ(1, MatPoly1.mat_id());
    
    //Verify vertices
    const std::vector<Tangram::Point2>& matpoly_points = MatPoly1.points();
    ASSERT_EQ(3, MatPoly1.num_vertices());
    for (int ivrt = 0; ivrt < 3; ivrt++)
      ASSERT_TRUE(approxEq(points1[ivrt], matpoly_points[ivrt], 1.0e-15));
    
    std::vector<std::vector<int>> expected_mp1_faces = {
      {0, 1}, {1, 2}, {2, 0} };
    
    //Verify faces
    ASSERT_EQ(3, MatPoly1.num_faces());
    for (int iface = 0; iface < 3; iface++) {
      const std::vector<int>& face_vertices = MatPoly1.face_vertices(iface);
      ASSERT_EQ(2, face_vertices.size());
      ASSERT_EQ(expected_mp1_faces[iface][0], face_vertices[0]);
      ASSERT_EQ(expected_mp1_faces[iface][1], face_vertices[1]);
    }
  }
}
