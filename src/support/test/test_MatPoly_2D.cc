/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#include "tangram/support/MatPoly.h"

#include "gtest/gtest.h"
#include "tangram/support/tangram.h"

/// Test the MatPoly structure for a 2D polygon

TEST(MatPoly, Mesh2D) {
  int mat_id = 1;
  
  //Test for a unit square
  std::vector<Tangram::Point2> square_points = {
    Tangram::Point2(0.0, 0.0), Tangram::Point2(1.0, 0.0),
    Tangram::Point2(1.0, 1.0), Tangram::Point2(0.0, 1.0) };
  std::vector< std::vector<int> > square_faces = {
    {0, 1}, {1, 2}, {2, 3}, {3, 0} };
  std::vector<Tangram::Point2> face_centroids = {
    Tangram::Point2(0.5, 0.0), Tangram::Point2(1.0, 0.5),
    Tangram::Point2(0.5, 1.0), Tangram::Point2(0.0, 0.5) };
  
  //Check material ID correctness
  Tangram::MatPoly<2> square_matpoly;
  ASSERT_EQ(-1, square_matpoly.mat_id());
  square_matpoly.set_mat_id(mat_id);
  ASSERT_EQ(mat_id, square_matpoly.mat_id());
  square_matpoly.reset_mat_id();
  ASSERT_EQ(-1, square_matpoly.mat_id());
  
  //Initialization from ccw ordered vertices
  square_matpoly.initialize(square_points);
  
  //Verify coordinates
  const std::vector<Tangram::Point2>& matpoly_points = square_matpoly.matpoly_points();
  ASSERT_EQ(square_points.size(), square_matpoly.nvertices());
  for (int ivrt = 0; ivrt < square_points.size(); ivrt++)
    ASSERT_TRUE(approxEq(square_points[ivrt], matpoly_points[ivrt], 1.0e-15));
  
  //Verify faces
  ASSERT_EQ(square_faces.size(), square_matpoly.nfaces());
  for (int iface = 0; iface < square_faces.size(); iface++) {
    const std::vector<int>& face_vertices = square_matpoly.face_vertices(iface);
    ASSERT_EQ(2, face_vertices.size());
    ASSERT_EQ(square_faces[iface][0], face_vertices[0]);
    ASSERT_EQ(square_faces[iface][1], face_vertices[1]);
  }
  
  //Verify centroids
  for (int iface = 0; iface < square_faces.size(); iface++)
    ASSERT_TRUE(approxEq(face_centroids[iface],
                         square_matpoly.face_centroid(iface), 1.0e-15));
}
