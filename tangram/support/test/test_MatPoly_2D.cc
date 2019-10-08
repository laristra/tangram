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
  double dst_tol = sqrt(2)*std::numeric_limits<double>::epsilon();
  
  //Test for a unit square
  std::vector<Tangram::Point2> square_points = {
    Tangram::Point2(0.0, 0.0), Tangram::Point2(1.0, 0.0),
    Tangram::Point2(1.0, 1.0), Tangram::Point2(0.0, 1.0)
  };

  std::vector< std::vector<int> > square_faces = {
    {0, 1}, {1, 2}, {2, 3}, {3, 0}
  };

  std::vector<Tangram::Point2> face_centroids = {
    Tangram::Point2(0.5, 0.0), Tangram::Point2(1.0, 0.5),
    Tangram::Point2(0.5, 1.0), Tangram::Point2(0.0, 0.5)
  };
  
  //Check material ID correctness
  Tangram::MatPoly<2> square_matpoly;
  ASSERT_EQ(-1, square_matpoly.mat_id());
  square_matpoly.set_mat_id(mat_id);
  ASSERT_EQ(mat_id, square_matpoly.mat_id());
  square_matpoly.reset_mat_id();
  ASSERT_EQ(-1, square_matpoly.mat_id());
  
  //Initialization from ccw ordered vertices
  square_matpoly.initialize(square_points, dst_tol);

  auto const& matpoly_points = square_matpoly.points();
  int const nb_regu_points   = square_points.size();
  int const nb_regu_faces    = square_faces.size();
  int const nb_poly_points   = square_matpoly.num_vertices();
  int const nb_poly_faces    = square_matpoly.num_faces();

  ASSERT_EQ(nb_regu_points, nb_poly_points);
  ASSERT_EQ(nb_regu_faces, nb_poly_faces);

  // Verify coordinates
  for (int ivrt = 0; ivrt < nb_regu_points; ivrt++) {
    auto const& regu_point = square_points[ivrt];
    auto const& poly_point = matpoly_points[ivrt];
    ASSERT_TRUE(approxEq(regu_point, poly_point, 1.0e-15));
  }

  // Verify faces
  for (int iface = 0; iface < nb_regu_faces; iface++) {
    auto const& face_vertices = square_matpoly.face_vertices(iface);
    ASSERT_EQ(unsigned(2), face_vertices.size());
    ASSERT_EQ(square_faces[iface][0], face_vertices[0]);
    ASSERT_EQ(square_faces[iface][1], face_vertices[1]);
  }
  
  // Verify centroids
  for (int iface = 0; iface < nb_regu_faces; iface++) {
    auto const& regu_centroid = face_centroids[iface];
    auto const& poly_centroid = square_matpoly.face_centroid(iface);
    ASSERT_TRUE(approxEq(regu_centroid, poly_centroid, 1.0e-15));
  }

  // Test moments for a non-convex pentagon
  std::vector<Tangram::Point2> ncv_poly_points = {
    Tangram::Point2(0.0, 0.0), Tangram::Point2(1.0, 0.0),
    Tangram::Point2(1.0, 1.0), Tangram::Point2(0.5, 0.5),
    Tangram::Point2(0.0, 1.0)
  };

  std::vector<double> ncv_poly_moments = { 0.75, 0.375, 21.0/72.0 };
  Tangram::MatPoly<2> ncv_matpoly(mat_id);
  ncv_matpoly.initialize(ncv_poly_points, dst_tol);

  // Verify moments
  auto matpoly_moments = ncv_matpoly.moments();
  for (int im = 0; im < 3; im++)
    ASSERT_NEAR(ncv_poly_moments[im], matpoly_moments[im], 1.0e-15);
}
