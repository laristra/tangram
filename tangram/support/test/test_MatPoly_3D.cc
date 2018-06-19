/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#include "tangram/support/MatPoly.h"

#include "gtest/gtest.h"
#include "tangram/support/tangram.h"

/// Test the MatPoly structure for a 3D polyhedron

TEST(MatPoly, Mesh3D) {
  int mat_id = 1;
  
  //Test for a right triangular prism
  std::vector<Tangram::Point3> prism_points = {
    Tangram::Point3(1.0, 0.0, 0.0), Tangram::Point3(0.0, 1.0, 0.0),
    Tangram::Point3(0.0, 0.0, 0.0), Tangram::Point3(0.0, 1.0, 1.0),
    Tangram::Point3(1.0, 0.0, 1.0), Tangram::Point3(0.0, 0.0, 1.0) };
  std::vector< std::vector<int> > prism_faces = {
    {0, 2, 1}, {2, 0, 4, 5}, {3, 4, 0, 1}, {5, 3, 1, 2}, {3, 5, 4} };
  std::vector<Tangram::Point3> face_centroids = {
    Tangram::Point3(1.0/3.0, 1.0/3.0, 0.0), Tangram::Point3(0.5, 0.0, 0.5),
    Tangram::Point3(0.5, 0.5, 0.5), Tangram::Point3(0.0, 0.5, 0.5),
    Tangram::Point3(1.0/3.0, 1.0/3.0, 1.0) };
  
  //Check material ID correctness
  Tangram::MatPoly<3> prism_matpoly(mat_id);
  ASSERT_EQ(mat_id, prism_matpoly.mat_id());
  
  //Initialization
  prism_matpoly.initialize(prism_points, prism_faces);
  
  //Verify coordinates
  const std::vector<Tangram::Point3>& matpoly_points = prism_matpoly.points();
  ASSERT_EQ(prism_points.size(), prism_matpoly.num_vertices());
  for (int ivrt = 0; ivrt < prism_points.size(); ivrt++)
    ASSERT_TRUE(approxEq(prism_points[ivrt], matpoly_points[ivrt], 1.0e-15));
  
  //Verify faces
  ASSERT_EQ(prism_faces.size(), prism_matpoly.num_faces());
  for (int iface = 0; iface < prism_faces.size(); iface++) {
    const std::vector<int>& face_vertices = prism_matpoly.face_vertices(iface);
    ASSERT_EQ(prism_faces[iface].size(), face_vertices.size());
    for (int ivrt = 0; ivrt < prism_faces[iface].size(); ivrt++)
      ASSERT_EQ(prism_faces[iface][ivrt], face_vertices[ivrt]);
  }
  
  //Verify centroids
  for (int iface = 0; iface < prism_faces.size(); iface++)
    ASSERT_TRUE(approxEq(face_centroids[iface],
                         prism_matpoly.face_centroid(iface), 1.0e-15));
  
  std::vector<Tangram::Point3> faceted_prism_points = {
    Tangram::Point3(1.0, 0.0, 0.0), Tangram::Point3(0.0, 1.0, 0.0),
    Tangram::Point3(0.0, 0.0, 0.0), Tangram::Point3(0.0, 1.0, 1.0),
    Tangram::Point3(1.0, 0.0, 1.0), Tangram::Point3(0.0, 0.0, 1.0),
    Tangram::Point3(0.5, 0.0, 0.5), Tangram::Point3(0.5, 0.5, 0.5),
    Tangram::Point3(0.0, 0.5, 0.5),
  };
  std::vector< std::vector<int> > faceted_prism_faces = {
    {0, 2, 1},
    {6, 2, 0}, {6, 0, 4}, {6, 4, 5}, {6, 5, 2},
    {7, 3, 4}, {7, 4, 0}, {7, 0, 1}, {7, 1, 3},
    {8, 5, 3}, {8, 3, 1}, {8, 1, 2}, {8, 2, 5},
    {3, 5, 4} };
  
  //Create faceted poly
  Tangram::MatPoly<3> faceted_prism_matpoly;
  prism_matpoly.faceted_matpoly(&faceted_prism_matpoly);
  
  //Check material ID correctness
  ASSERT_EQ(mat_id, faceted_prism_matpoly.mat_id());
  
  //Verify facetization
  //Verify node coordinates
  const std::vector<Tangram::Point3>& faceted_matpoly_points =
    faceted_prism_matpoly.points();
  ASSERT_EQ(faceted_prism_points.size(), faceted_prism_matpoly.num_vertices());
  for (int ivrt = 0; ivrt < faceted_prism_points.size(); ivrt++)
    ASSERT_TRUE(approxEq(faceted_prism_points[ivrt],
                         faceted_matpoly_points[ivrt], 1.0e-15));
  
  //Verify facets
  ASSERT_EQ(faceted_prism_faces.size(), faceted_prism_matpoly.num_faces());
  for (int iface = 0; iface < faceted_prism_faces.size(); iface++) {
    const std::vector<int>& face_vertices = faceted_prism_matpoly.face_vertices(iface);
    ASSERT_EQ(faceted_prism_faces[iface].size(), face_vertices.size());
    for (int ivrt = 0; ivrt < faceted_prism_faces[iface].size(); ivrt++)
      ASSERT_EQ(faceted_prism_faces[iface][ivrt], face_vertices[ivrt]);
  }
  
  //Non-convex distorted prism
  std::vector<Tangram::Point3> ncv_prism_points = {
    Tangram::Point3(1.0, 0.0, 0.0), Tangram::Point3(0.4, 0.8, 0.2),
    Tangram::Point3(1.0, 1.0, 0.0), Tangram::Point3(0.0, 1.0, 0.0),
    Tangram::Point3(0.5, 0.1, 0.1), Tangram::Point3(0.0, 0.0, 0.0),
    Tangram::Point3(0.0, 0.0, 1.0), Tangram::Point3(1.0, 1.0, 1.0),
    Tangram::Point3(0.5, 0.9, 1.1), Tangram::Point3(1.0, 0.0, 1.0),
    Tangram::Point3(0.0, 1.0, 1.0), Tangram::Point3(0.6, 0.2, 0.8) };
  std::vector< std::vector<int> > ncv_prism_faces = {
    {3, 1, 2, 0, 4, 5},
    {5, 4, 11, 6}, {0, 9, 11, 4}, {7, 9, 0, 2},
    {2, 1, 8, 7}, {8, 1, 3, 10}, {5, 6, 10, 3},
    {7, 8, 10, 6, 11, 9} };
  std::vector<Tangram::Point3> ncv_prism_face_centroids = {
    Tangram::Point3(0.49982029799691324, 0.48793283780407681, 0.038845671672909921),
    Tangram::Point3(0.26077565588523149, 0.071880513370773835, 0.49415030095918006),
    Tangram::Point3(0.78729807432173626, 0.070061777159007799, 0.46574083274162237),
    Tangram::Point3(1.0, 0.5, 0.5),
    Tangram::Point3(0.72714164844017881, 0.92490808257951729, 0.55650182505587276),
    Tangram::Point3(0.22120243710721832, 0.92700420108481663, 0.58407100072352491),
    Tangram::Point3(0.0, 0.5, 0.5),
    Tangram::Point3(0.50063253857544732, 0.51180586493539493, 0.98662997695018306) };

  std::vector<double> ncv_prism_moments = {0.80774747387852563, 
    0.40401879823038706, 0.40692862723931095, 0.41298799905686806};  
  
  Tangram::MatPoly<3> ncv_prism_matpoly(mat_id);
  //Initialization
  ncv_prism_matpoly.initialize(ncv_prism_points, ncv_prism_faces);
  
  //Verify face centroids
  for (int iface = 0; iface < ncv_prism_faces.size(); iface++)
    ASSERT_TRUE(approxEq(ncv_prism_face_centroids[iface],
                         ncv_prism_matpoly.face_centroid(iface), 1.0e-15));

  //Verify moments
  std::vector<double> matpoly_moments = ncv_prism_matpoly.moments();
  for (int im = 0; im < 4; im++)
    ASSERT_NEAR(ncv_prism_moments[im], matpoly_moments[im], 1.0e-15);    

  //Moments of the faceted poly should be the same
  Tangram::MatPoly<3> faceted_ncv_prism;
  ncv_prism_matpoly.faceted_matpoly(&faceted_ncv_prism);
  std::vector<double> faceted_matpoly_moments = faceted_ncv_prism.moments();
  for (int im = 0; im < 4; im++)
    ASSERT_NEAR(ncv_prism_moments[im], faceted_matpoly_moments[im], 1.0e-15); 
}
