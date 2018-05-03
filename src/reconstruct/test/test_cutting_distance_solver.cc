/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#include "tangram/reconstruct/cutting_distance_solver.h"
#include "tangram/intersect/split_r3d.h"

#include "gtest/gtest.h"
#include "tangram/support/tangram.h"

// Test the method for finding the cutting distance in cutting_distance_solver.h

TEST(cutting_distance_solver, Mesh3D) {
  //Non-convex polyhedron with flat faces
  std::vector<Tangram::Point3> ncv_poly_points = {
    Tangram::Point3(0.0, 0.0, 0.0), Tangram::Point3(1.0, 0.0, 0.0),
    Tangram::Point3(1.0, 1.0, 0.0), Tangram::Point3(0.0, 1.0, 0.0),
    Tangram::Point3(0.4, 0.4, 0.2), Tangram::Point3(0.6, 0.4, 0.2),
    Tangram::Point3(0.6, 0.6, 0.2), Tangram::Point3(0.4, 0.6, 0.2),
    Tangram::Point3(0.0, 0.0, 1.0), Tangram::Point3(1.0, 0.0, 1.0),
    Tangram::Point3(1.0, 1.0, 1.0), Tangram::Point3(0.0, 1.0, 1.0),
    Tangram::Point3(0.4, 0.4, 0.8), Tangram::Point3(0.6, 0.4, 0.8),
    Tangram::Point3(0.6, 0.6, 0.8), Tangram::Point3(0.4, 0.6, 0.8) };
  std::vector< std::vector<int> > ncv_poly_faces = {
    {0, 1, 9, 8}, {1, 2, 10, 9}, {2, 3, 11, 10}, {3, 0, 8, 11},
    {4, 7, 6, 5}, {0, 4, 5, 1}, {1, 5, 6, 2}, {2, 6, 7, 3},
    {3, 7, 4, 0}, {12, 13, 14, 15}, {8, 9, 13, 12}, {9, 10, 14, 13},
    {10, 11, 15, 14}, {11, 8, 12, 15} };

  std::vector<Tangram::Point3> ref_plane_pts = {
    Tangram::Point3(0.0, 0.0, 0.4), Tangram::Point3(1.0, 0.0, 0.6),
    Tangram::Point3(0.0, 1.0, 0.5) };

  Tangram::Plane_t<3> ref_cutting_plane;
  //Oblique cutting plane
  ref_cutting_plane.normal = Tangram::Vector3(-0.2, -0.1, 1.0);
  ref_cutting_plane.normal.normalize();
  ref_cutting_plane.dist2origin = -dot(ref_cutting_plane.normal, ref_plane_pts[0].asV());

  double target_volume = 701.0/1500.0;

  Tangram::IterativeMethodTolerances_t tols = {
    .max_num_iter = 1000, .arg_eps = 1.0e-15, .fun_eps = 1.0e-15};  

  std::vector< Tangram::MatPoly<3> > ncv_matpolys(1);
  //Initialize non-convex MatPoly
  ncv_matpolys[0].initialize(ncv_poly_points, ncv_poly_faces);
  double ncv_poly_volume = ncv_matpolys[0].moments()[0];

  //Decomposition into convex subpolys
  std::vector< Tangram::MatPoly<3> > cv_subpolys;
  ncv_matpolys[0].decompose(cv_subpolys);

  //Create cutting distance solvers
  Tangram::CuttingDistanceSolver<3, Tangram::ClipR3D> 
    solve_cut_dst_ncv(ncv_matpolys, ref_cutting_plane.normal, tols, true);

  Tangram::CuttingDistanceSolver<3, Tangram::ClipR3D> 
    solve_cut_dst_cv(cv_subpolys, ref_cutting_plane.normal, tols, true);

  // Find distance to origin for the cutting plane, non-convex poly
  solve_cut_dst_ncv.set_target_volume(target_volume);
  std::vector<double> clip_res = solve_cut_dst_ncv();

  ASSERT_NEAR(ref_cutting_plane.dist2origin, clip_res[0], 1.0e-15);

  // Find distance to origin for the cutting plane, convex decomposition
  solve_cut_dst_cv.set_target_volume(target_volume);
  clip_res = solve_cut_dst_cv();

  ASSERT_NEAR(ref_cutting_plane.dist2origin, clip_res[0], 1.0e-15);

  //Cutting plane passing through the centroid
  ref_cutting_plane.normal = Tangram::Vector3(-1.0, 0.0, 0.0);
  ref_cutting_plane.dist2origin = 0.5;

  solve_cut_dst_ncv.set_target_volume(0.5*ncv_poly_volume);
  clip_res = solve_cut_dst_ncv();

  ASSERT_NEAR(ref_cutting_plane.dist2origin, clip_res[0], 1.0e-15);

  solve_cut_dst_cv.set_target_volume(0.5*ncv_poly_volume);
  clip_res = solve_cut_dst_cv();

  ASSERT_NEAR(ref_cutting_plane.dist2origin, clip_res[0], 1.0e-15);

  //Cutting plane passing through the centroid and edges
  ref_cutting_plane.normal = Tangram::Vector3(1.0, 1.0, 0.0);
  ref_cutting_plane.normal.normalize();
  ref_cutting_plane.dist2origin = -0.5*sqrt(2);

  solve_cut_dst_ncv.set_target_volume(0.5*ncv_poly_volume);
  clip_res = solve_cut_dst_ncv();

  ASSERT_NEAR(ref_cutting_plane.dist2origin, clip_res[0], 1.0e-15);

  solve_cut_dst_cv.set_target_volume(0.5*ncv_poly_volume);
  clip_res = solve_cut_dst_cv();

  ASSERT_NEAR(ref_cutting_plane.dist2origin, clip_res[0], 1.0e-15);
}