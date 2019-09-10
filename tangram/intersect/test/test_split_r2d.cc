/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#include <iostream>
#include <vector>

#include "gtest/gtest.h"

// tangram includes
#define OUTPUT_TO_GMV

#ifdef OUTPUT_TO_GMV
  #include "tangram/driver/CellMatPoly.h"
  #include "tangram/driver/write_to_gmv.h"
#endif

#include "tangram/intersect/split_r2d.h"
#include "tangram/support/tangram.h"

// wonton includes
#include "wonton/support/Vector.h"

enum ptype {
  CONVEX_SINGLEPOLY,
  NONCONVEX_SINGLEPOLY
};

void matpoly_cases(ptype POLYTYPE, std::vector<Tangram::MatPoly<2>>& matpolys,
                   double dst_tol)
{
  matpolys.clear();

  switch (POLYTYPE)
  {
    case CONVEX_SINGLEPOLY:
      {
	std::vector<Tangram::Point2> hexagon_pnts = {
          Tangram::Point2(1.0, 0.0), Tangram::Point2(2.0, 0.0),
          Tangram::Point2(3.0, 1.0), Tangram::Point2(2.0, 2.0),
          Tangram::Point2(1.0, 2.0), Tangram::Point2(0.0, 1.0) };

        Tangram::MatPoly<2> hexagon_matpoly;
        hexagon_matpoly.initialize(hexagon_pnts, dst_tol);
        matpolys.emplace_back(hexagon_matpoly);

	break;
      }
    case NONCONVEX_SINGLEPOLY:
      {
        std::vector<Tangram::Point2> cat_pnts = {
          Tangram::Point2(1.0, 0.0), Tangram::Point2(4.0, 0.0),
          Tangram::Point2(3.0, 2.0), Tangram::Point2(4.0, 4.0),
          Tangram::Point2(1.0, 4.0)};

        Tangram::MatPoly<2> cat_matpoly;
        cat_matpoly.initialize(cat_pnts, dst_tol);
        matpolys.emplace_back(cat_matpoly);

	break;
      }
    default:
      {
        std::cerr<<"Requesting test matpolys for non-supported configuration"<<std::endl;
      }

  }
}

Tangram::Plane_t<2> get_cutting_plane(Tangram::Point<2>& plane_pt, Tangram::Vector<2>& normal)
{
  Tangram::Plane_t<2> cutting_plane;
  cutting_plane.normal = normal;
  cutting_plane.normal.normalize();
  cutting_plane.dist2origin = - Wonton::dot(plane_pt.asV(), cutting_plane.normal);
  return cutting_plane;
}

void reference_matpolys(ptype POLYTYPE,
                        std::vector<std::vector<Tangram::Point2>>& ref_pnts_lower,
                        std::vector<std::vector<Tangram::Point2>>& ref_pnts_upper)
{

  switch (POLYTYPE)
  {
    case CONVEX_SINGLEPOLY:
      {
        // Upper polygons
	ref_pnts_upper[0] = {
          Tangram::Point2(1.0,0.0), Tangram::Point2(1.5,0.0),
          Tangram::Point2(0.5,1.5), Tangram::Point2(0.0,1.0) };

        // Lower polygons
	ref_pnts_lower[0] = {
          Tangram::Point2(2.0,0.0), Tangram::Point2(3.0,1.0),
          Tangram::Point2(2.0,2.0), Tangram::Point2(1.0,2.0),
          Tangram::Point2(0.5,1.5), Tangram::Point2(1.5,0.0) };

	break;
      }
    case NONCONVEX_SINGLEPOLY:
      {
        // Upper polygons
	ref_pnts_upper[0] = {
          Tangram::Point2(1.0,0.0), Tangram::Point2(3.5,0.0),
          Tangram::Point2(3.5,0.57692307692), Tangram::Point2(2.2666666667,2.0)};

	ref_pnts_upper[1] = {
          Tangram::Point2(3.0,2.0), Tangram::Point2(2.2666666667,2.0),
          Tangram::Point2(3.5,0.57692307692), Tangram::Point2(3.5,1.0) };

	ref_pnts_upper[2] = {
          Tangram::Point2(3.0,2.0), Tangram::Point2(3.5,3.0),
          Tangram::Point2(3.5,3.4230769231), Tangram::Point2(2.2666666667,2.0) };

	ref_pnts_upper[3] = {
          Tangram::Point2(1.0,4.0), Tangram::Point2(2.2666666667,2.0),
          Tangram::Point2(3.5,3.4230769231) , Tangram::Point2(3.5,4.0) };

	ref_pnts_upper[4] = {
          Tangram::Point2(1.0,4.0), Tangram::Point2(1.0,0.0), Tangram::Point2(2.2666666667,2.0)};

	// Lower polygons
	ref_pnts_lower[0] = {
          Tangram::Point2(4.0,0.0), Tangram::Point2(3.5,0.57692307692), Tangram::Point2(3.5,0.0) };

	ref_pnts_lower[1] = {
          Tangram::Point2(4.0,0.0), Tangram::Point2(3.5,1.0), Tangram::Point2(3.5,0.57692307692) };

	ref_pnts_lower[2] = {
          Tangram::Point2(4.0,4.0), Tangram::Point2(3.5,3.4230769231), Tangram::Point2(3.5,3.0) };

	ref_pnts_lower[3] = {
          Tangram::Point2(4.0,4.0), Tangram::Point2(3.5,4.0), Tangram::Point2(3.5,3.4230769231) };

        break;
      }
    default:
      {
        std::cerr<<"Requesting test matpolys for non-supported configuration"<<std::endl;
      }
  }
}

// Test methods in split_r2d.h
// If OUTPUT_TO_GMV is defined, material poly's resulting from splitting
// will be written to a gmv file.

TEST(split_r2d, ConvexPoly) {
  double dst_tol = sqrt(2)*std::numeric_limits<double>::epsilon();
  double vol_tol = pow(dst_tol, 2);

  //Create a single convex polygon
  std::vector<Tangram::MatPoly<2>> convex_singlepoly;
  matpoly_cases(CONVEX_SINGLEPOLY, convex_singlepoly, dst_tol);

  //Cutting plane
  Tangram::Point2 plane_pt(1.5,0.0);
  Tangram::Vector2 normal(-1.5,-1.0);
  Tangram::Plane_t<2> cutting_plane = get_cutting_plane(plane_pt, normal);

  //Construct reference split polygons
  std::vector<std::vector<Tangram::Point2>> ref_cp_pnts_lower(1);
  std::vector<std::vector<Tangram::Point2>> ref_cp_pnts_upper(1);
  reference_matpolys(CONVEX_SINGLEPOLY, ref_cp_pnts_lower, ref_cp_pnts_upper);

  //Split using SplitR2D class
  Tangram::SplitR2D split(convex_singlepoly, cutting_plane, vol_tol, dst_tol, true);
  Tangram::HalfSpaceSets_t<2> hsp_sets = split();

  //Check
  ASSERT_EQ(hsp_sets.lower_halfspace_set.matpolys.size(), 1);
  ASSERT_EQ(hsp_sets.upper_halfspace_set.matpolys.size(), 1);

  std::vector<Tangram::MatPoly<2>> hsp_lower_matpolys = hsp_sets.lower_halfspace_set.matpolys;
  std::vector<Tangram::MatPoly<2>> hsp_upper_matpolys = hsp_sets.upper_halfspace_set.matpolys;

  //Check that we obtained the same MatPoly's as before
  ASSERT_EQ(ref_cp_pnts_lower[0].size(), hsp_lower_matpolys[0].num_vertices());
  for (int ivrt = 0; ivrt < ref_cp_pnts_lower[0].size(); ivrt++)
    ASSERT_TRUE(approxEq(ref_cp_pnts_lower[0][ivrt],
                         hsp_lower_matpolys[0].vertex_point(ivrt), 1.0e-15));
  ASSERT_EQ(ref_cp_pnts_upper[0].size(), hsp_upper_matpolys[0].num_vertices());
  for (int ivrt = 0; ivrt < ref_cp_pnts_upper[0].size(); ivrt++)
    ASSERT_TRUE(approxEq(ref_cp_pnts_upper[0][ivrt],
                         hsp_upper_matpolys[0].vertex_point(ivrt), 1.0e-15));

#ifdef OUTPUT_TO_GMV
  std::vector<std::shared_ptr<Tangram::CellMatPoly<2>>> cellmatpoly_list;
  cellmatpoly_list.push_back(std::make_shared< Tangram::CellMatPoly<2> >(0));
  for (int ihs = 0; ihs < hsp_lower_matpolys.size(); ihs++) {
    int nverts = hsp_lower_matpolys[ihs].num_vertices();
    std::vector<Tangram::Point2> vertices = hsp_lower_matpolys[ihs].points();

    cellmatpoly_list[0]->add_matpoly(ihs, nverts, vertices.data(), dst_tol,
                                     nullptr, nullptr,
                                     nullptr, nullptr);
  }
  Tangram::write_to_gmv(cellmatpoly_list, "convex_sp_lower.gmv");
  cellmatpoly_list.clear();
  cellmatpoly_list.push_back(std::make_shared< Tangram::CellMatPoly<2> >(0));
  for (int ihs = 0; ihs < hsp_upper_matpolys.size(); ihs++) {
    int nverts = hsp_upper_matpolys[ihs].num_vertices();
    std::vector<Tangram::Point2> vertices = hsp_upper_matpolys[ihs].points();

    cellmatpoly_list[0]->add_matpoly(ihs, nverts, vertices.data(), dst_tol,
                                     nullptr, nullptr,
                                     nullptr, nullptr);
  }
  Tangram::write_to_gmv(cellmatpoly_list, "convex_sp_upper.gmv");

#endif
}

TEST(split_r2d, NonConvexPoly) {
  double dst_tol = sqrt(2)*std::numeric_limits<double>::epsilon();
  double vol_tol = pow(dst_tol, 2);

  //Create a single convex polygon
  std::vector<Tangram::MatPoly<2>> cpmatpolys;
  matpoly_cases(NONCONVEX_SINGLEPOLY, cpmatpolys, dst_tol);

  //Cutting plane
  Tangram::Point2 plane_pt(3.5,0.0);
  Tangram::Vector2 normal(-1.0,0.0);
  Tangram::Plane_t<2> cutting_plane = get_cutting_plane(plane_pt, normal);

  //Construct reference split polygons
  std::vector<std::vector<Tangram::Point2>> ref_ncp_pnts_upper(5);
  std::vector<std::vector<Tangram::Point2>> ref_ncp_pnts_lower(4);
  reference_matpolys(NONCONVEX_SINGLEPOLY, ref_ncp_pnts_lower, ref_ncp_pnts_upper);

  //Split using SplitR2D class
  Tangram::SplitR2D split(cpmatpolys, cutting_plane, vol_tol, dst_tol, false);
  Tangram::HalfSpaceSets_t<2> hsp_sets = split();

  //Check
  ASSERT_EQ(hsp_sets.lower_halfspace_set.matpolys.size(), 4);
  ASSERT_EQ(hsp_sets.upper_halfspace_set.matpolys.size(), 5);

  std::vector<Tangram::MatPoly<2>> hsp_lower_matpolys = hsp_sets.lower_halfspace_set.matpolys;
  std::vector<Tangram::MatPoly<2>> hsp_upper_matpolys = hsp_sets.upper_halfspace_set.matpolys;

  //Check
  for (int ihs = 0; ihs < hsp_lower_matpolys.size(); ihs++) {
    ASSERT_EQ(ref_ncp_pnts_lower[ihs].size(), hsp_lower_matpolys[ihs].num_vertices());
    for (int ivrt = 0; ivrt < ref_ncp_pnts_lower[ihs].size(); ivrt++) {
      ASSERT_TRUE(approxEq(ref_ncp_pnts_lower[ihs][ivrt],
                           hsp_lower_matpolys[ihs].vertex_point(ivrt), 1.0e-10));
    }
  }
  for (int ihs = 0; ihs < hsp_upper_matpolys.size(); ihs++) {
    ASSERT_EQ(ref_ncp_pnts_upper[ihs].size(), hsp_upper_matpolys[ihs].num_vertices());
    for (int ivrt = 0; ivrt < ref_ncp_pnts_upper[ihs].size(); ivrt++){
      ASSERT_TRUE(approxEq(ref_ncp_pnts_upper[ihs][ivrt],
                           hsp_upper_matpolys[ihs].vertex_point(ivrt), 1.0e-10));
    }
  }

#ifdef OUTPUT_TO_GMV
  std::vector<std::shared_ptr<Tangram::CellMatPoly<2>>> cellmatpoly_list;
  cellmatpoly_list.push_back(std::make_shared< Tangram::CellMatPoly<2> >(0));
  for (int ihs = 0; ihs < hsp_lower_matpolys.size(); ihs++) {
    int nverts = hsp_lower_matpolys[ihs].num_vertices();
    std::vector<Tangram::Point2> vertices = hsp_lower_matpolys[ihs].points();

    cellmatpoly_list[0]->add_matpoly(ihs, nverts, vertices.data(), dst_tol,
                                     nullptr, nullptr,
                                     nullptr, nullptr);
  }
  Tangram::write_to_gmv(cellmatpoly_list, "nonconvex_sp_lower.gmv");
  cellmatpoly_list.clear();
  cellmatpoly_list.push_back(std::make_shared< Tangram::CellMatPoly<2> >(0));
  for (int ihs = 0; ihs < hsp_upper_matpolys.size(); ihs++) {
    int nverts = hsp_upper_matpolys[ihs].num_vertices();
    std::vector<Tangram::Point2> vertices = hsp_upper_matpolys[ihs].points();

    cellmatpoly_list[0]->add_matpoly(ihs, nverts, vertices.data(), dst_tol,
                                     nullptr, nullptr,
                                     nullptr, nullptr);
  }
  Tangram::write_to_gmv(cellmatpoly_list, "nonconvex_sp_upper.gmv");
#endif
}

TEST(clip_r2d, ConvexPoly) {
  double dst_tol = sqrt(2)*std::numeric_limits<double>::epsilon();

  //Create a single convex polygon
  std::vector<Tangram::MatPoly<2>> cpmatpolys;
  matpoly_cases(CONVEX_SINGLEPOLY, cpmatpolys, dst_tol);

  //Cutting plane
  Tangram::Point2 plane_pt(1.5,0.0);
  Tangram::Vector2 normal(1.5,1.0);
  Tangram::Plane_t<2> cutting_plane = get_cutting_plane(plane_pt, normal);

  //Reference
  std::vector<double> ref_moments = {0.875, 0.625, 0.604166666666666667};

  //Clip using ClipR2D class
  Tangram::ClipR2D clip(dst_tol);
  clip.set_matpolys(cpmatpolys, true);
  clip.set_plane(cutting_plane);

  std::vector<double> agmoments = clip();
  ASSERT_EQ(agmoments.size(),3);

  //Check
  for (int i = 0; i < 3; i++)
    ASSERT_NEAR(agmoments[i], ref_moments[i], 1.0e-15);
}

TEST(clip_r2d, NonConvexPoly) {
  double dst_tol = sqrt(2)*std::numeric_limits<double>::epsilon();

  //Create a single convex polygon
  std::vector<Tangram::MatPoly<2>> cpmatpolys;
  matpoly_cases(NONCONVEX_SINGLEPOLY, cpmatpolys, dst_tol);

  //Cutting plane
  Tangram::Point2 plane_pt(3.5,0.0);
  Tangram::Vector2 normal(-1.0,0.0);
  Tangram::Plane_t<2> cutting_plane = get_cutting_plane(plane_pt, normal);

  //Reference
  std::vector<double> ref_moments = {0.5, 1.83333333333333333, 1.0};

  //Clip using ClipR2D class
  Tangram::ClipR2D clip(dst_tol);
  clip.set_matpolys(cpmatpolys, true);
  clip.set_plane(cutting_plane);

  std::vector<double> agmoments = clip();
  ASSERT_EQ(agmoments.size(),3);

  //Check
  for (int i = 0; i < 3; i++)
    ASSERT_NEAR(agmoments[i], ref_moments[i], 1.0e-15);
}
