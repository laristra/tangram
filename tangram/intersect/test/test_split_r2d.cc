/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#include "tangram/intersect/split_r2d.h"

#include "gtest/gtest.h"
#include "tangram/support/tangram.h"

#define OUTPUT_TO_GMV

#ifdef OUTPUT_TO_GMV
  #include "tangram/driver/CellMatPoly.h"
  #include "tangram/driver/write_to_gmv.h"
#endif

#include <iostream>

// Test methods in split_r2d.h
// If OUTPUT_TO_GMV is defined, material poly's resulting from splitting
// will be written to a gmv file.

TEST(split_r2d, ConvexPoly) {
  //Create a single convex polygon
  std::vector<Tangram::Point2> hexagon_pnts = {
  Tangram::Point2(1.0, 0.0), Tangram::Point2(2.0, 0.0),
  Tangram::Point2(3.0, 1.0), Tangram::Point2(2.0, 2.0),
  Tangram::Point2(1.0, 2.0), Tangram::Point2(0.0, 1.0) };
   
  //Cutting plane
  Tangram::Plane_t<2> cutting_plane; 
  cutting_plane.normal = Tangram::Vector2(-1.5,-1.0);
  cutting_plane.dist2origin = 2.25; 

  //Construct reference split polygons
  std::vector<std::vector<Tangram::Point2>> ref_cp_pnts(2);
  ref_cp_pnts[0] = {
  Tangram::Point2(1.0,0.0), Tangram::Point2(1.5,0.0),
  Tangram::Point2(0.5,1.5), Tangram::Point2(0.0,1.0) };

  ref_cp_pnts[1] = {
  Tangram::Point2(1.5,0.0), Tangram::Point2(2.0,0.0),
  Tangram::Point2(3.0,1.0), Tangram::Point2(2.0,2.0),
  Tangram::Point2(1.0,2.0), Tangram::Point2(0.5,1.5) };

  //Create Matpoly corresponding to the poly
  Tangram::MatPoly<2> hexagon_matpoly; 
  hexagon_matpoly.initialize(hexagon_pnts);
  std::vector<Tangram::MatPoly<2>> cpmatpolys = {hexagon_matpoly};

  //Split using the split routine directly
  Tangram::MatPoly<2> convex_polys[2];
  std::vector<double> cp_moments[2];
  split_convex_matpoly_r2d(hexagon_matpoly, cutting_plane,
                           convex_polys[0], convex_polys[1],
                           cp_moments[0], cp_moments[1]);

  //Check that vertices of polygons below and above the plane match
  for (int ihs = 0; ihs < 2; ihs++) {
    std::cout<<"CP "<<ihs<<" : NV =  "<<convex_polys[ihs].num_vertices()<<std::endl;    
    
    for (int ivrt = 0; ivrt < convex_polys[ihs].num_vertices(); ivrt++)
    {
       Tangram::Point2 pt = convex_polys[ihs].vertex_point(ivrt); 
       std::cout<<"pts = {"<<pt[0]<<", "<<pt[1]<<"}"<<std::endl;
    }

    //ASSERT_EQ(ref_cp_pnts[ihs].size(), convex_polys[ihs].num_vertices());
   // for (int ivrt = 0; ivrt < ref_cp_pnts[ihs].size(); ivrt++)
   //   ASSERT_TRUE(approxEq(ref_cp_pnts[ihs][ivrt],
   //                        convex_polys[ihs].vertex_point(ivrt), 1.0e-15));
  }

  //Split using SplitR2D class 
  Tangram::SplitR2D split(cpmatpolys, cutting_plane, true);
  Tangram::HalfSpaceSets_t<2> hsp_sets = split();

  //Check
  ASSERT_EQ(hsp_sets.lower_halfspace_set.matpolys.size(), 1);
  ASSERT_EQ(hsp_sets.upper_halfspace_set.matpolys.size(), 1);

  Tangram::MatPoly<2>* hsp_ptrs[2] = {&hsp_sets.lower_halfspace_set.matpolys[0],
                                          &hsp_sets.upper_halfspace_set.matpolys[0]};
  
  //Check that we obtained the same MatPoly's as before                                        
  for (int ihs = 0; ihs < 2; ihs++) {
    std::cout<<"SP-CP "<<ihs<<" : NV =  "<<convex_polys[ihs].num_vertices()<<std::endl;    

    for (int ivrt = 0; ivrt < hsp_ptrs[ihs]->num_vertices(); ivrt++)
    {
       Tangram::Point2 pt = hsp_ptrs[ihs]->vertex_point(ivrt); 
       std::cout<<"pts = {"<<pt[0]<<", "<<pt[1]<<"}"<<std::endl;
    }
   /*
    ASSERT_EQ(ref_cp_pnts[ihs].size(), hsp_ptrs[ihs]->num_vertices());
    for (int ivrt = 0; ivrt < ref_cp_pnts[ihs].size(); ivrt++)
      ASSERT_TRUE(approxEq(ref_cp_pnts[ihs][ivrt],
                           hsp_ptrs[ihs]->vertex_point(ivrt), 1.0e-15));
   */
   }

#ifdef OUTPUT_TO_GMV
  std::vector<std::shared_ptr<Tangram::CellMatPoly<2>>> cellmatpoly_list;
  cellmatpoly_list.push_back(std::make_shared< Tangram::CellMatPoly<2> >(0));
  for (int ihs = 0; ihs < 2; ihs++) {
    int nverts = convex_polys[ihs].num_vertices();
    std::vector<Tangram::Point2> vertices = convex_polys[ihs].points(); 
   
    cellmatpoly_list[0]->add_matpoly(ihs, convex_polys[ihs].num_vertices(),
                                     &convex_polys[ihs].points()[0],
                                     nullptr, nullptr,
                                     nullptr, nullptr);
  }
  Tangram::write_to_gmv(cellmatpoly_list, "hexagon_matpolys.gmv");
#endif


  //Check 
  //
}

TEST(split_r2d, NonConvexPoly) {

  //Create a single non-convex polygon
  std::vector<Tangram::Point2> cat_pnts = {
  Tangram::Point2(1.0, 0.0), Tangram::Point2(4.0, 0.0),
  Tangram::Point2(2.0, 2.0), Tangram::Point2(4.0, 4.0),
  Tangram::Point2(1.0, 4.0)};
   
  //Cutting plane
  Tangram::Plane_t<2> cutting_plane; 
  cutting_plane.normal = Tangram::Vector2(-1.0,0.0);
  cutting_plane.dist2origin = 3.0; 

  //Construct reference split polygons
  //FIX IT: NEED THE CORRECT REFERENCES 
  /*std::vector<std::vector<Tangram::Point2>> ref_cp_pnts(2);
  ref_cp_pnts[0] = {
  Tangram::Point2(1.0,0.0), Tangram::Point2(1.5,0.0),
  Tangram::Point2(0.5,1.5), Tangram::Point2(0.0,1.0) };

  ref_cp_pnts[1] = {
  Tangram::Point2(1.5,0.0), Tangram::Point2(2.0,0.0),
  Tangram::Point2(3.0,1.0), Tangram::Point2(2.0,2.0),
  Tangram::Point2(1.0,2.0), Tangram::Point2(0.5,1.5) };
*/
  //Create Matpoly corresponding to the poly
  Tangram::MatPoly<2> cat_matpoly; 
  cat_matpoly.initialize(cat_pnts);
  std::vector<Tangram::MatPoly<2>> cpmatpolys = {cat_matpoly};

  //Split using SplitR2D class 
  Tangram::SplitR2D split(cpmatpolys, cutting_plane, false);
  Tangram::HalfSpaceSets_t<2> hsp_sets = split();

  //Check
  //ASSERT_EQ(hsp_sets.lower_halfspace_set.matpolys.size(), 1);
  //ASSERT_EQ(hsp_sets.upper_halfspace_set.matpolys.size(), 1);

  std::vector<Tangram::MatPoly<2>> hsp_lower_matpolys = hsp_sets.lower_halfspace_set.matpolys;
  std::vector<Tangram::MatPoly<2>> hsp_upper_matpolys = hsp_sets.upper_halfspace_set.matpolys;

  //Check that we obtained the same MatPoly's as before                                        
  std::cout<<"HSP_LOWER_MATPOLYS"<<std::endl;
  for (int ihs = 0; ihs < hsp_lower_matpolys.size(); ihs++) {
    std::cout<<"CP "<<ihs<<" : NV =  "<<hsp_lower_matpolys[ihs].num_vertices()<<std::endl;    

    for (int ivrt = 0; ivrt < hsp_lower_matpolys[ihs].num_vertices(); ivrt++)
    {
       Tangram::Point2 pt = hsp_lower_matpolys[ihs].vertex_point(ivrt); 
       std::cout<<"pts = {"<<pt[0]<<", "<<pt[1]<<"}"<<std::endl;
    }
  }
  std::cout<<"HSP_UPPER_MATPOLYS"<<std::endl;
  for (int ihs = 0; ihs < hsp_upper_matpolys.size(); ihs++) {
    std::cout<<"CP "<<ihs<<" : NV =  "<<hsp_upper_matpolys[ihs].num_vertices()<<std::endl;    

    for (int ivrt = 0; ivrt < hsp_upper_matpolys[ihs].num_vertices(); ivrt++)
    {
       Tangram::Point2 pt = hsp_upper_matpolys[ihs].vertex_point(ivrt); 
       std::cout<<"pts = {"<<pt[0]<<", "<<pt[1]<<"}"<<std::endl;
    }
   /*
    ASSERT_EQ(ref_cp_pnts[ihs].size(), hsp_ptrs[ihs]->num_vertices());
    for (int ivrt = 0; ivrt < ref_cp_pnts[ihs].size(); ivrt++)
      ASSERT_TRUE(approxEq(ref_cp_pnts[ihs][ivrt],
                           hsp_ptrs[ihs]->vertex_point(ivrt), 1.0e-15));
   */
   }

#ifdef OUTPUT_TO_GMV
/*
  std::vector<std::shared_ptr<Tangram::CellMatPoly<2>>> cellmatpoly_list;
  cellmatpoly_list.push_back(std::make_shared< Tangram::CellMatPoly<2> >(0));
  for (int ihs = 0; ihs < 2; ihs++) {
    int nverts = convex_polys[ihs].num_vertices();
    std::vector<Tangram::Point2> vertices = convex_polys[ihs].points(); 
   
    cellmatpoly_list[0]->add_matpoly(ihs, convex_polys[ihs].num_vertices(),
                                     &convex_polys[ihs].points()[0],
                                     nullptr, nullptr,
                                     nullptr, nullptr);
  }
  Tangram::write_to_gmv(cellmatpoly_list, "hexagon_matpolys.gmv");*/
#endif


  //Check 
  //
}
