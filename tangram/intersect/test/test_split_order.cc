/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#include "tangram/intersect/split_r2d.h"
#include "tangram/intersect/split_r3d.h"

#include "gtest/gtest.h"
#include "tangram/support/tangram.h"

#define OUTPUT_TO_GMV

#ifdef OUTPUT_TO_GMV
  #include "tangram/driver/CellMatPoly.h"
  #include "tangram/driver/write_to_gmv.h"
#endif

#include <iostream>

template<int Dim>
Tangram::Plane_t<Dim> get_cutting_plane(Tangram::Point<Dim>& plane_pt, Tangram::Vector<Dim>& normal) 
{
  Tangram::Plane_t<Dim> cutting_plane; 
  cutting_plane.normal = normal;
  cutting_plane.normal.normalize();
  cutting_plane.dist2origin = - Tangram::dot(plane_pt.asV(), cutting_plane.normal);
  return cutting_plane; 
}

//This test illustrates that the notions of lower and upper half-spaces
//is inverted between R2D and R3D. 
/* 2D Test
*
*  *-----|-----*
*  |     |     |
*  | low |  up |
*  |     |     |
*  *-----|-----* 
*  O <---| Normal Direction
*/

TEST(split_order, ConvexPoly2D) {

   //Cutting plane in 2d
   Tangram::Point2 plane_pt2d(0.5,0.0);
   Tangram::Vector2 normal2d(-0.5,0.0);
   Tangram::Plane_t<2> cutting_plane2d = get_cutting_plane<2>(plane_pt2d, normal2d);
   std::cout<<"----SQUARE----"<<std::endl;
   std::cout<<"    Cutting Plane: normal = {"<<cutting_plane2d.normal[0]
            <<", "<<cutting_plane2d.normal[1]<<"}"<<std::endl;
   std::cout<<"    Cutting Plane: dist2origin = "<<cutting_plane2d.dist2origin<<std::endl;

   //r2d cutting plane
   r2d_plane r2d_cut_plane;
   for (int ixy = 0; ixy < 2; ixy++)
     r2d_cut_plane.n.xy[ixy] = cutting_plane2d.normal[ixy];
   r2d_cut_plane.d = cutting_plane2d.dist2origin;
  
   //Square Matpoly
   std::vector<Tangram::Point2> square_pnts = {
   Tangram::Point2(0.0, 0.0), Tangram::Point2(1.0, 0.0),
   Tangram::Point2(1.0, 1.0), Tangram::Point2(0.0, 1.0)};
   
   Tangram::MatPoly<2> square_matpoly; 
   square_matpoly.initialize(square_pnts);

   //Convert matpoly to r2d_poly and split
   r2d_poly r2dized_poly;
   r2d_poly r2d_subpolys[2];
   Tangram::matpoly_to_r2dpoly(square_matpoly, r2dized_poly);
   r2d_split(&r2dized_poly, 1, r2d_cut_plane, &r2d_subpolys[1], &r2d_subpolys[0]);

   //Reference splitted polygons and check
   std::vector<Tangram::Point2> square_ref_up ={
   Tangram::Point2(1.0, 0.0), Tangram::Point2(1.0, 1.0), 
   Tangram::Point2(0.5, 1.0), Tangram::Point2(0.5, 0.0)};
   
   std::vector<Tangram::Point2> square_ref_low ={
   Tangram::Point2(0.0, 0.0), Tangram::Point2(0.5, 0.0), 
   Tangram::Point2(0.5, 1.0), Tangram::Point2(0.0, 1.0)};
    
   //Get a MatPoly for a r2d subpoly
   Tangram::MatPoly<2> square_low, square_up;
   Tangram::r2dpoly_to_matpoly(r2d_subpolys[0], square_low);
   Tangram::r2dpoly_to_matpoly(r2d_subpolys[1], square_up);

   ///Check lower and upper subpolys
   for (int i = 0; i < square_low.num_vertices(); i++)
      ASSERT_TRUE(approxEq(square_low.vertex_point(i), square_ref_low[i], 1.0e-15));

   for (int i = 0; i < square_up.num_vertices(); i++)
     ASSERT_TRUE(approxEq(square_up.vertex_point(i), square_ref_up[i], 1.0e-15));
}

/* 3D Test
*  
*  *-----\---------*
*  | \     \        \
*  |   \     \        \
*  |    *----- \-------*
*  |    |       |      |
*  |    |       |      |
*  *    |       |      |
*   \   |  up   | low  |
*     \ |       |      |
*       *-------|------*
*      O       \|
*         <-----Normal Direction
*             
*/

TEST(split_order, ConvexPoly3D) {

   //Cutting plane in 3d
   Tangram::Point3 plane_pt3d(0.5, 0.0, 0.0);
   Tangram::Vector3 normal3d(-0.5, 0.0, 0.0);
   Tangram::Plane_t<3> cutting_plane3d = get_cutting_plane<3>(plane_pt3d, normal3d);
   std::cout<<"\n----CUBE----"<<std::endl;
   std::cout<<"    Cutting Plane: normal = {"<<cutting_plane3d.normal[0]
            <<", "<<cutting_plane3d.normal[1]<<", "<<cutting_plane3d.normal[2]<<"}"<<std::endl;
   std::cout<<"    Cutting Plane: dist2origin = "<<cutting_plane3d.dist2origin<<std::endl;

   //r3d cutting plane
   r3d_plane r3d_cut_plane;
   for (int ixyz = 0; ixyz < 3; ixyz++)
     r3d_cut_plane.n.xyz[ixyz] = cutting_plane3d.normal[ixyz];
   r3d_cut_plane.d = cutting_plane3d.dist2origin;
  
   //Cube Matpoly
   std::vector<Tangram::Point3> cube_pnts = {
   Tangram::Point3(0.0, 0.0, 0.0), Tangram::Point3(1.0, 0.0, 0.0),
   Tangram::Point3(1.0, 1.0, 0.0), Tangram::Point3(0.0, 1.0, 0.0),
   Tangram::Point3(0.0, 0.0, 1.0), Tangram::Point3(1.0, 0.0, 1.0),
   Tangram::Point3(1.0, 1.0, 1.0), Tangram::Point3(0.0, 1.0, 1.0)};
 
   std::vector<std::vector<int>> cube_faces = {{0,3,2,1},{0,1,5,4},{1,2,6,5},
                                               {2,3,7,6},{0,4,7,3},{4,5,6,7}}; 
 
   Tangram::MatPoly<3> cube_matpoly; 
   cube_matpoly.initialize(cube_pnts, cube_faces);

   //Convert matpoly to r3d_poly and split
   r3d_poly r3dized_poly;
   r3d_poly r3d_subpolys[2];
   Tangram::matpoly_to_r3dpoly(cube_matpoly, r3dized_poly);
   r3d_split(&r3dized_poly, 1, r3d_cut_plane, &r3d_subpolys[1], &r3d_subpolys[0]);

   //Reference splitted polygons and check
   std::vector<Tangram::Point3> cube_ref_low ={
   Tangram::Point3(1.0, 0.0, 0.0), Tangram::Point3(1.0, 1.0, 0.0), 
   Tangram::Point3(1.0, 0.0, 1.0), Tangram::Point3(1.0, 1.0, 1.0), 
   Tangram::Point3(0.5, 0.0, 0.0), Tangram::Point3(0.5, 1.0, 0.0), 
   Tangram::Point3(0.5, 0.0, 1.0), Tangram::Point3(0.5, 1.0, 1.0)}; 
   
   std::vector<Tangram::Point3> cube_ref_up ={
   Tangram::Point3(0.0, 0.0, 0.0), Tangram::Point3(0.0, 1.0, 0.0), 
   Tangram::Point3(0.0, 0.0, 1.0), Tangram::Point3(0.0, 1.0, 1.0), 
   Tangram::Point3(0.5, 0.0, 0.0), Tangram::Point3(0.5, 1.0, 0.0), 
   Tangram::Point3(0.5, 0.0, 1.0), Tangram::Point3(0.5, 1.0, 1.0)}; 
    
   //Get a MatPoly for a r2d subpoly
   std::vector<Tangram::MatPoly<3>> cube_low, cube_up;
   Tangram::r3dpoly_to_matpolys(r3d_subpolys[0], cube_low);
   Tangram::r3dpoly_to_matpolys(r3d_subpolys[1], cube_up);

   ASSERT_EQ(cube_low.size(), 1);
   ASSERT_EQ(cube_up.size(), 1);

   ///Check lower and upper subpolys
   for (int i = 0; i < cube_low[0].num_vertices(); i++)
      ASSERT_TRUE(approxEq(cube_low[0].vertex_point(i), cube_ref_low[i], 1.0e-15));

   for (int i = 0; i < cube_up[0].num_vertices(); i++)
     ASSERT_TRUE(approxEq(cube_up[0].vertex_point(i), cube_ref_up[i], 1.0e-15));
}
   
