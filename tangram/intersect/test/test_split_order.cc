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

// Test methods in split_r2d.h
// If OUTPUT_TO_GMV is defined, material poly's resulting from splitting
// will be written to a gmv file.

TEST(split_order, ConvexPoly) {

   //Square
   std::vector<Tangram::Point2> square_pnts = {
   Tangram::Point2(0.0, 0.0), Tangram::Point2(1.0, 0.0),
   Tangram::Point2(1.0, 1.0), Tangram::Point2(0.0, 1.0)};
   
   Tangram::MatPoly<2> square_matpoly; 
   square_matpoly.initialize(square_pnts);
   std::vector<Tangram::MatPoly<2>> square_matpoly_vec{square_matpoly}; 

   //Cutting plane
   Tangram::Point2 plane_pt2d(0.5,0.0);
   Tangram::Vector2 normal2d(-0.5,0.0);
   Tangram::Plane_t<2> cutting_plane2d = get_cutting_plane<2>(plane_pt2d, normal2d);
   std::cout<<"\n----Using SplitR2D/ SplitR3D Classes----"<<std::endl;
   std::cout<<"----SQUARE----"<<std::endl;
   std::cout<<"    Cutting Plane: normal = {"<<cutting_plane2d.normal[0]
            <<", "<<cutting_plane2d.normal[1]<<"}"<<std::endl;
   std::cout<<"    Cutting Plane: dist2origin = "<<cutting_plane2d.dist2origin<<std::endl;
   
   //Split using SplitR2D class 
   Tangram::SplitR2D split2d(square_matpoly_vec, cutting_plane2d, true);
   Tangram::HalfSpaceSets_t<2> hs_sets2d = split2d();

   ///Print lower and upper subpolys
   ASSERT_EQ(hs_sets2d.lower_halfspace_set.matpolys.size(), 1);
   ASSERT_EQ(hs_sets2d.upper_halfspace_set.matpolys.size(), 1);
   std::cout<<"     Lower matpoly vertices"<<std::endl;
   for (int i = 0; i < hs_sets2d.lower_halfspace_set.matpolys[0].num_vertices(); i++)
  	std::cout<<"pt["<<i<<"] = "<<hs_sets2d.lower_halfspace_set.matpolys[0].vertex_point(i)<<std::endl;
   std::cout<<"     Upper matpoly vertices"<<std::endl;
   for (int i = 0; i < hs_sets2d.upper_halfspace_set.matpolys[0].num_vertices(); i++)
  	std::cout<<"pt["<<i<<"] = "<<hs_sets2d.upper_halfspace_set.matpolys[0].vertex_point(i)<<std::endl;

   //Cube
   std::vector<Tangram::Point3> cube_pnts = {
   Tangram::Point3(0.0, 0.0, 0.0), Tangram::Point3(1.0, 0.0, 0.0),
   Tangram::Point3(1.0, 1.0, 0.0), Tangram::Point3(0.0, 1.0, 0.0),
   Tangram::Point3(0.0, 0.0, 1.0), Tangram::Point3(1.0, 0.0, 1.0),
   Tangram::Point3(1.0, 1.0, 1.0), Tangram::Point3(0.0, 1.0, 1.0)};
 
   std::vector<std::vector<int>> cube_faces = {{0,3,2,1},{0,1,5,4},{1,2,6,5},
                                               {2,3,7,6},{0,4,7,3},{4,5,6,7}}; 
 
   Tangram::MatPoly<3> cube_matpoly; 
   cube_matpoly.initialize(cube_pnts, cube_faces);
   std::vector<Tangram::MatPoly<3>> cube_matpoly_vec{cube_matpoly}; 

   //Cutting plane
   Tangram::Point3 plane_pt3d(0.5, 0.0, 0.0);
   Tangram::Vector3 normal3d(-0.5, 0.0, 0.0);
   Tangram::Plane_t<3> cutting_plane3d = get_cutting_plane<3>(plane_pt3d, normal3d);
   std::cout<<"\n----CUBE----"<<std::endl;
   std::cout<<"    Cutting Plane: normal = {"<<cutting_plane3d.normal[0]
            <<", "<<cutting_plane3d.normal[1]<<", "<<cutting_plane3d.normal[2]<<"}"<<std::endl;
   std::cout<<"    Cutting Plane: dist2origin = "<<cutting_plane3d.dist2origin<<std::endl;
   
   //Split using SplitR2D class 
   Tangram::SplitR3D split3d(cube_matpoly_vec, cutting_plane3d, true);
   Tangram::HalfSpaceSets_t<3> hs_sets3d = split3d();

   ///Print lower and upper subpolys
   ASSERT_EQ(hs_sets3d.lower_halfspace_set.matpolys.size(), 1);
   ASSERT_EQ(hs_sets3d.upper_halfspace_set.matpolys.size(), 1);
   std::cout<<"     Lower matpoly vertices"<<std::endl;
   for (int i = 0; i < hs_sets3d.lower_halfspace_set.matpolys[0].num_vertices(); i++)
  	std::cout<<"pt["<<i<<"] = "<<hs_sets3d.lower_halfspace_set.matpolys[0].vertex_point(i)<<std::endl;
   std::cout<<"     Upper matpoly vertices"<<std::endl;
   for (int i = 0; i < hs_sets3d.upper_halfspace_set.matpolys[0].num_vertices(); i++)
  	std::cout<<"pt["<<i<<"] = "<<hs_sets3d.upper_halfspace_set.matpolys[0].vertex_point(i)<<std::endl;


}

TEST(split_order, r2dsplit) {

   //Square
   std::vector<Tangram::Point2> square_pnts = {
   Tangram::Point2(0.0, 0.0), Tangram::Point2(1.0, 0.0),
   Tangram::Point2(1.0, 1.0), Tangram::Point2(0.0, 1.0)};
   
   Tangram::MatPoly<2> square_matpoly; 
   square_matpoly.initialize(square_pnts);

   //Cutting plane
   Tangram::Point2 plane_pt2d(0.5,0.0);
   Tangram::Vector2 normal2d(-0.5,0.0);
   Tangram::Plane_t<2> cutting_plane2d = get_cutting_plane<2>(plane_pt2d, normal2d);
   std::cout<<"\n----Using split_convex_matpoly_r2d/r3d routines----"<<std::endl;
   std::cout<<"----SQUARE----"<<std::endl;
   std::cout<<"    Cutting Plane: normal = {"<<cutting_plane2d.normal[0]
            <<", "<<cutting_plane2d.normal[1]<<"}"<<std::endl;
   std::cout<<"    Cutting Plane: dist2origin = "<<cutting_plane2d.dist2origin<<std::endl;
   
   //Split using SplitR2D class 
   Tangram::MatPoly<2> square_low, square_up; 
   std::vector<double> square_lowmo, square_upmo;
   Tangram::split_convex_matpoly_r2d(square_matpoly, cutting_plane2d, square_low, 
                                     square_up, square_lowmo, square_upmo); 

   ///Print lower and upper subpolys
   std::cout<<"     Lower matpoly vertices"<<std::endl;
   for (int i = 0; i < square_low.num_vertices(); i++)
  	std::cout<<"pt["<<i<<"] = "<<square_low.vertex_point(i)<<std::endl;
   std::cout<<"     Upper matpoly vertices"<<std::endl;
   for (int i = 0; i < square_up.num_vertices(); i++)
  	std::cout<<"pt["<<i<<"] = "<<square_up.vertex_point(i)<<std::endl;

   //Cube
   std::vector<Tangram::Point3> cube_pnts = {
   Tangram::Point3(0.0, 0.0, 0.0), Tangram::Point3(1.0, 0.0, 0.0),
   Tangram::Point3(1.0, 1.0, 0.0), Tangram::Point3(0.0, 1.0, 0.0),
   Tangram::Point3(0.0, 0.0, 1.0), Tangram::Point3(1.0, 0.0, 1.0),
   Tangram::Point3(1.0, 1.0, 1.0), Tangram::Point3(0.0, 1.0, 1.0)};
 
   std::vector<std::vector<int>> cube_faces = {{0,3,2,1},{0,1,5,4},{1,2,6,5},
                                               {2,3,7,6},{0,4,7,3},{4,5,6,7}}; 
 
   Tangram::MatPoly<3> cube_matpoly; 
   cube_matpoly.initialize(cube_pnts, cube_faces);

   //Cutting plane
   Tangram::Point3 plane_pt3d(0.5, 0.0, 0.0);
   Tangram::Vector3 normal3d(-0.5, 0.0, 0.0);
   Tangram::Plane_t<3> cutting_plane3d = get_cutting_plane<3>(plane_pt3d, normal3d);
   std::cout<<"\n----CUBE----"<<std::endl;
   std::cout<<"    Cutting Plane: normal = {"<<cutting_plane3d.normal[0]
            <<", "<<cutting_plane3d.normal[1]<<", "<<cutting_plane3d.normal[2]<<"}"<<std::endl;
   std::cout<<"    Cutting Plane: dist2origin = "<<cutting_plane3d.dist2origin<<std::endl;
   
   //Split using SplitR2D class 
   Tangram::MatPoly<3> cube_low, cube_up; 
   std::vector<double> cube_lowmo, cube_upmo;
   Tangram::split_convex_matpoly_r3d(cube_matpoly, cutting_plane3d, cube_low, 
                                     cube_up, cube_lowmo, cube_upmo); 

   ///Print lower and upper subpolys
   std::cout<<"     Lower matpoly vertices"<<std::endl;
   for (int i = 0; i < cube_low.num_vertices(); i++)
  	std::cout<<"pt["<<i<<"] = "<<cube_low.vertex_point(i)<<std::endl;
   std::cout<<"     Upper matpoly vertices"<<std::endl;
   for (int i = 0; i < cube_up.num_vertices(); i++)
  	std::cout<<"pt["<<i<<"] = "<<cube_up.vertex_point(i)<<std::endl;
}
