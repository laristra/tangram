/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#include "tangram/intersect/split_r3d.h"

#include "gtest/gtest.h"
#include "tangram/support/tangram.h"

//#define OUTPUT_TO_GMV 

#ifdef OUTPUT_TO_GMV
  #include "tangram/driver/CellMatPoly.h"
  #include "tangram/driver/write_to_gmv.h"
#endif

// Test methods in split_r3d.h
// If OUTPUT_TO_GMV is defined, material poly's resulting from splitting
// will be written to a gmv file.

TEST(split_r3d, Mesh3D) {
  //Right triangular prism
  std::vector<Tangram::Point3> prism_points = {
    Tangram::Point3(0.0, 0.0, 0.0), Tangram::Point3(1.0, 0.0, 0.0),
    Tangram::Point3(0.0, 1.0, 0.0), Tangram::Point3(0.0, 0.0, 1.0),
    Tangram::Point3(1.0, 0.0, 1.0), Tangram::Point3(0.0, 1.0, 1.0) };
  std::vector< std::vector<int> > prism_faces = {
    {0, 2, 1}, {0, 1, 4, 3}, {1, 2, 5, 4}, {2, 0, 3, 5}, {3, 4, 5} };

  Tangram::Plane_t cutting_plane;
  cutting_plane.normal = Tangram::Vector3(0.0, -1.0/sqrt(2.0), 1.0/sqrt(2.0));
  cutting_plane.dist2origin = -0.35355339059327373;    
  
  std::vector< std::vector<Tangram::Point3> > ref_cp_points(2);
  std::vector< std::vector< std::vector<int> > > ref_cp_faces(2);
  ref_cp_points[0] = {
    Tangram::Point3(0.0, 0.0, 0.0), Tangram::Point3(1.0, 0.0, 0.0),
    Tangram::Point3(0.0, 1.0, 0.0), Tangram::Point3(0.0, 1.0, 1.0),
    Tangram::Point3(0.0, 0.0, 0.5), Tangram::Point3(0.0, 0.5, 1.0),
    Tangram::Point3(1.0, 0.0, 0.5), Tangram::Point3(0.5, 0.5, 1.0) };
  ref_cp_faces[0] = {
    {0, 2, 1}, {0, 1, 6, 4}, {0, 4, 5, 3, 2}, {1, 2, 3, 7, 6}, {3, 5, 7}, {4, 6, 7, 5} };

  ref_cp_points[1] = {
    Tangram::Point3(0.0, 0.0, 1.0), Tangram::Point3(1.0, 0.0, 1.0),
    Tangram::Point3(0.0, 0.0, 0.5), Tangram::Point3(0.0, 0.5, 1.0),
    Tangram::Point3(1.0, 0.0, 0.5), Tangram::Point3(0.5, 0.5, 1.0) };
  ref_cp_faces[1] = {
    {0, 2, 4, 1}, {0, 1, 5, 3}, {0, 3, 2}, {1, 4, 5}, {2, 3, 5, 4} };

  //Initialize MatPoly
  Tangram::MatPoly<3> prism_matpoly;
  prism_matpoly.initialize(prism_points, prism_faces);

  //Split the prism with the cutting plane
  std::vector< Tangram::MatPoly<3> > convex_polys;
  std::vector< std::vector<double> > cp_moments;
  split_convex_matpoly(prism_matpoly, cutting_plane, convex_polys, cp_moments);

  //Check that vertices and faces of polygons below and above the plane match
  for (int ihp = 0; ihp < 2; ihp++) {
    ASSERT_EQ(ref_cp_points[ihp].size(), convex_polys[ihp].num_vertices());
    for (int ivrt = 0; ivrt < ref_cp_points[ihp].size(); ivrt++)
      ASSERT_TRUE(approxEq(ref_cp_points[ihp][ivrt], convex_polys[ihp].vertex_point(ivrt), 1.0e-15));
    ASSERT_EQ(ref_cp_faces[ihp].size(), convex_polys[ihp].num_faces()); 
    for (int iface = 0; iface < ref_cp_faces[ihp].size(); iface++) {
      const std::vector<int>& face_vrts = convex_polys[ihp].face_vertices(iface);
      ASSERT_EQ(ref_cp_faces[ihp][iface].size(), face_vrts.size());
      for (int ivrt = 0; ivrt < ref_cp_faces[ihp][iface].size(); ivrt++)
        ASSERT_EQ(ref_cp_faces[ihp][iface][ivrt], face_vrts[ivrt]); 
    }
  }

#ifdef OUTPUT_TO_GMV
  std::vector<std::shared_ptr<Tangram::CellMatPoly<3>>> cellmatpoly_list;
  cellmatpoly_list.push_back(std::make_shared< Tangram::CellMatPoly<3> >(0));
  for (int ihp = 0; ihp < 2; ihp++) {
    int nfaces = convex_polys[ihp].num_faces();
    std::vector<int> nface_vrts(nfaces);
    std::vector<int> faces_vrts;
    for (int iface = 0; iface < nfaces; iface++) {
      const std::vector<int> face_ivrts = convex_polys[ihp].face_vertices(iface);
      nface_vrts[iface] = face_ivrts.size();
      faces_vrts.insert(faces_vrts.end(), face_ivrts.begin(), face_ivrts.end());
    }
    cellmatpoly_list[0]->add_matpoly(ihp, convex_polys[ihp].num_vertices(), 
                                     &convex_polys[ihp].points()[0],
                                     nullptr, nullptr,
                                     convex_polys[ihp].num_faces(), 
                                     &nface_vrts[0], &faces_vrts[0],
                                     nullptr, nullptr);
  }
  Tangram::write_to_gmv(cellmatpoly_list, "right_prism_matpolys.gmv");
#endif

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

  //Tangram::Plane_t cutting_plane;
  cutting_plane.normal = Tangram::Vector3(0.0, 0.0, 1.0);
  cutting_plane.dist2origin = -0.9;

  std::vector<double> above_plane_moments = {0.034666666666666651, 
    0.017333333333333326, 0.017333333333333326, 0.032399999999999984};

  std::vector< std::vector<double> > above_plane_components_moments = {
    {0.0046666666666666653, 0.0023333333333333322, 0.00014999999999999993, 0.0043583333333333312},
    {0.0046666666666666653, 0.0045166666666666645, 0.0023333333333333318, 0.0043583333333333312},
    {0.0046666666666666653, 0.0023333333333333322, 0.0045166666666666645, 0.0043583333333333312},
    {0.0046666666666666653, 0.00014999999999999993, 0.0023333333333333331, 0.0043583333333333312},
    {0.0039999999999999983, 0.0019999999999999992, 0.00038333333333333307, 0.003741666666666664},
    {0.0039999999999999983, 0.0036166666666666656, 0.0019999999999999987, 0.0037416666666666653},
    {0.0039999999999999983, 0.0019999999999999992, 0.0036166666666666656, 0.0037416666666666653},
    {0.0039999999999999992, 0.00038333333333333313, 0.0019999999999999996, 0.0037416666666666644} };  

  //Initialize non-convex MatPoly
  Tangram::MatPoly<3> ncv_matpoly;
  ncv_matpoly.initialize(ncv_poly_points, ncv_poly_faces);

  //Split non-convex MatPoly with the horizontal cutting plane
  std::vector< std::vector< Tangram::MatPoly<3> > > sub_polys;
  std::vector< std::vector<double> > halfplane_moments;
  Tangram::split_matpoly(ncv_matpoly, cutting_plane, sub_polys, halfplane_moments);

  //Check the number of MatPoly below and above the plane
  ASSERT_EQ(sub_polys[0].size(), 14);
  ASSERT_EQ(sub_polys[1].size(), 8);

  //Check consistency between moments values from r3d and MatPoly
  std::vector< std::vector<double> > acc_moments(2);
  for (int ihp = 0; ihp < 2; ihp++) {
    acc_moments[ihp].resize(4, 0.0);
    for(int ipoly = 0; ipoly < sub_polys[ihp].size(); ipoly++) {
      std::vector<double> cur_moments = sub_polys[ihp][ipoly].moments();
      for(int im = 0; im < 4; im++)
        acc_moments[ihp][im] += cur_moments[im];
    }
    for(int im = 0; im < 4; im++)
      ASSERT_NEAR(halfplane_moments[ihp][im], acc_moments[ihp][im], 1.0e-15);
  }

  //Confirm moments values for the part above the plane
  for(int im = 0; im < 4; im++)
    ASSERT_NEAR(halfplane_moments[1][im], above_plane_moments[im], 1.0e-15);

  //Confirm moments values for the components above the plane  
  for(int ipoly = 0; ipoly < sub_polys[1].size(); ipoly++) {
    std::vector<double> cur_moments = sub_polys[1][ipoly].moments();
    for(int im = 0; im < 4; im++)
      ASSERT_NEAR(cur_moments[im], above_plane_components_moments[ipoly][im], 1.0e-15);
  }

  //Get moments without decomposing into convex components
  std::vector<double> fast_below_plane_moments;
  Tangram::below_plane_moments(ncv_matpoly, cutting_plane, fast_below_plane_moments);
  ASSERT_EQ(fast_below_plane_moments.size(), 4);
  //Check that both methods return the same values
  for(int im = 0; im < 4; im++)
    ASSERT_NEAR(halfplane_moments[0][im], fast_below_plane_moments[im], 1.0e-15);

#ifdef OUTPUT_TO_GMV
  cellmatpoly_list.clear();
  cellmatpoly_list.push_back(std::make_shared< Tangram::CellMatPoly<3> >(0));
  for (int ihp = 0; ihp < 2; ihp++)
    for(int ipoly = 0; ipoly < sub_polys[ihp].size(); ipoly++) {
      int nfaces = sub_polys[ihp][ipoly].num_faces();
      std::vector<int> nface_vrts(nfaces);
      std::vector<int> faces_vrts;
      for (int iface = 0; iface < nfaces; iface++) {
        const std::vector<int> face_ivrts = sub_polys[ihp][ipoly].face_vertices(iface);
        nface_vrts[iface] = face_ivrts.size();
        faces_vrts.insert(faces_vrts.end(), face_ivrts.begin(), face_ivrts.end());
      }
      cellmatpoly_list[0]->add_matpoly(ihp, sub_polys[ihp][ipoly].num_vertices(), 
                                       &sub_polys[ihp][ipoly].points()[0],
                                       nullptr, nullptr,
                                       sub_polys[ihp][ipoly].num_faces(), 
                                       &nface_vrts[0], &faces_vrts[0],
                                       nullptr, nullptr);
    }
  Tangram::write_to_gmv(cellmatpoly_list, "ncv_poly_matpolys.gmv");
#endif  

}