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
  double dst_tol = sqrt(3)*std::numeric_limits<double>::epsilon();
  double vol_tol = pow(dst_tol, 3);

  //Right triangular prism
  std::vector<Tangram::Point3> prism_points = {
    Tangram::Point3(0.0, 0.0, 0.0), Tangram::Point3(1.0, 0.0, 0.0),
    Tangram::Point3(0.0, 1.0, 0.0), Tangram::Point3(0.0, 0.0, 1.0),
    Tangram::Point3(1.0, 0.0, 1.0), Tangram::Point3(0.0, 1.0, 1.0) };
  std::vector< std::vector<int> > prism_faces = {
    {0, 2, 1}, {0, 1, 4, 3}, {1, 2, 5, 4}, {2, 0, 3, 5}, {3, 4, 5} };

  Tangram::Plane_t<3> cutting_plane;
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
  prism_matpoly.initialize(prism_points, prism_faces, dst_tol);

  //Split the prism with the cutting plane
  Tangram::MatPoly<3> convex_polys[2];
  std::vector<double> cp_moments[2];
  split_convex_matpoly_r3d(prism_matpoly, cutting_plane, 
                           convex_polys[0], convex_polys[1], 
                           cp_moments[0], cp_moments[1], vol_tol, dst_tol);

  //Check that vertices and faces of polygons below and above the plane match
  for (int ihs = 0; ihs < 2; ihs++) {
    ASSERT_EQ(ref_cp_points[ihs].size(), convex_polys[ihs].num_vertices());
    std::vector<int> ref2res_vrt(ref_cp_points[ihs].size(), -1);
    for (int i = 0; i < ref_cp_points[ihs].size(); i++)
      for (int j = 0; j < convex_polys[ihs].num_vertices(); j++)
        if (Wonton::approxEq(ref_cp_points[ihs][i], convex_polys[ihs].vertex_point(j), dst_tol))
          ref2res_vrt[i] = j;
    for (int ivrt = 0; ivrt < ref_cp_points[ihs].size(); ivrt++)
      ASSERT_NE(ref2res_vrt[ivrt], -1);

    std::vector<int> ref2res_face(ref_cp_faces[ihs].size(), -1);
    ASSERT_EQ(ref_cp_faces[ihs].size(), convex_polys[ihs].num_faces());
    for (int i = 0; i < ref_cp_faces[ihs].size(); i++)
      for (int j = 0; j < convex_polys[ihs].num_faces(); j++) {
        const std::vector<int>& face_vrts = convex_polys[ihs].face_vertices(j);
        int nfvrts = static_cast<int>(face_vrts.size());
        if (ref_cp_faces[ihs][i].size() != nfvrts) continue;
        int iref_start = std::distance(face_vrts.begin(), 
          std::find(face_vrts.begin(), face_vrts.end(), ref2res_vrt[ref_cp_faces[ihs][i][0]]));
        if (iref_start == nfvrts) continue;

        bool matching_face = true;
        for (int iv = 0; iv < nfvrts; iv++)
          if (face_vrts[(iref_start + iv)%nfvrts] != ref2res_vrt[ref_cp_faces[ihs][i][iv]]) {
            matching_face = false;
            break;
          }

        if (matching_face) ref2res_face[i] = j;
      }
    
    for (int iface = 0; iface < ref_cp_faces[ihs].size(); iface++)
      ASSERT_NE(ref2res_face[iface], -1);
  }

  //Use the class instead
  std::vector< Tangram::MatPoly<3> > convex_matpolys = {prism_matpoly};
  Tangram::SplitR3D split_poly(convex_matpolys, cutting_plane, vol_tol, dst_tol, true);
  Tangram::HalfSpaceSets_t<3> hs_poly_sets = split_poly();

  ASSERT_EQ(hs_poly_sets.lower_halfspace_set.matpolys.size(), 1);
  ASSERT_EQ(hs_poly_sets.upper_halfspace_set.matpolys.size(), 1);                                          
  Tangram::MatPoly<3>* hs_poly_ptrs[2] = {&hs_poly_sets.lower_halfspace_set.matpolys[0],
                                          &hs_poly_sets.upper_halfspace_set.matpolys[0]};
  //Check that we obtained the same MatPoly's as before                                        
  for (int ihs = 0; ihs < 2; ihs++) {
    ASSERT_EQ(ref_cp_points[ihs].size(), hs_poly_ptrs[ihs]->num_vertices());
    std::vector<int> ref2res_vrt(ref_cp_points[ihs].size(), -1);
    for (int i = 0; i < ref_cp_points[ihs].size(); i++)
      for (int j = 0; j < hs_poly_ptrs[ihs]->num_vertices(); j++)
        if (Wonton::approxEq(ref_cp_points[ihs][i], hs_poly_ptrs[ihs]->vertex_point(j), dst_tol))
          ref2res_vrt[i] = j;
    for (int ivrt = 0; ivrt < ref_cp_points[ihs].size(); ivrt++)
      ASSERT_NE(ref2res_vrt[ivrt], -1);

    std::vector<int> ref2res_face(ref_cp_faces[ihs].size(), -1);
    ASSERT_EQ(ref_cp_faces[ihs].size(), hs_poly_ptrs[ihs]->num_faces());
    for (int i = 0; i < ref_cp_faces[ihs].size(); i++)
      for (int j = 0; j < hs_poly_ptrs[ihs]->num_faces(); j++) {
        const std::vector<int>& face_vrts = hs_poly_ptrs[ihs]->face_vertices(j);
        int nfvrts = static_cast<int>(face_vrts.size());
        if (ref_cp_faces[ihs][i].size() != nfvrts) continue;
        int iref_start = std::distance(face_vrts.begin(), 
          std::find(face_vrts.begin(), face_vrts.end(), ref2res_vrt[ref_cp_faces[ihs][i][0]]));
        if (iref_start == nfvrts) continue;

        bool matching_face = true;
        for (int iv = 0; iv < nfvrts; iv++)
          if (face_vrts[(iref_start + iv)%nfvrts] != ref2res_vrt[ref_cp_faces[ihs][i][iv]]) {
            matching_face = false;
            break;
          }

        if (matching_face) ref2res_face[i] = j;
      }
    
    for (int iface = 0; iface < ref_cp_faces[ihs].size(); iface++)
      ASSERT_NE(ref2res_face[iface], -1);
  }  

#ifdef OUTPUT_TO_GMV
  std::vector<std::shared_ptr<Tangram::CellMatPoly<3>>> cellmatpoly_list;
  cellmatpoly_list.push_back(std::make_shared< Tangram::CellMatPoly<3> >(0));
  for (int ihs = 0; ihs < 2; ihs++) {
    int nfaces = convex_polys[ihs].num_faces();
    std::vector<int> nface_vrts(nfaces);
    std::vector<int> faces_vrts;
    for (int iface = 0; iface < nfaces; iface++) {
      const std::vector<int> face_ivrts = convex_polys[ihs].face_vertices(iface);
      nface_vrts[iface] = face_ivrts.size();
      faces_vrts.insert(faces_vrts.end(), face_ivrts.begin(), face_ivrts.end());
    }
    cellmatpoly_list[0]->add_matpoly(ihs, convex_polys[ihs].num_vertices(), 
                                     &convex_polys[ihs].points()[0],
                                     nullptr, nullptr,
                                     convex_polys[ihs].num_faces(), 
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
  ncv_matpoly.initialize(ncv_poly_points, ncv_poly_faces, dst_tol);

  //Split non-convex MatPoly with the horizontal cutting plane
  std::vector< Tangram::MatPoly<3> > ncv_matpolys = {ncv_matpoly};
  Tangram::SplitR3D split_ncv_poly(ncv_matpolys, cutting_plane, vol_tol, dst_tol, false);
  hs_poly_sets = split_ncv_poly();

  std::vector< Tangram::MatPoly<3> >* hs_poly_sets_ptr[2] = {
    &hs_poly_sets.lower_halfspace_set.matpolys,
    &hs_poly_sets.upper_halfspace_set.matpolys };

  std::vector<double>* hs_moments_ptrs[2] = {
    &hs_poly_sets.lower_halfspace_set.moments, 
    &hs_poly_sets.upper_halfspace_set.moments }; 

  //Check the number of MatPoly below and above the plane
  ASSERT_EQ(hs_poly_sets_ptr[0]->size(), 14);
  ASSERT_EQ(hs_poly_sets_ptr[1]->size(), 8);

  //Check consistency between moments values from r3d and MatPoly
  std::vector< std::vector<double> > acc_moments(2);
  for (int ihs = 0; ihs < 2; ihs++) {
    acc_moments[ihs].resize(4, 0.0);
    for(int ipoly = 0; ipoly < hs_poly_sets_ptr[ihs]->size(); ipoly++) {
      std::vector<double> cur_moments = (*hs_poly_sets_ptr[ihs])[ipoly].moments();
      for(int im = 0; im < 4; im++)
        acc_moments[ihs][im] += cur_moments[im];
    }
    for(int im = 0; im < 4; im++)
      ASSERT_NEAR((*hs_moments_ptrs[ihs])[im], acc_moments[ihs][im], 1.0e-15);
  }

  //Confirm moments values for the part above the plane
  for(int im = 0; im < 4; im++)
    ASSERT_NEAR((*hs_moments_ptrs[1])[im], above_plane_moments[im], 1.0e-15);

  //Confirm moments values for the components above the plane  
  for(int ipoly = 0; ipoly < hs_poly_sets_ptr[1]->size(); ipoly++) {
    std::vector<double> cur_moments = (*hs_poly_sets_ptr[1])[ipoly].moments();
    for(int im = 0; im < 4; im++)
      ASSERT_NEAR(cur_moments[im], above_plane_components_moments[ipoly][im], 1.0e-15);
  }

  //Assign moments to non-convex poly
  Tangram::MatPoly<3> original_poly_copy = ncv_matpoly;
  const Tangram::MatPoly<3>& original_poly = ncv_matpoly;
  std::vector<double> combined_moments(4, 0.0);
  for (int ihs = 0; ihs < 2; ihs++)
    for(int im = 0; im < 4; im++)
      combined_moments[im] += (*hs_moments_ptrs[ihs])[im];
  original_poly.assign_moments(combined_moments);
  for(int im = 0; im < 4; im++)
    ASSERT_NEAR(combined_moments[im], original_poly_copy.moments()[im], 1.0e-15);

  //Decompose manually
  std::vector< Tangram::MatPoly<3> > cv_matpolys;
  ncv_matpoly.decompose(cv_matpolys);
  Tangram::SplitR3D split_cv_poly(cv_matpolys, cutting_plane, vol_tol, dst_tol, true);
  Tangram::HalfSpaceSets_t<3> alt_hs_poly_sets = split_cv_poly();

  std::vector< Tangram::MatPoly<3> >* alt_hs_poly_sets_ptr[2] = {
    &alt_hs_poly_sets.lower_halfspace_set.matpolys,
    &alt_hs_poly_sets.upper_halfspace_set.matpolys };

  //Check that we obtained the same MatPoly's as before                                        
  for (int ihs = 0; ihs < 2; ihs++) {                          
    ASSERT_EQ(hs_poly_sets_ptr[0]->size(), alt_hs_poly_sets_ptr[0]->size());
    for(int ipoly = 0; ipoly < hs_poly_sets_ptr[ihs]->size(); ipoly++) {
      ASSERT_EQ((*hs_poly_sets_ptr[ihs])[ipoly].num_vertices(), 
                (*alt_hs_poly_sets_ptr[ihs])[ipoly].num_vertices());
      for (int ivrt = 0; ivrt < (*hs_poly_sets_ptr[ihs])[ipoly].num_vertices(); ivrt++)
        ASSERT_TRUE(approxEq((*hs_poly_sets_ptr[ihs])[ipoly].vertex_point(ivrt),
                             (*alt_hs_poly_sets_ptr[ihs])[ipoly].vertex_point(ivrt), 1.0e-15));
      ASSERT_EQ((*hs_poly_sets_ptr[ihs])[ipoly].num_faces(), 
                (*alt_hs_poly_sets_ptr[ihs])[ipoly].num_faces()); 
      for (int iface = 0; iface < (*hs_poly_sets_ptr[ihs])[ipoly].num_faces(); iface++) {
        const std::vector<int>& face_vrts = 
          (*hs_poly_sets_ptr[ihs])[ipoly].face_vertices(iface);
        const std::vector<int>& alt_face_vrts = 
          (*alt_hs_poly_sets_ptr[ihs])[ipoly].face_vertices(iface);
        ASSERT_EQ(face_vrts.size(), alt_face_vrts.size());
        for (int ivrt = 0; ivrt < face_vrts.size(); ivrt++)
          ASSERT_EQ(face_vrts[ivrt], alt_face_vrts[ivrt]); 
      }
    }
  }  

  //Get moments using ClipR3D class and facetizing faces
  Tangram::ClipR3D clip_poly(dst_tol);
  clip_poly.set_matpolys(ncv_matpolys, false);
  clip_poly.set_plane(cutting_plane);

  std::vector<double> clipper_moments = clip_poly();
  ASSERT_EQ(clipper_moments.size(), 4);

  //Check that methods return the same values
  for(int im = 0; im < 4; im++)
    ASSERT_NEAR((*hs_moments_ptrs[0])[im], clipper_moments[im], 1.0e-15);

#ifdef OUTPUT_TO_GMV
  cellmatpoly_list.clear();
  cellmatpoly_list.push_back(std::make_shared< Tangram::CellMatPoly<3> >(0));
  for (int ihs = 0; ihs < 2; ihs++)
    for(int ipoly = 0; ipoly < sub_polys[ihs].size(); ipoly++) {
      int nfaces = sub_polys[ihs][ipoly].num_faces();
      std::vector<int> nface_vrts(nfaces);
      std::vector<int> faces_vrts;
      for (int iface = 0; iface < nfaces; iface++) {
        const std::vector<int> face_ivrts = sub_polys[ihs][ipoly].face_vertices(iface);
        nface_vrts[iface] = face_ivrts.size();
        faces_vrts.insert(faces_vrts.end(), face_ivrts.begin(), face_ivrts.end());
      }
      cellmatpoly_list[0]->add_matpoly(ihs, sub_polys[ihs][ipoly].num_vertices(), 
                                       &sub_polys[ihs][ipoly].points()[0],
                                       nullptr, nullptr,
                                       sub_polys[ihs][ipoly].num_faces(), 
                                       &nface_vrts[0], &faces_vrts[0],
                                       nullptr, nullptr);
    }
  Tangram::write_to_gmv(cellmatpoly_list, "ncv_poly_matpolys.gmv");
#endif  

}
