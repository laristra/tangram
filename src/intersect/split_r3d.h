/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#ifndef SPLIT_R3D_H
#define SPLIT_R3D_H

#include <vector>
#include <algorithm>
#include <numeric>
#include <memory>
#include "tangram/support/tangram.h"
#include "tangram/support/MatPoly.h"

extern "C" {
#include "r3d.h"
}

/*!
 @file split_r3d.h
 @brief Routines for splitting convex or non-convex MatPoly's
 by a cutting plane. Uses R3D.
 A convex MatPoly will be split into two MatPoly's.
 A non-convex MatPoly will be decomposed into convex polyhedra using
 its centroid, several MatPoly's could be returned for each half-plane.
 All the faces of the input MatPoly are assumed to be planar and no
 additional facetization is performed, so when sequentially splitting of a
 faceted MatPoly, no additional facets will be introduced.
 */

namespace Tangram {

/*!
  @brief Converts a MatPoly to a polyhedron in R3D format.
  @param[in] mat_poly MatPoly object to convert
  @param[out] r3dpoly Corresponding R3D polyhedron
*/
void
matpoly_to_r3dpoly(const MatPoly<3>& mat_poly,
                   r3d_poly& r3dpoly) {
  //Translate coordinates of vertices to R3D format
  const std::vector<Point3>& matpoly_vrts = mat_poly.points();
  int nvrts = mat_poly.num_vertices();
  r3d_rvec3* r3dized_poly_vrts;
  r3dized_poly_vrts = new r3d_rvec3 [nvrts];
  for (int ivrt = 0; ivrt < nvrts; ivrt++)
    for (int ixyz = 0; ixyz < 3; ixyz++)
      r3dized_poly_vrts[ivrt].xyz[ixyz] = matpoly_vrts[ivrt][ixyz];
  
  //Get faces of the polyhedron in terms of local indices of its vertices
  int nfaces = mat_poly.num_faces();
  r3d_int* nface_vrts;
  nface_vrts = new r3d_int [nfaces];
  r3d_int** r3dized_poly_faces;
  r3dized_poly_faces = new r3d_int* [nfaces];
  for (int iface = 0; iface < nfaces; iface++) {
    const std::vector<int>& matpoly_faces = mat_poly.face_vertices(iface);
    int face_nvrts = (int) matpoly_faces.size();
    nface_vrts[iface] = face_nvrts;
    r3dized_poly_faces[iface] = new r3d_int [face_nvrts];
    for (int ivrt = 0; ivrt < face_nvrts; ivrt++)
      r3dized_poly_faces[iface][ivrt] = matpoly_faces[ivrt];
  }
  
  r3d_init_poly(&r3dpoly, r3dized_poly_vrts, nvrts, r3dized_poly_faces, nface_vrts, nfaces);

  delete [] r3dized_poly_vrts;
  delete [] nface_vrts;
  for (int iface = 0; iface < nfaces; iface++)
    delete [] r3dized_poly_faces[iface];
  delete [] r3dized_poly_faces;
}

/*!
  @brief Converts a polyhedron in R3D format to MatPoly's.
  @param[in] r3dpoly R3D polyhedron to convert
  @param[out] mat_polys Corresponding vector of MatPoly's: 
  will contain as many MatPoly's as there are components 
  in the R3D polyhedron
*/
void
r3dpoly_to_matpolys(r3d_poly& r3dpoly,
                    std::vector< MatPoly<3> >& mat_polys) {
  r3d_brep* poly_brep;
  r3d_int ncomponents;
  r3d_init_brep(&r3dpoly, &poly_brep, &ncomponents);

  mat_polys.clear();
  if (ncomponents == 0) {
    r3d_free_brep(&poly_brep, 0);
    return;
  }

  mat_polys.resize(ncomponents);
  for (int ipoly = 0; ipoly < ncomponents; ipoly++) {
    int nvrts = poly_brep[ipoly].numvertices;
    std::vector<Point3> curpoly_vrts;
    curpoly_vrts.reserve(nvrts);
    for (int ivrt = 0; ivrt < nvrts; ivrt++) {
      Point3 cur_vrt;
      for (int ixyz = 0; ixyz < 3; ixyz++)
        cur_vrt[ixyz] = poly_brep[ipoly].vertices[ivrt].xyz[ixyz];
      curpoly_vrts.push_back(cur_vrt);
    }

    int nfaces = poly_brep[ipoly].numfaces;
    std::vector< std::vector<int> > curpoly_face(nfaces);
    for (int iface = 0; iface < nfaces; iface++) {
      int face_nverts = poly_brep[ipoly].numvertsperface[iface];
      curpoly_face[iface].assign(poly_brep[ipoly].faceinds[iface],
                                 poly_brep[ipoly].faceinds[iface] + face_nverts);
    }

    mat_polys[ipoly].initialize(curpoly_vrts, curpoly_face);
  }

  r3d_free_brep(&poly_brep, ncomponents);
}

/*!
  @brief Splits a convex MatPoly into two (convex) MatPoly's
  with a cutting plane.
  @param[in] mat_poly Convex MatPoly to split
  @param[in] cutting_plane Cutting plane to split with
  @param[out] lower_halfspace_poly MatPoly below the plane
  @param[out] upper_halfspace_poly MatPoly above the plane
  @param[out] lower_halfspace_moments Moments of MatPoly below the plane
  @param[out] upper_halfspace_moments Moments of MatPoly above the plane
*/
void
split_convex_matpoly_r3d(const MatPoly<3>& mat_poly,
                         const Plane_t& cutting_plane,
                         MatPoly<3>& lower_halfspace_poly,
                         MatPoly<3>& upper_halfspace_poly,
                         std::vector<double>& lower_halfspace_moments,
                         std::vector<double>& upper_halfspace_moments) {
  //Translate the cutting plane to R3D format
  r3d_plane r3d_cut_plane;
  for (int ixyz = 0; ixyz < 3; ixyz++)
    r3d_cut_plane.n.xyz[ixyz] = cutting_plane.normal[ixyz];
  r3d_cut_plane.d = cutting_plane.dist2origin;

  r3d_poly r3dized_poly;
  r3d_poly r3d_subpolys[2];
  matpoly_to_r3dpoly(mat_poly, r3dized_poly);
  r3d_split(&r3dized_poly, 1, r3d_cut_plane, &r3d_subpolys[1], &r3d_subpolys[0]);

  MatPoly<3>* subpoly_ptrs[2] = {&lower_halfspace_poly, &upper_halfspace_poly};
  std::vector<double>* subpoly_moments_ptrs[2] = {&lower_halfspace_moments, 
                                                  &upper_halfspace_moments};

  const int POLY_ORDER = 1;
  r3d_real r3d_moments[R3D_NUM_MOMENTS(POLY_ORDER)];
  for (int isp = 0; isp < 2; isp++) {
    //Check if the subpoly is empty
    if (r3d_subpolys[isp].nverts == 0) {
      subpoly_ptrs[isp]->clear();
      subpoly_moments_ptrs[isp]->clear();
      continue;
    }

    //Find the moments for a subpoly
    r3d_reduce(&r3d_subpolys[isp], r3d_moments, POLY_ORDER);
    subpoly_moments_ptrs[isp]->assign(r3d_moments, r3d_moments + 4);

    //Get a MatPoly for a subpoly
    std::vector< MatPoly<3> > sub_matpoly;
    r3dpoly_to_matpolys(r3d_subpolys[isp], sub_matpoly);
    int ncomponents = (int) sub_matpoly.size();
    if (ncomponents > 1) 
      throw std::runtime_error("Non-convex MatPoly is split using the method for convex MatPoly's!");
    if (ncomponents == 1) 
      *subpoly_ptrs[isp] = sub_matpoly[0];
  }  
}

/*!
  @brief Decomposes a MatPoly into convex MatPoly's using its centroid.
  @param[in] mat_poly MatPoly to decompose
  @param[out] convex_matpolys Corresponding vector of MatPoly's: 
  will contain as many MatPoly's as mat_poly has faces
*/
void
decompose_matpoly(MatPoly<3>& mat_poly,
                  std::vector< MatPoly<3> >& convex_matpolys) {
  std::vector<double> matpoly_moments = mat_poly.moments();
  Point3 matpoly_cen;
  for (int ixyz = 0; ixyz < 3; ixyz++)
    matpoly_cen[ixyz] = matpoly_moments[ixyz + 1]/matpoly_moments[0];

  int nfaces = mat_poly.num_faces();
  convex_matpolys.clear();
  convex_matpolys.resize(nfaces);

  for (int icpoly = 0; icpoly < nfaces; icpoly++) {
    const std::vector<int>& face_vrts = mat_poly.face_vertices(icpoly);
    int face_nvrts = (int) face_vrts.size();

    std::vector<Point3> convex_matpoly_vrts(face_nvrts + 1);
    std::vector< std::vector<int> > convex_matpoly_faces(face_nvrts + 1);
    for (int ivrt = 0; ivrt < face_nvrts; ivrt++) {
      convex_matpoly_vrts[ivrt] = mat_poly.vertex_point(face_vrts[ivrt]);
      convex_matpoly_faces[ivrt] = {face_nvrts, (ivrt + 1)%face_nvrts, ivrt};
    }
    convex_matpoly_vrts[face_nvrts] = matpoly_cen;
    convex_matpoly_faces[face_nvrts].resize(face_nvrts);
    std::iota(convex_matpoly_faces[face_nvrts].begin(), convex_matpoly_faces[face_nvrts].end(), 0);      
    
    convex_matpolys[icpoly].initialize(convex_matpoly_vrts, convex_matpoly_faces);
  }
}

/*!
  @brief Splits a non-convex MatPoly into two collections of convex MatPoly's
  with a cutting plane.
  @param[in] mat_poly MatPoly to split
  @param[in] cutting_plane Cutting plane to split with
  @param[out] lower_halfspace_polys Vector of MatPoly's below the plane
  @param[out] upper_halfspace_polys Vector of MatPoly's above the plane
  @param[out] lower_halfspace_moments Aggregated moments of all MatPoly's below the plane
  @param[out] upper_halfspace_moments Aggregated moments of all MatPoly's above the plane
*/
void
split_nonconvex_matpoly_r3d(MatPoly<3>& mat_poly,
                            const Plane_t& cutting_plane,
                            std::vector< MatPoly<3> >& lower_halfspace_polys,
                            std::vector< MatPoly<3> >& upper_halfspace_polys,
                            std::vector<double>& lower_halfspace_moments,
                            std::vector<double>& upper_halfspace_moments) {

  std::vector< MatPoly<3> >* halfspace_subpoly_ptrs[2] = {&lower_halfspace_polys, 
                                                          &upper_halfspace_polys};
  std::vector<double>* halfspace_moments_ptrs[2] = {&lower_halfspace_moments, 
                                                    &upper_halfspace_moments};                              
  for (int ihs = 0; ihs < 2; ihs++)
    halfspace_moments_ptrs[ihs]->assign(4, 0.0);

  std::vector< MatPoly<3> > convex_matpolys;
  decompose_matpoly(mat_poly, convex_matpolys);

  int ncpolys = (int) convex_matpolys.size();
  for (int icpoly = 0; icpoly < ncpolys; icpoly++) {
    MatPoly<3> cur_polys[2];
    std::vector<double> cur_moments[2];
    split_convex_matpoly_r3d(convex_matpolys[icpoly], cutting_plane, 
                             cur_polys[0], cur_polys[1],
                             cur_moments[0], cur_moments[1]);
    for (int ihs = 0; ihs < 2; ihs++)
      if (cur_polys[ihs].num_vertices() != 0) {
        halfspace_subpoly_ptrs[ihs]->emplace_back(cur_polys[ihs]);
        for (int im = 0; im < 4; im++)
          (*halfspace_moments_ptrs[ihs])[im] += cur_moments[ihs][im];
      }
  }
}

/*!
  @brief Moments of MatPoly's components below the cutting plane.
  For a non-convex MatPoly no decomposition is performed.
  @param[in] mat_poly MatPoly to split
  @param[in] cutting_plane Cutting plane to split with
  @param[out] lower_halfspace_moments Computed moments
*/
void
lower_halfspace_moments_r3d(const MatPoly<3>& mat_poly,
                            const Plane_t& cutting_plane,
                            std::vector<double>& lower_halfspace_moments) {
  //Translate the cutting plane to R3D format
  r3d_plane r3d_cut_plane;
  for (int ixyz = 0; ixyz < 3; ixyz++)
    r3d_cut_plane.n.xyz[ixyz] = -cutting_plane.normal[ixyz];
  r3d_cut_plane.d = -cutting_plane.dist2origin;

  r3d_poly r3dized_poly;
  matpoly_to_r3dpoly(mat_poly, r3dized_poly);
  r3d_clip(&r3dized_poly, &r3d_cut_plane, 1);

  const int POLY_ORDER = 1;
  r3d_real r3d_moments[R3D_NUM_MOMENTS(POLY_ORDER)];
  //Check if the subpoly is empty
  if (r3dized_poly.nverts == 0)
    lower_halfspace_moments.clear();
  else {
    //Find the moments for a halfspace
    r3d_reduce(&r3dized_poly, r3d_moments, POLY_ORDER);
    lower_halfspace_moments.assign(r3d_moments, r3d_moments + 4);
  }
}

/*!
 * \class ChopR3D  3-D chopping algorithm
 *
 * In this routine, we will utilize the fact that R3D can clip a 
 * non-convex polyhedron with a plane. We are given a target material
 * polyhedron with planar faces and a plane to split it. By default, 
 * material polyhedron is assumed to have non-planar faces which 
 * are facetized.
 *
 * The moments of the part of the polyhedron that is below the plane
 * are computed. This can be interpreted as moments of a chopped off part
 * in the nested dissections algorithm.
 */

class ChopR3D {
 public:
  ChopR3D(const std::vector< MatPoly<3> >& mat_polys,
          const Plane_t& cutting_plane, 
          const bool planar_faces = false) : 
          mat_polys_(mat_polys), 
          cutting_plane_(cutting_plane),
          planar_faces_(planar_faces) {}

  /*! \brief Computes moments of a chopped off poly: the part of MatPoly
  below the plane
   * \param[in] matpoly_ID ID of a MatPoly to chop
   * \return Weights_t structure containing moments of chopped off part
   */
  Weights_t operator() (const int matpoly_ID) const {
    Weights_t curpoly_weights;
    curpoly_weights.entityID = matpoly_ID;

    if (planar_faces_)
      lower_halfspace_moments_r3d(mat_polys_[matpoly_ID], cutting_plane_, 
                                  curpoly_weights.weights);
    else {
      MatPoly<3> faceted_poly;
      mat_polys_[matpoly_ID].faceted_matpoly(&faceted_poly);
      lower_halfspace_moments_r3d(faceted_poly, cutting_plane_, 
                                  curpoly_weights.weights);
    }
    
    return curpoly_weights;
  }

//! Default constructor (disabled)
  ChopR3D() = delete;

  //! Assignment operator (disabled)
  ChopR3D& operator = (const ChopR3D&) = delete;

 private:
  const std::vector< MatPoly<3> >& mat_polys_;
  const Plane_t& cutting_plane_;
  bool planar_faces_;
};

} // namespace Tangram

#endif