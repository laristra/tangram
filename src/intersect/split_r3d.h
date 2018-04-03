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
 @brief Routines for splitting a convex of a non-convex MatPoly
 by a cutting plane. Uses R3D.
 A convex MatPoly will be split into two MatPoly's.
 A non-convex MatPoly will be decomposed into convex polyhedra using
 its centroid, several MatPoly's will be returned for each half-plane.
 All the faces of the input MatPoly are assumed to be planar and no
 additional facetization is performed, so when sequentially splitting of a
 faceted MatPoly, no additional facets will be introduced.
 */

namespace Tangram {

/*!
  @brief Converts a MatPoly to a polyhedron in R3D format.
  @param[in] mat_poly MatPoly object to convert
  @param[out] r3dized_poly Corresponding R3D polyhedron
*/
void
r3dize_matpoly(const MatPoly<3>& mat_poly,
               r3d_poly& r3dized_poly) {
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
  
  r3d_init_poly(&r3dized_poly, r3dized_poly_vrts, nvrts, r3dized_poly_faces, nface_vrts, nfaces);

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
matpolize_r3dpoly(r3d_poly& r3dpoly,
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
  @brief Moments of MatPoly's components below the cutting plane.
  For a non-convex MatPoly no decomposition is performed.
  @param[in] mat_poly MatPoly to split
  @param[in] cutting_plane Cutting plane to split with
  @param[out] belowplane_moments Computed moments
*/
void
below_plane_moments(const MatPoly<3>& mat_poly,
                    const Plane_t& cutting_plane,
                    std::vector<double>& belowplane_moments) {
  //Translate the cutting plane to R3D format
  r3d_plane r3d_cut_plane;
  for (int ixyz = 0; ixyz < 3; ixyz++)
    r3d_cut_plane.n.xyz[ixyz] = -cutting_plane.normal[ixyz];
  r3d_cut_plane.d = -cutting_plane.dist2origin;

  r3d_poly r3dized_poly;
  r3dize_matpoly(mat_poly, r3dized_poly);
  r3d_clip(&r3dized_poly, &r3d_cut_plane, 1);

  const int POLY_ORDER = 1;
  r3d_real r3d_moments[R3D_NUM_MOMENTS(POLY_ORDER)];
  //Check if the subpoly is empty
  if (r3dized_poly.nverts == 0)
    belowplane_moments.clear();
  else {
    //Find the moments for a halfplane
    r3d_reduce(&r3dized_poly, r3d_moments, POLY_ORDER);
    belowplane_moments.assign(r3d_moments, r3d_moments + 4);
  }
}

/*!
  @brief Splits a convex MatPoly into two (convex) MatPoly's
  with a cutting plane.
  @param[in] mat_poly Convex MatPoly to split
  @param[in] cutting_plane Cutting plane to split with
  @param[out] sub_polys Vector of resulting MatPoly's with two elements: 
  sub_polys[0] is a MatPoly below the plane and 
  sub_polys[1] is a MatPoly above the plane
  @param[out] subpoly_moments Corresponding moments, 
  subpoly_moments[0] for sub_polys[0] and subpoly_moments[1] for sub_polys[1].
  subpoly_moments[0][0] is the volume of the MatPoly below the plane, and i'th
  coordinate of its centroid is subpoly_moments[0][i]/subpoly_moments[0][0].
*/
void
split_convex_matpoly(const MatPoly<3>& mat_poly,
                     const Plane_t& cutting_plane,
                     std::vector< MatPoly<3> >& sub_polys,
                     std::vector< std::vector<double> >& subpoly_moments) {
  //Translate the cutting plane to R3D format
  r3d_plane r3d_cut_plane;
  for (int ixyz = 0; ixyz < 3; ixyz++)
    r3d_cut_plane.n.xyz[ixyz] = cutting_plane.normal[ixyz];
  r3d_cut_plane.d = cutting_plane.dist2origin;

  r3d_poly r3dized_poly;
  r3d_poly r3d_subpolys[2];
  r3dize_matpoly(mat_poly, r3dized_poly);
  r3d_split(&r3dized_poly, 1, r3d_cut_plane, &r3d_subpolys[1], &r3d_subpolys[0]);

  sub_polys.resize(2);
  subpoly_moments.resize(2);

  const int POLY_ORDER = 1;
  r3d_real r3d_moments[R3D_NUM_MOMENTS(POLY_ORDER)];
  for (int isp = 0; isp < 2; isp++) {
    //Check if the subpoly is empty
    if (r3d_subpolys[isp].nverts == 0) {
      sub_polys[isp].clear();
      subpoly_moments[isp].clear();
      continue;
    }

    //Find the moments for a subpoly
    r3d_reduce(&r3d_subpolys[isp], r3d_moments, POLY_ORDER);
    subpoly_moments[isp].assign(r3d_moments, r3d_moments + 4);

    //Get a MatPoly for a subpoly
    std::vector< MatPoly<3> > sub_matpoly;
    matpolize_r3dpoly(r3d_subpolys[isp], sub_matpoly);
    int ncomponents = (int) sub_matpoly.size();
    if (ncomponents > 1) 
      throw std::runtime_error("Non-convex MatPoly is split using the method for convex MatPoly's!");
    if (ncomponents == 1) 
      sub_polys[isp] = sub_matpoly[0];
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
  @param[out] sub_polys Vector of two vectors of MatPoly's: 
  sub_polys[0] is a vector of MatPoly's below the plane and 
  sub_polys[1] is a vector of MatPoly's above the plane
  @param[out] halfplane_moments Aggregated moments of MatPoly's, 
  halfplane_moments[0] for all MatPoly's below the plane and 
  halfplane_moments[1] for all MatPoly's above the plane.
  halfplane_moments[0][0] is the total volume of materials below the plane, 
  and i'th coordinate of their centroid is 
  halfplane_moments[0][i]/halfplane_moments[0][0].
*/
void
split_matpoly(MatPoly<3>& mat_poly,
              const Plane_t& cutting_plane,
              std::vector< std::vector< MatPoly<3> > >& sub_polys,
              std::vector< std::vector<double> >& halfplane_moments) {
  sub_polys.clear();
  sub_polys.resize(2);
  halfplane_moments.clear();
  halfplane_moments.resize(2);
  for (int ihp = 0; ihp < 2; ihp++)
    halfplane_moments[ihp].resize(4, 0.0);

  std::vector< MatPoly<3> > convex_matpolys;
  decompose_matpoly(mat_poly, convex_matpolys);

  int ncpolys = (int) convex_matpolys.size();
  for (int icpoly = 0; icpoly < ncpolys; icpoly++) {
    std::vector< MatPoly<3> > cur_polys;
    std::vector< std::vector<double> > cur_moments;
    split_convex_matpoly(convex_matpolys[icpoly], cutting_plane, cur_polys, cur_moments);
    for (int ihp = 0; ihp < 2; ihp++)
      if (cur_polys[ihp].num_vertices() != 0) {
        sub_polys[ihp].push_back(cur_polys[ihp]);
        for (int im = 0; im < 4; im++)
          halfplane_moments[ihp][im] += cur_moments[ihp][im];
      }
  }
}

} // namespace Tangram

#endif