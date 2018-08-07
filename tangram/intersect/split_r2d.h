/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#ifndef INTERSECT_SPLIT_R2D_H
#define INTERSECT_SPLIT_R2D_H

#include <vector>
#include <algorithm>
#include <memory>
#include "tangram/support/tangram.h"
#include "tangram/support/MatPoly.h"

extern "C" {
#include "r2d.h"
}

/*!
 @file split_r2d.h
 @brief Routines for splitting convex or non-convex MatPoly's
 by a cutting plane. Uses R2D.
 A convex MatPoly will be split into two MatPoly's.
 A non-convex MatPoly will be decomposed into convex polyhedra using
 its centroid, several MatPoly's could be returned for each half-plane.
 All the faces of the input MatPoly are assumed to be planar and no
 additional facetization is performed, so when sequentially splitting of a
 faceted MatPoly, no additional facets will be introduced.
 */

namespace Tangram {

/*!
  @brief Converts a MatPoly to a polygon in R2D format.
  @param[in] mat_poly MatPoly object to convert
  @param[out] r2dpoly Corresponding R2D polygon
*/
void
matpoly_to_r2dpoly(const MatPoly<2>& mat_poly,
                   r2d_poly& r2dpoly) {

  //Translate coordinates of vertices to R2D format
  const std::vector<Point2>& matpoly_vrts = mat_poly.points();
  int nvrts = mat_poly.num_vertices();

  //Convert matpoly vertices to r2d_rvec2 array 
  r2d_rvec2* r2dized_poly_vrts = new r2d_rvec2[nvrts];

  for (int ivrt = 0; ivrt < nvrts; ivrt++)
    for (int ixy = 0; ixy < 2; ixy++)
      r2dized_poly_vrts[ivrt].xy[ixy] = matpoly_vrts[ivrt][ixy];

  //Initialize the r2dpoly 
  r2d_init_poly(&r2dpoly, r2dized_poly_vrts, nvrts);

  delete [] r2dized_poly_vrts;
}//matpoly_to_r2dpoly

/*!
  @brief Computes moments of the intersection of a MatPoly and an r2d_poly. 
  If MatPoly is convex, r2d_poly will be clipped with lines containing the
  faces of MatPoly. Otherwise, MatPoly will be decomposed into triangular
  MatPolys, each of which will be intersected with r2d_poly. 
  Note that r2d_poly does not need to be convex.
  @param[in] mat_poly MatPoly object to intersect with
  @param[in] r2dpoly r2d_poly that is intersected with MatPoly
  @param[out] intersection_moments Moments of the intersection
  @param[in] convex_matpoly flag indicating if MatPoly is convex: if not
  it will be decomposed into triangles
*/
void get_intersection_moments(const MatPoly<2>& mat_poly,
                              const r2d_poly& r2dpoly,
                              std::vector<double>& intersection_moments,
                              bool convex_matpoly = false) {   
  const int POLY_ORDER = 1;
  int nmoments = R2D_NUM_MOMENTS(POLY_ORDER);
  r2d_real r2d_moments[R2D_NUM_MOMENTS(POLY_ORDER)];

  if (convex_matpoly) {                                
    std::vector< Plane_t<2> > face_lines;
    mat_poly.face_planes(face_lines);
    r2d_poly intersection = r2dpoly;

    int nlines = static_cast<int>(face_lines.size());
    r2d_plane* r2d_face_lines = new r2d_plane [nlines];
    for (int iline = 0; iline < nlines; iline++) {
      for (int ixy = 0; ixy < 2; ixy++)
        r2d_face_lines[iline].n.xy[ixy] = -face_lines[iline].normal[ixy];
      r2d_face_lines[iline].d = -face_lines[iline].dist2origin;
    }

    r2d_clip(&intersection, r2d_face_lines, (r2d_int) nlines);
    
    delete [] r2d_face_lines;

    r2d_reduce(&intersection, r2d_moments, POLY_ORDER);
    intersection_moments.resize(nmoments);
    for (int im = 0; im < nmoments; im++)
      intersection_moments[im] = r2d_moments[im];
  }
  else {
    std::vector< MatPoly<2> > mat_poly_tris;
    mat_poly.decompose(mat_poly_tris);

    int ntris = static_cast<int>(mat_poly_tris.size());
    intersection_moments.assign(nmoments, 0.0);
    for (int itri = 0; itri < ntris; itri++) {
      std::vector< Plane_t<2> > face_lines;
      mat_poly_tris[itri].face_planes(face_lines);
      assert(face_lines.size() == 3);

      r2d_poly intersection = r2dpoly; 
      r2d_plane r2d_face_lines[3];
      for (int iline = 0; iline < 3; iline++) {
        for (int ixy = 0; ixy < 2; ixy++)
          r2d_face_lines[iline].n.xy[ixy] = -face_lines[iline].normal[ixy];
        r2d_face_lines[iline].d = -face_lines[iline].dist2origin;
      }

      r2d_clip(&intersection, r2d_face_lines, 3);

      r2d_reduce(&intersection, r2d_moments, POLY_ORDER);
      for (int im = 0; im < nmoments; im++)
        intersection_moments[im] += r2d_moments[im];
    }
  }
}

} // namespace Tangram

#endif