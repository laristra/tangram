/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#ifndef INTERSECT_SPLIT_R2D_H
#define INTERSECT_SPLIT_R2D_H

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

#include "tangram/support/tangram.h"
#include "tangram/support/MatPoly.h"

extern "C" {
#include "wonton/intersect/r3d/r2d.h"
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
  @brief Converts a R2D poly to a MatPoly where 
  the vertices in the Matpoly are ordered counter 
  clock-wise. 
  @param[in] r2dpoly Corresponding R2D polygon
  @param[out] mat_poly MatPoly object to convert
  @param[in] dst_tol Distance tolerance
  @param[in] reference_pts Preferred points to snap vertices to
*/
void r2dpoly_to_matpoly(const r2d_poly& r2dpoly, MatPoly<2>& mat_poly,
                        const double dst_tol,
                        const std::vector< Point<2> >* reference_pts = nullptr)
{
    int nvrts = r2dpoly.nverts;
    std::vector< Point<2> > r2d_poly_vrts(nvrts);

    // Walk the r2d graph to collect vertices in the ccw order
    int icur_node = 0;
    for (int ivrt = 0; ivrt < nvrts; ivrt++) {
      for (int ixy = 0; ixy < 2; ixy++)
        r2d_poly_vrts[ivrt][ixy] = r2dpoly.verts[icur_node].pos.xy[ixy];

      icur_node = r2dpoly.verts[icur_node].pnbrs[0];
    }

    mat_poly = natural_selection(r2d_poly_vrts, dst_tol, reference_pts);
}

/*!
  @brief Converts a MatPoly to a polygon in R2D format.
  @param[in] mat_poly MatPoly object to convert
  @param[out] r2dpoly Corresponding R2D polygon
*/
 void matpoly_to_r2dpoly(const MatPoly<2>& mat_poly,
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
                                bool convex_matpoly) {   
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
void split_convex_matpoly_r2d(const MatPoly<2>& mat_poly,
                              const Plane_t<2>& cutting_plane,
                              MatPoly<2>& lower_halfspace_poly,
                              MatPoly<2>& upper_halfspace_poly,
                              std::vector<double>& lower_halfspace_moments,
                              std::vector<double>& upper_halfspace_moments,
                              double vol_tol, double dst_tol) {

  //Translate the cutting plane to R2D format
  //Here the direction of the cutting plane normal
  //and distance are inverted compared to the same
  //for 3D split. The notions of what is lower and
  //upper halfspace in inverted between R2D and R3D. 
  r2d_plane r2d_cut_plane;
  for (int ixy = 0; ixy < 2; ixy++)
    r2d_cut_plane.n.xy[ixy] = -cutting_plane.normal[ixy];
  r2d_cut_plane.d = -cutting_plane.dist2origin;

  //Convert matpoly to r2d_poly and split 
  r2d_poly r2dized_poly;
  r2d_poly r2d_subpolys[2];
  matpoly_to_r2dpoly(mat_poly, r2dized_poly);
  r2d_split(&r2dized_poly, 1, r2d_cut_plane, &r2d_subpolys[1], &r2d_subpolys[0]);

  MatPoly<2>* subpoly_ptrs[2] = {&lower_halfspace_poly, &upper_halfspace_poly};
  std::vector<double>* subpoly_moments_ptrs[2] = {&lower_halfspace_moments,
                                                  &upper_halfspace_moments};

  //Compute moments for subpoly's 
  const int POLY_ORDER = 1;
  r2d_real r2d_moments[R2D_NUM_MOMENTS(POLY_ORDER)];
  int iempty_subpoly = -1;
  for (int isp = 0; isp < 2; isp++) {
    int nverts = r2d_subpolys[isp].nverts;

    //Check if the subpoly is empty
    if (nverts == 0) {
	    subpoly_ptrs[isp]->clear();
	    subpoly_moments_ptrs[isp]->assign(3, 0.0);
      iempty_subpoly = isp;
	    break;
    }

    //Find the moments for a subpoly
    r2d_reduce(&r2d_subpolys[isp], r2d_moments, POLY_ORDER);

    if (r2d_moments[0] < vol_tol) {
	    subpoly_ptrs[isp]->clear();
	    subpoly_moments_ptrs[isp]->assign(3, 0.0);
      iempty_subpoly = isp;
	    break;
    }

    subpoly_moments_ptrs[isp]->assign(r2d_moments, r2d_moments + 3);
  }

  if (iempty_subpoly == -1) {
    for (int isp = 0; isp < 2; isp++) {
      //Get a MatPoly for a r2d subpoly
      //We snap vertices to original positions if they are close enough
      MatPoly<2> sub_matpoly; 
      r2dpoly_to_matpoly(r2d_subpolys[isp], sub_matpoly, dst_tol, &mat_poly.points());
      
      *subpoly_ptrs[isp] = sub_matpoly; 
    }
  }
  else {
    int ifull_poly = (iempty_subpoly + 1)%2;
    matpoly_to_r2dpoly(mat_poly, r2d_subpolys[ifull_poly]);
    r2d_reduce(&r2d_subpolys[ifull_poly], r2d_moments, POLY_ORDER);
    subpoly_moments_ptrs[ifull_poly]->assign(r2d_moments, r2d_moments + 3);

    *subpoly_ptrs[ifull_poly] = mat_poly;
  }
} //split_convex_matpoly_r2d


/*!
 * \class SplitR2D  2-D splitting algorithm
 *
 * Splits a vector of MatPoly's with a cutting plane. The result is HalfSpaceSets_t
 * structure that contains two sets of convex MatPoly's: one corresponds to components below,
 * and the other to components above the plane. Each set has associated moments, which
 * are aggregated moments of all components in the set.
 * 
 * Important: input MatPoly's are assumed to have planar faces.
 * By default, input MatPoly's are assumed to be non-covex and are decomposed into
 * convex components. If all_convex is set to true in the constructor, all input
 * MatPoly's are assumed to be convex and will NOT be decomposed.
 */

class SplitR2D {
 public:
  SplitR2D(const std::vector< MatPoly<2> >& matpolys,
           const Plane_t<2>& cutting_plane,
           const double vol_tol,
           const double dst_tol, 
           const bool all_convex) : 
           matpolys_(matpolys), 
           cutting_plane_(cutting_plane),
           vol_tol_(vol_tol),
           dst_tol_(dst_tol),
           all_convex_(all_convex) {}

  /*! 
    @brief Splits a MatPoly into two sets of convex MatPoly's
    @return Structure containing two sets of MatPoly's, one for each half-space,
    with respective aggregated moments.
  */
  HalfSpaceSets_t<2> operator() () const {
    HalfSpaceSets_t<2> hs_sets;

    std::vector< MatPoly<2> >* hs_subpolys_ptrs[2] = {
      &hs_sets.lower_halfspace_set.matpolys, &hs_sets.upper_halfspace_set.matpolys};
    std::vector<double>* hs_moments_ptrs[2] = {
      &hs_sets.lower_halfspace_set.moments, &hs_sets.upper_halfspace_set.moments}; 

    std::vector< MatPoly<2> > convex_components;
    const std::vector< MatPoly<2> >* convex_polys;
    if (all_convex_)
      convex_polys = &matpolys_;
    else {
      int nncpolys = (int) matpolys_.size();
      for (int incp = 0; incp < nncpolys; incp++)
        matpolys_[incp].decompose(convex_components);

      convex_polys = &convex_components;
    }

    for (int ihs = 0; ihs < 2; ihs++)
      hs_moments_ptrs[ihs]->assign(3, 0.0);

    int hs_poly_count[2] = {0, 0};
    int npolys = (int) convex_polys->size();

    for (int icp = 0; icp < npolys; icp++) {
      MatPoly<2> cur_subpolys[2];
      std::vector<double> cur_moments[2];
      split_convex_matpoly_r2d((*convex_polys)[icp], cutting_plane_,
                               cur_subpolys[0], cur_subpolys[1],
                               cur_moments[0], cur_moments[1],
                               vol_tol_, dst_tol_);

      for (int ihs = 0; ihs < 2; ihs++)
        if (cur_subpolys[ihs].num_vertices() != 0) {
          hs_subpolys_ptrs[ihs]->emplace_back(cur_subpolys[ihs]);
          for (int im = 0; im < 3; im++)
            (*hs_moments_ptrs[ihs])[im] += cur_moments[ihs][im];
          hs_poly_count[ihs]++;
        }
    }

    for (int ihs = 0; ihs < 2; ihs++)
      if (hs_poly_count[ihs] == 0)
        hs_moments_ptrs[ihs]->assign(3, 0.0);

    return hs_sets;
  }

  //! Default constructor (disabled)
  SplitR2D() = delete;

  //! Assignment operator (disabled)
  SplitR2D& operator = (const SplitR2D&) = delete;

 private:
  const std::vector< MatPoly<2> >& matpolys_;
  const Plane_t<2>& cutting_plane_;
  double vol_tol_;
  double dst_tol_;
  bool all_convex_;
};

/*!
 * \class ClipR2D  2-D clipping algorithm
 *
 * For a vector of MatPoly's, computes the aggregated moments of their 
 * parts below the cutting line.
 * This can be interpreted as moments of chopped off parts in the nested 
 * dissections algorithm.
 * 
 * In this routine, we utilize the fact that R2D can clip a non-convex 
 * polygon with a plane.
 * 
 * For compatibility with the 3D functor we use the planar_faces flag, 
 * which does not have an effect in 2D.
 *
 */

class ClipR2D {
 public:
  ClipR2D(double vol_tol) : 
          vol_tol_(vol_tol) {}

  /*! 
    @brief Set the cutting line used by the functor. 
    @param[in] cutting_line Cutting line to clip with. 
  */
  void set_plane(const Plane_t<2>& cutting_line) {
    for (int ixy = 0; ixy < 2; ixy++)
      r2d_cut_line_.n.xy[ixy] = -cutting_line.normal[ixy];
    r2d_cut_line_.d = -cutting_line.dist2origin;
  }

  /*! 
    @brief Set the material polygons to be used by the functor.
    @param[in] matpolys Material polygons to be clipped.
    @param[in] planar_faces Flag kept for compatibility with the
    3D version.
  */
  void set_matpolys(const std::vector< MatPoly<2> >& matpolys,
                    const bool planar_faces = true) {
    int npolys = static_cast<int>(matpolys.size());
    r2d_polys_.resize(npolys);

    for(int ipoly = 0; ipoly < npolys; ipoly++)
      matpoly_to_r2dpoly(matpolys[ipoly], r2d_polys_[ipoly]);
  }

  /*! 
    @brief Computes moments of provided material polygons 
    @return Vector containing aggregated moments of all
    provided material polygons
  */
  std::vector<double> aggregated_moments() {
    std::vector<double> agg_moments;
    if (!r2d_polys_.empty()) {
      agg_moments.assign(3, 0.0);

      const int POLY_ORDER = 1;
      r2d_real r2d_moments[R2D_NUM_MOMENTS(POLY_ORDER)];
      for(int ipoly = 0; ipoly < r2d_polys_.size(); ipoly++)
        if (r2d_polys_[ipoly].nverts != 0) {
          r2d_reduce(&r2d_polys_[ipoly], r2d_moments, POLY_ORDER);
          for (int im = 0; im < 3; im++)
            agg_moments[im] += r2d_moments[im];
        }
    }

    if (agg_moments[0] < vol_tol_)
      agg_moments.assign(3, 0.0);

    return agg_moments;
  }

  /*! 
    @brief Computes moments of chopped off components: 
    the parts of MatPoly's below the line
    @return Vector containing aggregated moments 
    of the chopped off part of all MatPoly's
   */
  std::vector<double> operator() () const {
    std::vector<double> below_plane_moments(3, 0.0);

    const int POLY_ORDER = 1;
    r2d_real r2d_moments[R2D_NUM_MOMENTS(POLY_ORDER)];

    for (int ipoly = 0; ipoly < r2d_polys_.size(); ipoly++) {
      if (r2d_polys_[ipoly].nverts != 0) {
        //r2d does in-place clipping and does not have the const modifier
        //for the line: we need to make copies
        r2d_poly clipped_poly = r2d_polys_[ipoly];
        r2d_plane line_copy = r2d_cut_line_;
        r2d_clip(&clipped_poly, &line_copy, 1);
        if (clipped_poly.nverts != 0) {
          //Find the moments for the part in the lower halfspace
          r2d_reduce(&clipped_poly, r2d_moments, POLY_ORDER);
          for (int im = 0; im < 3; im++)
            below_plane_moments[im] += r2d_moments[im];
        }
      }
    }

    if (below_plane_moments[0] < vol_tol_)
      below_plane_moments.assign(3, 0.0);
      
    return below_plane_moments;
  }

  //! Assignment operator (disabled)
  ClipR2D& operator = (const ClipR2D&) = delete;

 private:
  std::vector<r2d_poly> r2d_polys_;
  r2d_plane r2d_cut_line_;
  double vol_tol_;
};

} // namespace Tangram

#endif
