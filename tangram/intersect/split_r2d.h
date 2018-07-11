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
split_convex_matpoly_r2d(const MatPoly<2>& mat_poly,
                         const Plane_t<2>& cutting_plane,
                         MatPoly<2>& lower_halfspace_poly,
                         MatPoly<2>& upper_halfspace_poly,
                         std::vector<double>& lower_halfspace_moments,
                         std::vector<double>& upper_halfspace_moments) {

  //Translate the cutting plane to R2D format
  r2d_plane r2d_cut_plane;
  for (int ixy = 0; ixy < 2; ixy++)
    r2d_cut_plane.n.xy[ixy] = cutting_plane.normal[ixy];
  r2d_cut_plane.d = cutting_plane.dist2origin;

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
  for (int isp = 0; isp < 2; isp++) {
     int nverts = r2d_subpolys[isp].nverts;

    //Check if the subpoly is empty
    if (nverts == 0) {
      subpoly_ptrs[isp]->clear();
      subpoly_moments_ptrs[isp]->clear();
      continue;
    }

    //Find the moments for a subpoly
    r2d_reduce(&r2d_subpolys[isp], r2d_moments, POLY_ORDER);

    if (r2d_moments[0] <= std::numeric_limits<double>::epsilon()) {
      subpoly_ptrs[isp]->clear();
      subpoly_moments_ptrs[isp]->clear();
      continue;
    }

    subpoly_moments_ptrs[isp]->assign(r2d_moments, r2d_moments + 3);

    //Get a MatPoly for a r2d subpoly
    const std::vector<Point2> matpoly_verts;
    matpoly_verts.resize(nverts);
    
    for (int v = 0; v < nverts; v++)
     for (int ixy = 0; ixy < 2; ixy++)
      matpoly_vrts[v][ixy] = r2d_subpolys[isp].verts[v].pos.xy[ixy];

    Matpoly<2> sub_matpoly;
    sub_matpoly.initialize(matpoly_verts);
    *subpoly_ptrs[isp] = sub_matpoly; 
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
           const bool all_convex = false) : 
           matpolys_(matpolys), 
           cutting_plane_(cutting_plane),
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
                               cur_moments[0], cur_moments[1]);

      for (int ihs = 0; ihs < 2; ihs++)
        if (!cur_moments[ihs].empty()) {
          hs_subpolys_ptrs[ihs]->emplace_back(cur_subpolys[ihs]);
          for (int im = 0; im < 3; im++)
            (*hs_moments_ptrs[ihs])[im] += cur_moments[ihs][im];
          hs_poly_count[ihs]++;
        }
    }

    for (int ihs = 0; ihs < 2; ihs++)
      if (hs_poly_count[ihs] == 0)
        hs_moments_ptrs[ihs]->clear();

    return hs_sets;
  }

  //! Default constructor (disabled)
  SplitR2D() = delete;

  //! Assignment operator (disabled)
  SplitR2D& operator = (const SplitR2D&) = delete;

 private:
  const std::vector< MatPoly<2> >& matpolys_;
  const Plane_t<2>& cutting_plane_;
  bool all_convex_;
};


/*!
  @brief Moments of MatPoly's components below the cutting plane.
  For a non-convex MatPoly no decomposition is performed.
  @param[in] mat_poly MatPoly to split
  @param[in] cutting_plane Cutting plane to split with
  @param[out] lower_halfspace_moments Computed moments
*/
void
lower_halfspace_moments_r2d(const MatPoly<2>& mat_poly,
                            const Plane_t<2>& cutting_plane,
                            std::vector<double>& lower_halfspace_moments) {
 
  //Translate the cutting plane to R2D format
  r2d_plane r2d_cut_plane;
  for (int ixy = 0; ixy < 2; ixy++)
    r2d_cut_plane.n.xy[ixy] = -cutting_plane.normal[ixy];
  r2d_cut_plane.d = -cutting_plane.dist2origin;

  //Convert matpoly to r2dpoly and clip
  r2d_poly r2dized_poly;
  matpoly_to_r2dpoly(mat_poly, r2dized_poly);
  r2d_clip(&r2dized_poly, &r2d_cut_plane, 1);

  //Compute moments for subpoly
  const int POLY_ORDER = 1;
  r2d_real r2d_moments[R2D_NUM_MOMENTS(POLY_ORDER)];
  //Check if the subpoly is empty
  if (r2dized_poly.nverts == 0)
    lower_halfspace_moments.clear();
  else {
    //Find the moments for a halfspace
    r2d_reduce(&r2dized_poly, r2d_moments, POLY_ORDER);
    lower_halfspace_moments.assign(r2d_moments, r2d_moments + 3);
  }
}//lower_halfspace_moments_r2d

/*!
 * \class ClipR2D  2-D clipping algorithm
 *
 * For a vector of MatPoly's, computes the aggregated moments of their 
 * parts below the cutting plane.
 * This can be interpreted as moments of chopped off parts in the nested 
 * dissections algorithm.
 * 
 * In this routine, we utilize the fact that R2D can clip a non-convex 
 * polyhedron with a plane.
 * 
 * By default, the faces of MatPoly's are assumed to be non-planar
 * and are facetized. If planar_faces is set to true,
 * faces of all MatPoly's are assumed to be planar and are NOT facetized.
 *
 */

class ClipR2D {
 public:
  ClipR2D(const std::vector< MatPoly<2> >& matpolys,
          const Plane_t<2>& cutting_plane, 
          const bool planar_faces = false) : 
          matpolys_(matpolys), 
          cutting_plane_(cutting_plane),
          planar_faces_(planar_faces) {}

  /*! 
    @brief Computes moments of chopped off components: 
    the parts of MatPoly's below the plane
    @return Vector containing aggreagted moments 
    of the chopped off part of all MatPoly's
   */
  std::vector<double> operator() () const {
    std::vector<double> below_plane_moments(3, 0.0);

    int poly_counter = 0;
    if (planar_faces_)
      for (int ipoly = 0; ipoly < matpolys_.size(); ipoly++) {
        std::vector<double> cur_moments;
        lower_halfspace_moments_r2d(matpolys_[ipoly], cutting_plane_, 
                                    cur_moments);

        if (!cur_moments.empty()) {
          for (int im = 0; im < 3; im++)
            below_plane_moments[im] += cur_moments[im];
          poly_counter++;
        }
      }
    else 
      for (int ipoly = 0; ipoly < matpolys_.size(); ipoly++) {
        MatPoly<2> faceted_poly;
        matpolys_[ipoly].faceted_matpoly(&faceted_poly);
        std::vector<double> cur_moments;
        lower_halfspace_moments_r2d(faceted_poly, cutting_plane_, 
                                    cur_moments);

        if (!cur_moments.empty()) {
          for (int im = 0; im < 3; im++)
            below_plane_moments[im] += cur_moments[im];
          poly_counter++;
        }                            
      }
    
    if (poly_counter == 0)
      below_plane_moments.clear();
      
    return below_plane_moments;
  }

//! Default constructor (disabled)
  ClipR2D() = delete;

  //! Assignment operator (disabled)
  ClipR2D& operator = (const ClipR2D&) = delete;

 private:
  const std::vector< MatPoly<2> >& matpolys_;
  const Plane_t<2>& cutting_plane_;
  bool planar_faces_;
};

} // namespace Tangram

#endif
