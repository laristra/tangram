/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#ifndef CUTTING_DISTANCE_SOLVER_H
#define CUTTING_DISTANCE_SOLVER_H

#include <cfloat>
#include "tangram/support/tangram.h"
#include "tangram/support/MatPoly.h"


/*!
 @file cutting_distance_solver.h
 @brief Given an orientation of the cutting plane, uses secant-bisection 
 algorithm to solve for the cutting distance that matches the target 
 volume fraction to the given tolerance. 
 All the provided material polygons are cut with the same plane.
 The unique position of the plane is determined based on the aggregate 
 moments of the components in the lower half-space (below the plane).
 */

namespace Tangram {

/*!
 * \class CuttingDistanceSolver Secant-bisection algorithm
 *
 * For a vector of MatPoly's and orientation of the cutting plane, 
 * finds the position of the plane for which the aggregated moments of their 
 * parts in the lower half-space match the specified volume fraction to the
 * given tolerance.
 * 
 * By default, the faces of MatPoly's are assumed to be planar
 * and are NOT facetized. If planar_faces is set to false,
 * faces of all MatPoly's are assumed to be non-planar and are facetized.
 *
 */

template <int D, class MatPoly_Clipper>
class CuttingDistanceSolver {
public:  
  CuttingDistanceSolver(const std::vector< MatPoly<D> >& matpolys,
                        const Vector<D>& plane_normal,
                        const IterativeMethodTolerances_t& tolerances,
                        const bool planar_faces = true) : 
                        matpolys_(matpolys), plane_normal_(plane_normal),
                        tolerances_(tolerances), planar_faces_(planar_faces) {}
  /*! 
    @brief Sets the reference volume fraction of the components below the plane
    @param target_vfraction Value of the volume fraction to match
  */
  void set_volume_fraction(double const target_vfraction) {
#ifdef DEBUG
    assert(target_vfraction >= 0.0);
#endif    
    target_vfrac_ = target_vfraction;
  }

  /*! 
    @brief Computes the cutting distance for the plane of the given orientation
    that matches the reference volume fraction
    @return Vector containing the cutting distance followed by aggregated moments 
    of the MatPoly's components below the plane
   */
  std::vector<double> operator() () const {
    std::vector<double> full_moments(D + 1, 0.0);
    double dst_bnd[2] = { DBL_MAX, 0.0 };
    Point<D> nearest_pt;
    for (int ipoly = 0; ipoly < matpolys_.size(); ipoly++) {
      const std::vector<double>& poly_moments = matpolys_[ipoly].moments();
      for (int im = 0; im < D + 1; im++)
        full_moments[im] += poly_moments[im];

      const std::vector< Point<D> >& poly_pts = matpolys_[ipoly].points();
      for (int ivrt = 0; ivrt < matpolys_[ipoly].num_vertices(); ivrt++) {
        double dst_vrt = -dot(plane_normal_, poly_pts[ivrt].asV());
        if (dst_vrt < dst_bnd[0]) {
          dst_bnd[0] = dst_vrt; nearest_pt = poly_pts[ivrt];
        }
        if (dst_vrt > dst_bnd[1]) dst_bnd[1] = dst_vrt;
      }
    }
    if (dot(plane_normal_, nearest_pt.asV()) + dst_bnd[1] > 0.0)
      std::reverse(std::begin(dst_bnd), std::end(dst_bnd));

    if (target_vfrac_ < tolerances_.fun_eps) {
      // Size of the clipped polys is zero: we return the distance corresponding 
      // to the nearest vertex and zero moments
      std::vector<double> empty(D + 2, 0.0);
      empty[0] = dst_bnd[0];
      return empty;
    }
    else if (target_vfrac_ > 1.0 - tolerances_.fun_eps) {
      // All the poly's are below the plane: we return the distance corresponding
      // to the farthest vertex and aggregate moments
      full_moments.insert(full_moments.begin(), dst_bnd[1]);
      return full_moments;
    }
 
    double d2orgn = 0.5*(dst_bnd[0] + dst_bnd[1]);
    double vfrac_bnd[2] = { 0.0, 1.0 };
    double dvfrac = 1.0;
    double cur_dst_bnd[2] = { dst_bnd[0], dst_bnd[1] };
    double cur_vfrac;
    std::vector<double> cur_moments;

    bool use_secant = true;
    int iter_count = 0;
    
    Plane_t<D> cutting_plane = {.normal = plane_normal_, .dist2origin = d2orgn};
    MatPoly_Clipper clip_matpolys(matpolys_, cutting_plane, planar_faces_);

    while (dvfrac > tolerances_.fun_eps) {
      double sec_coef;
      if (use_secant) {
        use_secant = (std::fabs(vfrac_bnd[1] - vfrac_bnd[0]) > tolerances_.fun_eps);
        if (use_secant) {
          sec_coef = (cur_dst_bnd[1] - cur_dst_bnd[0])/(vfrac_bnd[1] - vfrac_bnd[0]);
          use_secant = (std::fabs(sec_coef) > tolerances_.fun_eps);
          if (use_secant) {
            d2orgn = cur_dst_bnd[1] + sec_coef*(target_vfrac_ - vfrac_bnd[1]);
            use_secant = ((d2orgn > dst_bnd[0]) == (d2orgn < dst_bnd[1]));          
          }
        }

        if(!use_secant) {
          // Secant method failed: reset variables and fall back to the bisection algorithm
          cur_dst_bnd[0] = dst_bnd[0]; cur_dst_bnd[1] = dst_bnd[1];
          d2orgn = 0.5*(dst_bnd[0] + dst_bnd[1]);
        }
      }
      else
        d2orgn = 0.5*(cur_dst_bnd[0] + cur_dst_bnd[1]);
      
      cutting_plane.dist2origin = d2orgn;
      cur_moments = clip_matpolys();  
      cur_vfrac = cur_moments[0]/full_moments[0];
      
      if (!use_secant) {
        if (cur_vfrac > target_vfrac_)
          cur_dst_bnd[1] = d2orgn;
        else
          cur_dst_bnd[0] = d2orgn;
      }
      else {
        cur_dst_bnd[0] = cur_dst_bnd[1];
        cur_dst_bnd[1] = d2orgn;
        vfrac_bnd[0] = vfrac_bnd[1];
        vfrac_bnd[1] = cur_vfrac;
      }

      dvfrac = std::fabs(cur_vfrac - target_vfrac_);
      iter_count++;
      if (iter_count > tolerances_.max_num_iter) {
        // Max number of iterations reached: we return the current distance.
        // Because the moments are also returned, the caller can calculate the
        // discrepancy in volumes and take appropriate action
        cur_moments.insert(cur_moments.begin(), d2orgn);
        return cur_moments;
      }
    }
    
    cur_moments.insert(cur_moments.begin(), d2orgn);
    return cur_moments;
  }


  //! Default constructor (disabled)
  CuttingDistanceSolver() = delete;

  //! Assignment operator (disabled)
  CuttingDistanceSolver& operator = (const CuttingDistanceSolver&) = delete;

 private:
  const std::vector< MatPoly<D> >& matpolys_;
  const Vector<D>& plane_normal_;
  const IterativeMethodTolerances_t& tolerances_;
  double target_vfrac_ = 0.0;
  bool planar_faces_;
};

} // namespace Tangram

#endif // CUTTING_DISTANCE_SOLVER_H