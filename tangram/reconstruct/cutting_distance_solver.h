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
    @brief Sets the reference volume of the components below the plane
    @param target_vol Volume to match
  */
  void set_target_volume(double const target_vol) {
#ifdef DEBUG
    assert(target_vol >= 0.0);
#endif    
    target_vol_ = target_vol;
  }

  /*! 
    @brief Computes the cutting distance for the plane of the given orientation
    that matches the reference volume
    @return Vector containing the cutting distance followed by aggregated moments 
    of the MatPoly's components below the plane
   */
  std::vector<double> operator() () const {
    double vol_eps = tolerances_.fun_eps;

    double dst_bnd[2] = { -DBL_MAX, DBL_MAX };
    // Find the cutting distance corresponding to the planes passing through
    // the nearest and the farthest vertices
    for (int ipoly = 0; ipoly < matpolys_.size(); ipoly++) {
      const std::vector< Point<D> >& poly_pts = matpolys_[ipoly].points();
      for (int ivrt = 0; ivrt < matpolys_[ipoly].num_vertices(); ivrt++) {
        // dist2origin for Plane_t is defined as dot(PO, normal), where PO
        // is the vector from the point to the origin and is the negative of P.asV()
        double dst_vrt = -dot(plane_normal_, poly_pts[ivrt].asV());
        // dst_bnd is associated with the volume below the plane with dst_bnd[0] 
        // corresponding to the zero volume, so the sign is reversed
        if (dst_vrt > dst_bnd[0]) dst_bnd[0] = dst_vrt;
        if (dst_vrt < dst_bnd[1]) dst_bnd[1] = dst_vrt;
      }
    }

    if (target_vol_ < vol_eps) {
      // Size of the clipped polys is zero: we return the distance corresponding 
      // to the nearest vertex and zero moments
      std::vector<double> empty(D + 2, 0.0);
      empty[0] = dst_bnd[0];
      return empty;
    }

    // Aggregate moments of the cut MatPoly's
    std::vector<double> full_moments(D + 1, 0.0);
    for (int ipoly = 0; ipoly < matpolys_.size(); ipoly++) {
      const std::vector<double>& poly_moments = matpolys_[ipoly].moments();
      for (int im = 0; im < D + 1; im++)
        full_moments[im] += poly_moments[im];
    }

    if (target_vol_ > full_moments[0] - vol_eps) {
      // All the poly's are below the plane: we return the distance corresponding
      // to the farthest vertex and aggregate moments
      full_moments.insert(full_moments.begin(), dst_bnd[1]);
      return full_moments;
    }

    double d2orgn = 0.5*(dst_bnd[0] + dst_bnd[1]);
    // As distance changes from dst_bnd[0] to dst_bnd[1], the volume changes
    // from 0 to the aggregated volume of matpolys
    double vol_bnd[2] = { 0.0, full_moments[0] };
    double vol_err = full_moments[0];  // Current volume mismatch
    double dvol = full_moments[0]; // Change in volume over secant method's iteration
    double cur_dst_bnd[2] = { dst_bnd[0], dst_bnd[1] };
    std::vector<double> cur_moments;

    bool use_secant = true;
    int iter_count = 0;
    
    Plane_t<D> cutting_plane = {.normal = plane_normal_, .dist2origin = d2orgn};
    MatPoly_Clipper clip_matpolys(matpolys_, cutting_plane, planar_faces_);

    // We iterative solve for the cutting distance that corresponds to the target volume
    // below the cutting plane. By default, the secant method is used.
    // If the secant method fails, we discard its results, fall back to the initial guess 
    // and bounds, then use the bisection algorithm instead.
    while (vol_err > vol_eps) {
      double sec_coef;
      if (use_secant) {
        // Both the change in volume on the previous step and the max possible change
        // in volume on this step should be above volume tolerance to keep using 
        // the secant method
        use_secant = (std::fabs(dvol) > vol_eps) &&
                     (std::fabs(vol_bnd[1] - vol_bnd[0]) > vol_eps);
        if (use_secant) {
          sec_coef = (cur_dst_bnd[1] - cur_dst_bnd[0])/(vol_bnd[1] - vol_bnd[0]);
          d2orgn = cur_dst_bnd[1] + sec_coef*(target_vol_ - vol_bnd[1]);
          use_secant = ((d2orgn > dst_bnd[0]) == (d2orgn < dst_bnd[1]));          
        }

        if(!use_secant) {
          // Secant method failed: fall back to the bisection algorithm
          bool valid_lower_bnd = (vol_bnd[0] < target_vol_);
          if ( valid_lower_bnd != (vol_bnd[1] > target_vol_) ) {
            // Target volume fraction is no longer within the current bounds: 
            // reset current bounds
            cur_dst_bnd[0] = dst_bnd[0]; cur_dst_bnd[1] = dst_bnd[1];
          }
          else if (!valid_lower_bnd) {
            // vfrac_bnd[0] is the upper bound: swap the current bounds
            std::reverse(std::begin(cur_dst_bnd), std::end(cur_dst_bnd));
          }

          d2orgn = 0.5*(cur_dst_bnd[0] + cur_dst_bnd[1]);
        }
      }
      else
        d2orgn = 0.5*(cur_dst_bnd[0] + cur_dst_bnd[1]);

      if (std::fabs(cur_dst_bnd[1] - cur_dst_bnd[0]) <= tolerances_.arg_eps) {
        // The currently used method is valid and possible change in 
        // cutting distance is negligible, we can terminate
        cur_moments.insert(cur_moments.begin(), cutting_plane.dist2origin);
        return cur_moments;
      } 

      cutting_plane.dist2origin = d2orgn;
      cur_moments = clip_matpolys();  
      
      vol_err = std::fabs(cur_moments[0] - target_vol_);
      iter_count++;
      if (iter_count > tolerances_.max_num_iter) {
        // Max number of iterations reached: we return the current distance.
        // Because the moments are also returned, the caller can calculate the
        // discrepancy in volumes and take appropriate action
        break;
      }

      if (!use_secant) {
        if (cur_moments[0] > target_vol_)
          cur_dst_bnd[1] = d2orgn;
        else
          cur_dst_bnd[0] = d2orgn;
      }
      else {
        // Compute the change in volume over this step
        dvol = cur_moments[0] - vol_bnd[1];

        cur_dst_bnd[0] = cur_dst_bnd[1];
        cur_dst_bnd[1] = d2orgn;
        vol_bnd[0] = vol_bnd[1];
        vol_bnd[1] = cur_moments[0];
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
  double target_vol_ = 0.0;
  bool planar_faces_;
};

} // namespace Tangram

#endif // CUTTING_DISTANCE_SOLVER_H