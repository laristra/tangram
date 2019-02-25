/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#ifndef TANGRAM_RECONSTRUCT_SLIC_H_
#define TANGRAM_RECONSTRUCT_SLIC_H_

#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <memory>
#include <float.h>

#include "tangram/support/tangram.h"
#include "tangram/driver/CellMatPoly.h"
#include "tangram/support/MatPoly.h"
#include "tangram/reconstruct/cutting_distance_solver.h"

/*!
 @file SLIC.h
 @brief Simple implemenation of the crude Piecewise Linear Interface
 Reconstruction algorithm.

 Here, we use only vertical, "y"-aligned interfaces
 */

namespace Tangram {

  /*!
   @class SLIC "SLIC.h"
   @brief Calculates the interface and constructs CellMatPoly for the SLIC
   algorithm.

   @tparam Mesh_Wrapper A lightweight wrapper to a specific input mesh
   implementation that provides certain functionality
   @tparam Dim The spatial dimension of the problem.
   */

  template <class Mesh_Wrapper, int Dim, class MatPoly_Splitter, class MatPoly_Clipper>
  class SLIC {
  public:
    /*!
     @brief Constructor performing a SLIC algorithm for interface reconstruction.
     */
    explicit SLIC(const Mesh_Wrapper & Mesh,
                  const std::vector<IterativeMethodTolerances_t>& ims_tols,
                  const bool all_convex = false) :
                  mesh_(Mesh), ims_tols_(ims_tols), all_convex_(all_convex) {
      // For now
      if (ims_tols.empty())
        throw std::runtime_error(
          "SLIC uses 0-order moments and needs tolerances for the related iterative method!");
    }

    /*!
     @brief Pass in the volume fraction data for use in the reconstruction.
     @param[in] cell_num_mats A vector of length (num_cells) specifying the
     number of materials in each cell.
     @param[in] cell_mat_ids A vector of length (sum(cell_num_mats)) specifying
     the ID of each material in each cell
     @param[in] cell_mat_volfracs A vector of length(sum(cell_num_mats))
     specifying the volume fraction of each material in each cell.
     */
    void set_volume_fractions(std::vector<int> const& cell_num_mats,
                              std::vector<int> const& cell_mat_ids,
                              std::vector<double> const& cell_mat_volfracs,
                              std::vector<Point<Dim>>
                              const& cell_mat_centroids = {}) {
      cell_num_mats_ = cell_num_mats;
      cell_mat_ids_ = cell_mat_ids;
      cell_mat_volfracs_ = cell_mat_volfracs;
      int nc = mesh_.num_owned_cells() + mesh_.num_ghost_cells();
      cell_mat_offsets_.resize(nc);
      cell_mat_offsets_[0] = 0;
      for (int c(1); c < nc; ++c)
        cell_mat_offsets_[c] = cell_mat_offsets_[c-1] + cell_num_mats_[c-1];
    }

    /*!
      @brief Used iterative methods tolerances
      @return  Tolerances for iterative methods,
      here ims_tols_[0] correspond to methods for volumes
      and ims_tols_[1] are NOT used.
      In particular, ims_tols_[0].arg_eps is a negligible
      change in cutting distance, ims_tols_[0].fun_eps is a
      negligible discrepancy in volume.
    */
    const std::vector<IterativeMethodTolerances_t>&
    iterative_methods_tolerances() const {
      return ims_tols_;
    }

    void set_cell_indices_to_operate_on(std::vector<int> const& cellIDs_to_op_on) {
      icells_to_reconstruct = cellIDs_to_op_on;
    }

    /*!
     @brief Given a cell index, calculate the CellMatPoly for this reconstruction
     */
    std::shared_ptr<CellMatPoly<Dim>> operator()(const int cell_op_ID) const {
      double vol_tol = ims_tols_[0].fun_eps;  // Volume tolerance

      int cellID = icells_to_reconstruct[cell_op_ID];
      auto numMats = cell_num_mats_[cellID];

      CellMatPoly<Dim>* cellpoly = new CellMatPoly<Dim>(cellID);

      auto iStart = cell_mat_offsets_[cellID];

      //Sets of MatPoly's on two sides of the cutting plane
      HalfSpaceSets_t<Dim> hs_sets;

      // For every chopped off material, we split MatPoly's above the plane
      hs_sets.upper_halfspace_set.matpolys.resize(1);
      Tangram::cell_get_matpoly(mesh_, cellID, &hs_sets.upper_halfspace_set.matpolys[0]);

      // Just going along x-direction
      Plane_t<Dim> cutting_plane;
      for (int i = 0 ; i < Dim; i++)
       cutting_plane.normal[i] = 0.0;
      cutting_plane.normal[0] = 1.0;  

      //Create Splitter instance
      MatPoly_Splitter split_matpolys(hs_sets.upper_halfspace_set.matpolys,
                                      cutting_plane, all_convex_);

      //Create cutting distance solver: if not all cells are convex, we assume that
      //faces are non-planar
      CuttingDistanceSolver<Dim, MatPoly_Clipper>
        solve_cut_dst(hs_sets.upper_halfspace_set.matpolys,
                      cutting_plane.normal, ims_tols_[0], all_convex_);

      //Cutting from left to right
      double cell_volume = mesh_.cell_volume(cellID);
      for (int iMat(0); iMat < numMats; iMat++) {
        double target_vol = cell_mat_volfracs_[iStart + iMat]*cell_volume;
        // If the target volume is too small, skip it
        if (target_vol < vol_tol) continue;

        MatPolySet_t<Dim>* single_mat_set_ptr;
        
	//On the last iteration the remaining part is single-material,
        //so we don't need to split it
        if (iMat == numMats - 1)
          single_mat_set_ptr = &hs_sets.upper_halfspace_set;
        else {
          // Find distance to origin for the cutting plane
          solve_cut_dst.set_target_volume(target_vol);
          std::vector<double> clip_res = solve_cut_dst();
          cutting_plane.dist2origin = clip_res[0];

#ifdef DEBUG
          // Check if the resulting volume matches the reference value
          double cur_vol_err = std::fabs(clip_res[1] - target_vol);
          if (cur_vol_err > vol_tol)
            std::cerr << "SLIC for cell " << cellID << ": after " << ims_tols_[0].max_num_iter <<
              " iteration(s) achieved error in volume for material " <<
              cell_mat_ids_[iStart + iMat] << " is " << cur_vol_err <<
              ", volume tolerance is " << vol_tol << std::endl;
#endif
          hs_sets = split_matpolys();

          //Chopped off single-material MatPoly's are below the plane
          single_mat_set_ptr = &hs_sets.lower_halfspace_set;
        }
        //Add single-material MatPoly's to CellMatPoly
        for (int ismp = 0; ismp < single_mat_set_ptr->matpolys.size(); ismp++) {
          MatPoly<Dim>& cur_matpoly = single_mat_set_ptr->matpolys[ismp];
          cur_matpoly.set_mat_id(cell_mat_ids_[iStart+iMat]);
          cellpoly->add_matpoly(cur_matpoly);
    	}
      }

      return std::shared_ptr<CellMatPoly<Dim>>(cellpoly);
    }
  private:
    const Mesh_Wrapper & mesh_;
    const std::vector<IterativeMethodTolerances_t> ims_tols_;
    const bool all_convex_;
    std::vector<int> cell_num_mats_;
    std::vector<int> cell_mat_ids_;
    std::vector<double> cell_mat_volfracs_;
    std::vector<int> cell_mat_offsets_;
    std::vector<int> icells_to_reconstruct;
  };  // class SLIC
}  // namespace Tangram


#endif  // TANGRAM_RECONSTRUCT_SLIC_H_
