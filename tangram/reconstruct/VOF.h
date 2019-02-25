/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#ifndef TANGRAM_RECONSTRUCT_VOF_H_
#define TANGRAM_RECONSTRUCT_VOF_H_

#include <vector>

// tangram includes
#include "tangram/support/tangram.h"
#include "tangram/support/MatPoly.h"
#include "tangram/reconstruct/nested_dissections.h"
#include "tangram/reconstruct/cutting_distance_solver.h"

// wonton includes
#include "wonton/support/lsfits.h"


/*!
 @file VOF.h
  @brief Calculates the interface and constructs CellMatPoly using the VOF
  algorithm.

  @tparam Mesh_Wrapper A lightweight wrapper to a specific input mesh
  implementation that provides certain functionality
  @tparam Dim The spatial dimension of the problem
  @tparam MatPoly_Splitter An operator for splitting a vector of MatPoly's
  into half-space sets with a cutting plane
  @tparam MatPoly_Clipper An operator for computing moments of
  components of a vector of MatPoly's below a cutting plane
 */

namespace Tangram {

template <class Mesh_Wrapper, int Dim, class MatPoly_Splitter, class MatPoly_Clipper>
class VOF {
public:
  /*!
    @brief Constructor for a VOF interface reconstruction algorithm
    @param[in] Mesh A lightweight wrapper to a specific input mesh
    implementation that provides certain functionality
    @param[in] ims_tols Tolerances for iterative methods
    @param[in] all_convex Flag indicating whether all mesh cells are convex
  */
  explicit VOF(const Mesh_Wrapper& Mesh,
               const std::vector<IterativeMethodTolerances_t>& ims_tols,
               const bool all_convex = false) :
               mesh_(Mesh), ims_tols_(ims_tols), all_convex_(all_convex) {
    if (ims_tols.empty())
      throw std::runtime_error(
        "VOF uses 0-order moments and needs tolerances for the related iterative method!");
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
    cell_mat_ids_.clear();
    cell_mat_vfracs_.clear();
    int ncells = mesh_.num_owned_cells() + mesh_.num_ghost_cells();
    cell_mat_ids_.resize(ncells);
    cell_mat_vfracs_.resize(ncells);

    int offset = 0;
    for (int icell = 0; icell < ncells; icell++) {
      int nmats = cell_num_mats[icell];

      for (int icmat = 0; icmat < nmats; icmat++) {
        cell_mat_ids_[icell].push_back(cell_mat_ids[offset + icmat]);
        cell_mat_vfracs_[icell].push_back(cell_mat_volfracs[offset + icmat]);
      }
      offset += nmats;
    }
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

  /*!
    @brief Pass in indices of cells for which CellMatPoly objects
    are to be constructed. If the index is in the list, a CellMatPoly object will be
    created even for a single-material cell.
    @param[in] cellIDs_to_op_on A vector of length up to (num_cells)
    specifying the indices of cells for which CellMatPoly objects are requested.
  */
  void set_cell_indices_to_operate_on(std::vector<int> const& cellIDs_to_op_on) {
    icells_to_reconstruct = cellIDs_to_op_on;
  }

  /*!
    @brief Calculate the position of a plane that clips off a particular material.
    This method is used on every step of the nested dissections algorithm.
    Note that if MatPoly_Clipper can handle non-convex cells, this method
    does not require decomposion into tetrahedrons.
    @param[in] cellID Index of the multi-material cell to operate on
    @param[in] matID Index of the material to clip
    @param[in] mixed_polys Vector of material poly's that contain the material
    to clip and (possibly) other materials
    @param[out] cutting_plane The resulting cutting plane position
    @param[in] planar_faces Flag indicating whether the faces of all mixed_polys
    are planar
  */
  void get_plane_position(const int cellID,
                          const int matID,
                          const std::vector< MatPoly<Dim> >& mixed_polys,
                          Plane_t<Dim>& cutting_plane,
                          const bool planar_faces = true) const {
    double vol_tol = ims_tols_[0].fun_eps;

    std::vector<int> istencil_cells;
    mesh_.cell_get_node_adj_cells(cellID, Entity_type::ALL, &istencil_cells);
    istencil_cells.insert(istencil_cells.begin(), cellID);
    int nsc = static_cast<int>(istencil_cells.size());

    // Create stencil for the volume fractions gradient: the first entry corresponds
    // to the current cell, the rest correspond to all its neighbors through the nodes
    std::vector<double> stencil_vfracs(nsc);
    std::vector< Point<Dim> > stencil_centroids(nsc);
    for (int isc = 0; isc < nsc; isc++) {
      const std::vector<int>& cur_mat_ids = cell_mat_ids_[istencil_cells[isc]];
      int local_id = std::distance(cur_mat_ids.begin(),
        std::find(cur_mat_ids.begin(), cur_mat_ids.end(), matID));

      if (local_id != cur_mat_ids.size())
        stencil_vfracs[isc] = cell_mat_vfracs_[istencil_cells[isc]][local_id];

      mesh_.cell_centroid(istencil_cells[isc], &stencil_centroids[isc]);
    }

    double target_vol = stencil_vfracs[0]*mesh_.cell_volume(cellID);
    // Use least squares to compute the gradient
    cutting_plane.normal = -ls_gradient(stencil_centroids, stencil_vfracs);

    double grad_norm = cutting_plane.normal.norm();
    if (is_equal(grad_norm, 0.0)) {
      // Zero gradient: we choose SLIC-like plane orientation
      // Normal is set in the direction of the x-axis
      cutting_plane.normal.axis(0);
    }
    else
      cutting_plane.normal /= grad_norm;

    //Create cutting distance solver
    CuttingDistanceSolver<Dim, MatPoly_Clipper>
      solve_cut_dst(mixed_polys, cutting_plane.normal, ims_tols_[0], planar_faces);

    solve_cut_dst.set_target_volume(target_vol);
    std::vector<double> clip_res = solve_cut_dst();
    cutting_plane.dist2origin = clip_res[0];

#ifdef DEBUG
    // Check if the resulting volume matches the reference value
    double cur_vol_err = std::fabs(clip_res[1] - target_vol);
    if (cur_vol_err > vol_tol)
      std::cerr << "VOF for cell " << cellID << ": after " << ims_tols_[0].max_num_iter <<
        " iteration(s) achieved error in volume for material " <<
        matID << " is " << cur_vol_err << ", volume tolerance is " << vol_tol << std::endl;
#endif
  }

  /*!
    @brief Given a cell index, calculate the CellMatPoly using the VOF
    interface reconstruction method.
    Uses nested dissections algorithm.
  */
  std::shared_ptr<CellMatPoly<Dim>> operator()(const int cell_op_ID) const {
    int cellID = icells_to_reconstruct[cell_op_ID];

    // Check if the cell is single-material
    if (cell_mat_ids_[cellID].size() == 1) {
      std::shared_ptr< CellMatPoly<Dim> > cmp_ptr(new CellMatPoly<Dim>(cellID));
      MatPoly<Dim> cell_matpoly;
      cell_get_matpoly(mesh_, cellID, &cell_matpoly);
      cell_matpoly.set_mat_id(cell_mat_ids_[cellID][0]);
      cmp_ptr->add_matpoly(cell_matpoly);

      return cmp_ptr;
    }

    // Use the nested dissections algorithm for multi-material cells.
    // Note that nested dissections uses this instance of reconstructor
    // to invoke get_plane_position position method. Nested dissections
    // itself does not have its own MeshWrapper, MatPoly_Clipper, etc.,
    // they all are reconstructor specific
    NestedDissections<VOF, Dim, MatPoly_Splitter>
      nested_dissections(*this, cellID, all_convex_);

    // We clip material in the same order they are given for the cell
    // Note that this is the order of local materials, not material
    // indices. Nested dissections uses cell_materials method to get
    // actual material indices.
    nested_dissections.set_cell_materials_order(false);

    return nested_dissections();
  }

  /*!
    @brief Materials in the cell
    @param[in] cellID Cell index
    @return  Vector of indices of cell's materials
  */
  const std::vector<int>& cell_materials(const int cellID) const { return cell_mat_ids_[cellID]; }

  /*!
    @brief MatPoly corresponding to a mesh cell
    @param[in] cellID Cell index
    @return  MatPoly for this cell
  */
  MatPoly<Dim> cell_matpoly(const int cellID) const {
    MatPoly<Dim> mat_poly;
    cell_get_matpoly(mesh_, cellID, &mat_poly);

    return mat_poly;
  }

private:
  const Mesh_Wrapper& mesh_;
  const std::vector<IterativeMethodTolerances_t> ims_tols_;
  const bool all_convex_;
  std::vector< std::vector<int> > cell_mat_ids_;
  std::vector< std::vector<double> > cell_mat_vfracs_;
  std::vector<int> icells_to_reconstruct;
};  // class VOF

}  // namespace Tangram

#endif  // TANGRAM_RECONSTRUCT_VOF_H_
