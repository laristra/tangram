/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#ifndef VOF_H
#define VOF_H

#include <vector>
#include "tangram/support/tangram.h"
#include "tangram/support/MatPoly.h"
#include "tangram/support/lsfits.h"
#include "tangram/reconstruct/nested_dissections.h"
#include "tangram/reconstruct/cutting_distance_solver.h"

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
    @param[in] im_tols Tolerances for iterative methods
    @param[in] all_convex Flag indicating whether all mesh cells are convex
  */
  explicit VOF(const Mesh_Wrapper& Mesh, 
               const IterativeMethodTolerances_t& im_tols,
               const bool all_convex = false) : 
               mesh_(Mesh), im_tols_(im_tols), all_convex_(all_convex) {}
  
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
      cell_mat_ids_[icell].assign(cell_mat_ids.begin() + offset, 
                                  cell_mat_ids.begin() + offset + nmats);
      cell_mat_vfracs_[icell].assign(cell_mat_volfracs.begin() + offset, 
                                     cell_mat_volfracs.begin() + offset + nmats);
      offset += nmats;
    }
  }
  
  /*!
    @brief Used iterative methods tolerances
    @return  Tolerances for iterative methods, 
    here im_tols_.fun_eps is the volume tolerance
  */
  const IterativeMethodTolerances_t& iterative_method_tolerances() const {
    return im_tols_;
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
    @brief Calculate a position of a plane that clips off a particular material
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
    double vol_tol = im_tols_.fun_eps;

    std::vector<int> istencil_cells;
    mesh_.cell_get_node_adj_cells(cellID, Entity_type::ALL, &istencil_cells);
    istencil_cells.insert(istencil_cells.begin(), cellID);
    int nsc = (int) istencil_cells.size();

    // Create stencil for the volume fractions gradient: the first entry corresponds
    // to the current cell, the rest correspond to all its neighbors through the nodes
    std::vector<double> stencil_vfracs(nsc);
    std::vector< Point<Dim> > stencil_centroids(nsc);
    for (int isc = 0; isc < nsc; isc++) {
      const std::vector<int>& cur_mat_ids = cell_mat_ids_[istencil_cells[isc]];
      int local_id = (int) (std::find(cur_mat_ids.begin(), cur_mat_ids.end(), matID) -
                            cur_mat_ids.begin());
      if (local_id != cur_mat_ids.size())
        stencil_vfracs[isc] = cell_mat_vfracs_[istencil_cells[isc]][local_id];

      mesh_.cell_centroid(istencil_cells[isc], &stencil_centroids[isc]);
    }

    double target_vol = stencil_vfracs[0]*mesh_.cell_volume(cellID);

    // Use least squares to compute the gradient
    cutting_plane.normal = -ls_gradient(stencil_centroids, stencil_vfracs);
    cutting_plane.normal.normalize();

    //Create cutting distance solver: if not all cells are convex, we assume 
    //faces to be non-planar
    CuttingDistanceSolver<Dim, MatPoly_Clipper> 
      solve_cut_dst(mixed_polys, cutting_plane.normal, im_tols_, all_convex_);

    solve_cut_dst.set_target_volume(target_vol); 
    std::vector<double> clip_res = solve_cut_dst();
    cutting_plane.dist2origin = clip_res[0];

#ifdef DEBUG
    // Check if the resulting volume matches the reference value
    double cur_vol_err = std::fabs(clip_res[1] - target_vol);
    if (cur_vol_err > vol_tol) 
      std::cerr << "VOF for cell " << cellID << ": after " << im_tols_.max_num_iter <<
        " iteration(s) achieved error in volume for material " << 
        matID << " is " << cur_vol_err << ", volume tolerance is " << vol_tol << std::endl;
#endif
  }

  /*!
    @brief Given a cell index, calculate the CellMatPoly using the VOF 
    interface reconstruction method
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

    // Use the nested dissections algorithm for multi-material cells
    NestedDissections<VOF, Dim, MatPoly_Splitter> 
      nested_dissections(*this, cellID, all_convex_);

    // We clip material in the same order they are given for the cell
    int nmats = (int) cell_mat_ids_[cellID].size();
    std::vector<int> direct_order(nmats);
    std::iota(direct_order.begin(), direct_order.end(), 0);
    nested_dissections.set_cell_materials_order(direct_order);

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
  const IterativeMethodTolerances_t im_tols_;
  const bool all_convex_;
  std::vector< std::vector<int> > cell_mat_ids_;
  std::vector< std::vector<double> > cell_mat_vfracs_;
  std::vector<int> icells_to_reconstruct;
};  // class VOF

} // namespace Tangram

#endif  