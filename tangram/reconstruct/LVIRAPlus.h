/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#ifndef TANGRAM_RECONSTRUCT_LVIRA_PLUS_H
#define TANGRAM_RECONSTRUCT_LVIRA_PLUS_H

#include <vector>
#include <functional>

// tangram includes
#include "tangram/support/tangram.h"
#include "tangram/support/MatPoly.h"
#include "tangram/reconstruct/nested_dissections.h"
#include "tangram/reconstruct/cutting_distance_solver.h"
#include "tangram/support/bfgs.h"

// wonton includes
#include "wonton/support/lsfits.h"

/*!
 @file LVIRAPlus.h
  @brief Calculates the interfaces and constructs CellMatPoly's using the
  Least squares Volume-of-fluid Interface Reconstruction Algorithm (LVIRA)
  modified for use on domains with more than two materials.
  
  @tparam Mesh_Wrapper A lightweight wrapper to a specific input mesh
  implementation that provides certain functionality
  @tparam Dim The spatial dimension of the problem
  @tparam MatPoly_Splitter An operator for splitting a vector of MatPoly's 
  into half-space sets with a cutting plane
  @tparam MatPoly_Clipper An operator for computing moments of 
  components of a vector of MatPoly's below a cutting plane
*/

namespace Tangram {

constexpr BFGS_ALG lvira_plus_bfgs_alg = BFGS;

template <class Mesh_Wrapper, int Dim, class MatPoly_Splitter, class MatPoly_Clipper>
class LVIRAPlus {
public:
  /*!
    @brief Constructor for the LVIRA+ interface reconstruction
    @param[in] Mesh A lightweight wrapper to a specific input mesh
    implementation that provides certain functionality
    @param[in] ims_tols Tolerances for iterative methods
    @param[in] all_convex Flag indicating whether all mesh cells are convex
  */
  explicit LVIRAPlus(const Mesh_Wrapper& Mesh, 
                     const std::vector<IterativeMethodTolerances_t>& ims_tols,
                     const bool all_convex) : 
                     mesh_(Mesh), ims_tols_(ims_tols), all_convex_(all_convex) {
    if (ims_tols.size() < 2) {
      std::string err_msg = "LVIRA+ uses optimization procedure to find the orientation of the cutting plane ";
      err_msg += "and then numerically solves for its position that matches the given volume fraction.";
      err_msg += "Tolerances should be provided for both iterative methods!\n";

      throw std::runtime_error(err_msg);
    }
    lump_aliens_with_first_mat_ = true;
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
    cell_neighbors_ids_.clear();

    int ncells = mesh_.num_owned_cells() + mesh_.num_ghost_cells();
    cell_mat_ids_.resize(ncells);
    cell_mat_vfracs_.resize(ncells);
    cell_neighbors_ids_.resize(ncells);
    nghb_alien_mat_vfracs_.resize(ncells);

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
    and ims_tols_[1] correspond to the optimization methods.
    In particular, ims_tols_[0].arg_eps is a negligible 
    change in cutting distance, ims_tols_[0].fun_eps is a 
    negligible discrepancy in volume, ims_tols_[1].arg_eps
    is a negligible change in the cutting plane orientation,
    and ims_tols_[1].fun_eps is a negligible error in given
    material volumes over neighboring cells, computed according
    to the modified LVIRA. The change in cutting plane
    orientation is defined as the norm of change of the cutting
    plane's normal, which is expressed in polar coordinates (angles).
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
    int nIR_cells = static_cast<int>(icells_to_reconstruct.size());

    for (int iirc = 0; iirc < nIR_cells; iirc++) {
      int icell = icells_to_reconstruct[iirc];
      // Find neighbors of the cell
      mesh_.cell_get_node_adj_cells(icell, Entity_type::ALL, &cell_neighbors_ids_[icell]);

      // Find aggregated volume fractions of materials that are contained
      // in the neighboring cells, but not in the cell itself
      int nnc = static_cast<int>(cell_neighbors_ids_[icell].size());
      nghb_alien_mat_vfracs_[icell].assign(nnc, 0.0);
      for (int inc = 0; inc < nnc; inc++) {
        int inghb_cell = cell_neighbors_ids_[icell][inc];
        int nnghb_mats = static_cast<int>(cell_mat_ids_[inghb_cell].size());

        for (int incm = 0; incm < nnghb_mats; incm++) {
          if (cell_mat_ids_[inghb_cell].empty()) {
            std::string err_msg = "LVIRA+ requires material data in the cells "; 
            err_msg += "neighboring the multi-material cells!\n";
            err_msg += "No material data is provided for cell " + std::to_string(inghb_cell);
            err_msg += " when reconstruction is requested for cell " + std::to_string(icell) + "!\n";
            throw std::runtime_error(err_msg);
          }

          if (std::find(cell_mat_ids_[icell].begin(), 
                        cell_mat_ids_[icell].end(), cell_mat_ids_[inghb_cell][incm]) ==
              cell_mat_ids_[icell].end())
          nghb_alien_mat_vfracs_[icell][inc] += cell_mat_vfracs_[inghb_cell][incm];
        }
      }
    }
  }

  /*!
    @brief Calculate the position of a plane that clips off a particular material.
    This method is used on every step of the nested dissections algorithm.
    Note that if MatPoly_Clipper can handle non-convex cells, this method
    does not require decomposion into tetrahedrons.
    @param[in] cellID Index of the multi-material cell to operate on
    @param[in] matID Index of the material to clip
    @param[in] mixed_polys Vector of pointers to vectors of material poly's 
    that contain the material to clip and (possibly) other materials
    @param[out] cutting_plane The resulting cutting plane position
    @param[in] planar_faces Flag indicating whether the faces of all mixed_polys
    are planar
  */  
  void get_plane_position(const int cellID,
                          const int matID,
                          const std::vector< std::vector< MatPoly<Dim> >* >& mixed_polys_ptrs,
                          Plane_t<Dim>& cutting_plane,
                          const bool planar_faces) const {
    assert(Dim > 1);
    
    double vol_tol = ims_tols_[0].fun_eps;
    int cellMatID = std::distance(cell_mat_ids_[cellID].begin(),
      std::find(cell_mat_ids_[cellID].begin(), 
                cell_mat_ids_[cellID].end(), matID));

    //For the initial guess we use the standard VOF method
    std::vector<int> iVOF_stencil_cells = neighbor_cells_to_split(cellID);
    iVOF_stencil_cells.insert(iVOF_stencil_cells.begin(), cellID);
    int nsc = static_cast<int>(iVOF_stencil_cells.size());

    // Create stencil for the volume fractions gradient: the first entry corresponds
    // to the current cell, the rest correspond to all its neighbors through the nodes
    std::vector<double> VOF_stencil_vfracs(nsc);
    std::vector< Point<Dim> > VOF_stencil_centroids(nsc);
    for (int isc = 0; isc < nsc; isc++) {
      const std::vector<int>& cur_mat_ids = cell_mat_ids_[iVOF_stencil_cells[isc]];
      int local_id =
        std::distance(cur_mat_ids.begin(),
          std::find(cur_mat_ids.begin(), cur_mat_ids.end(), matID));

      if (static_cast<int>(cur_mat_ids.size()) != local_id)
        VOF_stencil_vfracs[isc] = cell_mat_vfracs_[iVOF_stencil_cells[isc]][local_id];

      mesh_.cell_centroid(iVOF_stencil_cells[isc], &VOF_stencil_centroids[isc]);
    }

    // Use least squares to compute the gradient
    cutting_plane.normal = -ls_gradient(VOF_stencil_centroids, VOF_stencil_vfracs);

    double grad_norm = cutting_plane.normal.norm();
    if (is_equal(grad_norm, 0.0)) {
      // Zero gradient: we choose SLIC-like plane orientation
      // Normal is set in the direction of the x-axis
      cutting_plane.normal.axis(0);
    }
    else
      cutting_plane.normal /= grad_norm;

    Vector<Dim - 1> init_guess = cartesian_to_polar(cutting_plane.normal);

    std::function<double(const Vector<Dim - 1>&)> bfgs_obj_fun = 
      [this, &cellID, &cellMatID, &mixed_polys_ptrs]
      (const Vector<Dim - 1>& cur_arg)->double {
      return neighbors_error(cellID, cellMatID, mixed_polys_ptrs, 
                             polar_to_cartesian(cur_arg));
    };    

    double bfgs_obj_fun_lbnd = 0.0;

    Vector<Dim - 1> ang_min;
    switch(lvira_plus_bfgs_alg) {
      case BFGS: ang_min = bfgs<Dim - 1>(bfgs_obj_fun, bfgs_obj_fun_lbnd, 
                                         init_guess, ims_tols_[1]);
        break;
      case DBFGS: ang_min = dbfgs<Dim - 1>(bfgs_obj_fun, bfgs_obj_fun_lbnd, 
                                           init_guess, ims_tols_[1]);
        break;
      default: throw std::runtime_error("Unknown BFGS algorithm is selected for LVIRA!");
    }

    cutting_plane.normal = polar_to_cartesian(ang_min);

    CuttingDistanceSolver<Dim, MatPoly_Clipper> 
      solve_cut_dst(*(mixed_polys_ptrs[0]), cutting_plane.normal, ims_tols_[0], planar_faces);

    double target_vol = cell_mat_vfracs_[cellID][cellMatID]*mesh_.cell_volume(cellID);
    solve_cut_dst.set_target_volume(target_vol); 
    std::vector<double> clip_res = solve_cut_dst();
    cutting_plane.dist2origin = clip_res[0];

    // Check if the resulting volume matches the reference value
    double cur_vol_err = std::fabs(clip_res[1] - target_vol);
    if (cur_vol_err > vol_tol) {
      std::cerr << "LVIRA+ for cell " << cellID << ": given a maximum of  " << ims_tols_[0].max_num_iter <<
        " iteration(s) achieved error in volume for material " <<
        matID << " is " << cur_vol_err << ", volume tolerance is " << vol_tol << std::endl;
      throw std::runtime_error("Target error in volume exceeded, terminating...");
    }
  }

  /*!
    @brief Given a cell index, calculate the CellMatPoly using the LVIRA+ 
    interface reconstruction method.
    Uses nested dissections algorithm.
  */
  std::shared_ptr<CellMatPoly<Dim>> operator()(const int cell_op_ID) const {
    int cellID = icells_to_reconstruct[cell_op_ID];

    double dst_tol = ims_tols_[0].arg_eps;
    // Check if the cell is single-material
    if (cell_mat_ids_[cellID].size() == 1) {
      std::shared_ptr< CellMatPoly<Dim> > cmp_ptr(new CellMatPoly<Dim>(cellID));
      MatPoly<Dim> cell_matpoly;
      cell_get_matpoly(mesh_, cellID, &cell_matpoly, dst_tol);
      cell_matpoly.set_mat_id(cell_mat_ids_[cellID][0]);
      cmp_ptr->add_matpoly(cell_matpoly);

      return cmp_ptr;
    }

    // Use the nested dissections algorithm for multi-material cells.
    // Note that nested dissections uses this instance of reconstructor
    // to invoke get_plane_position position method. Nested dissections
    // itself does not have its own MeshWrapper, MatPoly_Clipper, etc.,
    // they all are reconstructor specific
    NestedDissections<LVIRAPlus, Dim, MatPoly_Splitter> 
      nested_dissections(*this, cellID, all_convex_);

    // Normally, we test all the permutations of the materials order
    // for cell with three or more materials
    bool enable_permutations = (cell_mat_ids_[cellID].size() != 2);

    // But if the stencil of a two-material cell includes cells with
    // a material not present in it, we still permute materials
    const std::vector<int>& ineighbor_cells = cell_neighbors_ids_[cellID];
    int nnc = static_cast<int>(ineighbor_cells.size());
    if (!enable_permutations) {
      for (int inc = 0; inc < nnc; inc++) {
        if (nghb_alien_mat_vfracs_[cellID][inc] > 0.0) {
          enable_permutations = true;
          break;
        }
      }
    }

    nested_dissections.set_cell_materials_order(enable_permutations);

    int npermutations = nested_dissections.num_materials_orders();
    Wonton::vector< std::vector< std::shared_ptr< CellMatPoly<Dim> > > > 
      permutations_cellmatpolys(npermutations);

    Wonton::transform(Wonton::make_counting_iterator(0),
                      Wonton::make_counting_iterator(npermutations),
                      permutations_cellmatpolys.begin(), nested_dissections);

    int iopt_permutation = 0;
    if (enable_permutations) {
      double min_neighbors_error = DBL_MAX;
      for (int iperm = 0; iperm < npermutations; iperm++) {
        const std::vector< std::shared_ptr< CellMatPoly<Dim> > >& cur_perm_cmp_ptrs = 
          permutations_cellmatpolys[iperm];
        const std::vector<int>& cell_matids = cell_mat_ids_[cellID];
        int nmats = static_cast<int>(cell_matids.size());

        double cur_error = 0.0;
        for (int inc = 0; inc < nnc; inc++) {
          int inghb_cell = ineighbor_cells[inc];
          double nghb_cell_vol = mesh_.cell_volume(inghb_cell);
          std::shared_ptr< CellMatPoly<Dim> > nghb_cmp_ptr = cur_perm_cmp_ptrs[inc + 1];
          int nghb_cmp_nmats = nghb_cmp_ptr->num_materials();

          for (int imat = 0; imat < nmats; imat++) {
            int mat_id = cell_matids[imat];

            int nnghb_mats = static_cast<int>(cell_mat_ids_[inghb_cell].size());
            //Find the reference volume fraction of the material in the neighboring cell
            int nghb_cell_imat = std::distance(cell_mat_ids_[inghb_cell].begin(),
              std::find(cell_mat_ids_[inghb_cell].begin(), 
                        cell_mat_ids_[inghb_cell].end(), mat_id));
            double nghb_reference_vfrac = 0.0;
            double nghb_weight = 1.0 - nghb_alien_mat_vfracs_[cellID][inc];
            if (nghb_cell_imat != nnghb_mats) {
              //Material should extend to the neighboring cell
              nghb_reference_vfrac = cell_mat_vfracs_[inghb_cell][nghb_cell_imat];

              if (lump_aliens_with_first_mat_) {
                //If the neighboring cell has materials NOT present in the CellMatPoly,
                //they were lumped with the FIRST clipped material
                if (mat_id == nghb_cmp_ptr->cell_matids()[0]) 
                  nghb_reference_vfrac += nghb_alien_mat_vfracs_[cellID][inc];
              }
              else {
                //If the neighboring cell has materials NOT present in the CellMatPoly,
                //they were lumped with the LAST clipped material
                if (mat_id == nghb_cmp_ptr->cell_matids()[nghb_cmp_nmats - 1])
                  nghb_reference_vfrac += nghb_alien_mat_vfracs_[cellID][inc];
              }
            }

            //Get the volume of the material in the CellMatPoly
            double nghb_cmp_vol = nghb_cmp_ptr->is_cell_material(mat_id) ? 
              nghb_cmp_ptr->material_moments(mat_id)[0] : 0.0;
            cur_error += pow2(nghb_weight*(nghb_reference_vfrac*nghb_cell_vol - nghb_cmp_vol));
          }
        }
        cur_error = sqrt(cur_error);

        if (cur_error < min_neighbors_error) {
          iopt_permutation = iperm;
          min_neighbors_error = cur_error;
        }
      }                       
    }
    const std::vector< std::shared_ptr< CellMatPoly<Dim> > >& opt_perm_cmp_ptrs = 
      permutations_cellmatpolys[iopt_permutation];
    return opt_perm_cmp_ptrs[0];
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
    double dst_tol = ims_tols_[0].arg_eps;
    cell_get_matpoly(mesh_, cellID, &mat_poly, dst_tol);

    return mat_poly;
  }

  /*!
    @brief Indices of neigboring cell that are split when errors are computed
    @param[in] cellID Cell index
    @return  Vector of indices of cell's neighbors through the nodes
  */
  std::vector<int> neighbor_cells_to_split(const int cellID) const {
    return cell_neighbors_ids_[cellID];
  }  

  /*!
    @brief IDs of the cell's faces to be associated with MatPoly's
    in the cell's decomposition
    @param[in] cellID Cell index
    @return  IDs of the cell's faces
  */
  std::vector<int> cell_face_group_ids(const int cellID,
                                       const bool faceted_faces) const { 
    std::vector<int> cfaces, cfdirs;
    mesh_.cell_get_faces_and_dirs(cellID, &cfaces, &cfdirs);
    
    if (!faceted_faces || Dim == 2)
      return cfaces;

    std::vector<int> facets_group_ids;
    for (int & cface : cfaces) {
      std::vector<int> fnodes;
      mesh_.face_get_nodes(cface, &fnodes);
      facets_group_ids.insert(facets_group_ids.end(), fnodes.size(), cface);
    }

    return facets_group_ids;
  }

private:
  /*!
    @brief For a given material, normal, and collection of MatPoly's
    finds the position of the plane corresponding to the material's
    volume fraction and computes the error in accordance with the
    modified LVIRA
    @param[in] cellID Index of the multi-material cell
    @param[in] cellMatID Local index of the clipped material
    @param[in] mixed_polys Vector of pointers to vectors of material poly's 
    that contain the material to clip and (possibly) other materials
    @param[in] plane_normal Direction of the cutting plane
    @return Error computed according to the modified LVIRA. Finaly a hard 
    number estimate for the sins sommited by your neighbors. Fund LVIRA++ to
    see the numbers for your in-laws.
  */  
  double neighbors_error(const int cellID, 
                         const int cellMatID,
                         const std::vector< std::vector< MatPoly<Dim> >* >& mixed_polys_ptrs,
                         const Vector<Dim>& plane_normal) const {
    Plane_t<Dim> cutting_plane;
    cutting_plane.normal = plane_normal;
    const std::vector< MatPoly<Dim> >& cell_mixed_polys = *(mixed_polys_ptrs[0]);

    double target_vol = cell_mat_vfracs_[cellID][cellMatID]*mesh_.cell_volume(cellID);
    double vol_tol = ims_tols_[0].fun_eps;

    // Confirm that the clipped volume is smaller than the volume of MatPoly's
    double mixed_polys_vol = 0.0;
    int nmixed_polys = static_cast<int>(cell_mixed_polys.size());
    for (int ipoly = 0; ipoly < nmixed_polys; ipoly++)
      mixed_polys_vol += cell_mixed_polys[ipoly].moments()[0];
    assert(target_vol < mixed_polys_vol + vol_tol);

    double min_mat_vfrac = 0.0;
    bool first_cut = false;
    if (lump_aliens_with_first_mat_) {
      // Determine if a material has already been chopped off
      min_mat_vfrac = *std::min_element(cell_mat_vfracs_[cellID].begin(), cell_mat_vfracs_[cellID].end());
      first_cut = (mixed_polys_vol > mesh_.cell_volume(cellID)*(1.0 - 0.5*min_mat_vfrac));
    }

    // Create cutting distance solver
    CuttingDistanceSolver<Dim, MatPoly_Clipper> 
      solve_cut_dst(cell_mixed_polys, cutting_plane.normal, ims_tols_[0], true);

    solve_cut_dst.set_target_volume(target_vol); 
    std::vector<double> clip_res = solve_cut_dst();

    // Check if the resulting volume matches the reference value
    double cur_vol_err = std::fabs(clip_res[1] - target_vol);
    if (cur_vol_err > vol_tol) {
      std::cerr << "LVIRA+ for cell " << cellID << ", testing normal ( ";
      for (int idim = 0; idim < Dim; idim++)
        std::cerr << plane_normal[idim] << " ";
      std::cerr << "): given a maximum of " <<
        ims_tols_[0].max_num_iter <<
        " iteration(s) achieved error in volume for material " << 
        cell_mat_ids_[cellID][cellMatID] << " is " << cur_vol_err << 
        ", volume tolerance is " << vol_tol << 
        ", volume of the split MatPoly's is " << mixed_polys_vol << 
        ", target volume is " << target_vol << std::endl;
      throw std::runtime_error("Target error in volume exceeded, terminating...");
    }

    cutting_plane.dist2origin = clip_res[0];

    const std::vector<int>& ineighbor_cells = cell_neighbors_ids_[cellID];
    int nnc = static_cast<int>(ineighbor_cells.size());
    assert(mixed_polys_ptrs.size() == nnc + 1);

    double nghb_error = 0.0;
    int matID = cell_mat_ids_[cellID][cellMatID];
    for (int inc = 0; inc < nnc; inc++) {
      const std::vector< MatPoly<Dim> >& nghb_mat_polys = *(mixed_polys_ptrs[inc + 1]);
      // Find the reference volume fraction of the material in the neighboring cell
      int inghb_cell = ineighbor_cells[inc];
      int nnghb_mats = static_cast<int>(cell_mat_ids_[inghb_cell].size());
      int nghb_cell_imat = std::distance(cell_mat_ids_[inghb_cell].begin(),
        std::find(cell_mat_ids_[inghb_cell].begin(), 
                  cell_mat_ids_[inghb_cell].end(), matID));
      double nghb_reference_vfrac = 0.0;
      // The weight corresponds to the fraction of materials that extend from the cell
      // to its neighbor
      double nghb_weight = 1.0 - nghb_alien_mat_vfracs_[cellID][inc];
      if (nghb_cell_imat != nnghb_mats) {
        // Clipped material extends to the neighboring cell
        nghb_reference_vfrac = cell_mat_vfracs_[inghb_cell][nghb_cell_imat];

        if (lump_aliens_with_first_mat_ && first_cut) {
          // If the neighboring cell has materials NOT present in the current cell,
          // then we assume they are in the same half-space as the material chopped
          // off on the first cut
          nghb_reference_vfrac += nghb_alien_mat_vfracs_[cellID][inc];
        }
      }

      //Find the volume of the material in the neighboring cell that
      //corresponds to the extended material interface
      MatPoly_Clipper clip_neighbor(vol_tol);
      clip_neighbor.set_matpolys(nghb_mat_polys, true);
      clip_neighbor.set_plane(cutting_plane);
      std::vector<double> nghb_moments = clip_neighbor();
      double nghb_IR_vol = nghb_moments[0];

      //Compute error
      nghb_error += pow2(nghb_weight*(nghb_reference_vfrac*mesh_.cell_volume(inghb_cell) - nghb_IR_vol));
    }
    
    return sqrt(nghb_error);
  }

  const Mesh_Wrapper& mesh_;
  const std::vector<IterativeMethodTolerances_t> ims_tols_;
  const bool all_convex_;
  std::vector< std::vector<int> > cell_mat_ids_;
  std::vector< std::vector<double> > cell_mat_vfracs_;
  std::vector<int> icells_to_reconstruct;
  std::vector< std::vector<int> > cell_neighbors_ids_;
  std::vector< std::vector<double> > nghb_alien_mat_vfracs_;
  bool lump_aliens_with_first_mat_;
};  // class LVIRAPlus

} // namespace Tangram

#endif    
