/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#ifndef MOF_H
#define MOF_H

#include <vector>
#include <functional>
#include "tangram/support/tangram.h"
#include "tangram/support/MatPoly.h"
#include "tangram/reconstruct/nested_dissections.h"
#include "tangram/reconstruct/cutting_distance_solver.h"
#include "tangram/support/bfgs.h"

/*!
 @file MOF.h
  @brief Calculates the interfaces and constructs CellMatPoly using the MOF
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

Vector<1> cartesian_to_polar(const Vector<2> cartesian_vec) {
  return Vector<1>(atan2(cartesian_vec[1], cartesian_vec[0]));
}

Vector<2> cartesian_to_polar(const Vector<3> cartesian_vec) {
  double len = cartesian_vec.norm();
  return Vector<2>(atan2(cartesian_vec[1], cartesian_vec[0]), 
                   acos(cartesian_vec[2]/len));
}

Vector<2> polar_to_cartesian(const Vector<1> polar_vec) {
  return Vector<2>(cos(polar_vec[0]), sin(polar_vec[0]));
}

Vector<3> polar_to_cartesian(const Vector<2> polar_vec) {
  return Vector<3>(sin(polar_vec[1])*cos(polar_vec[0]),
                   sin(polar_vec[1])*sin(polar_vec[0]),
                   cos(polar_vec[1])); 
}

template <class Mesh_Wrapper, int Dim, class MatPoly_Splitter, class MatPoly_Clipper>
class MOF {
public:
  /*!
    @brief Constructor for a MOF interface reconstruction algorithm
    @param[in] Mesh A lightweight wrapper to a specific input mesh
    implementation that provides certain functionality
    @param[in] im_tols Tolerances for iterative methods
    @param[in] all_convex Flag indicating whether all mesh cells are convex
  */
  explicit MOF(const Mesh_Wrapper& Mesh, 
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
                              const& cell_mat_centroids) {
    cell_mat_ids_.clear();
    cell_mat_vfracs_.clear();       
    cell_mat_centroids_.clear();                               
    int ncells = mesh_.num_owned_cells() + mesh_.num_ghost_cells();
    cell_mat_ids_.resize(ncells);
    cell_mat_vfracs_.resize(ncells);  
    cell_mat_centroids_.resize(ncells);
    int offset = 0;
    for (int icell = 0; icell < ncells; icell++) {
      double cell_volume = mesh_.cell_volume(icell);  
      int nmats = cell_num_mats[icell];

      for (int icmat = 0; icmat < nmats; icmat++) {
        double mat_volume = cell_volume*cell_mat_volfracs[offset + icmat];
        if (mat_volume > im_tols_.fun_eps) {
          cell_mat_ids_[icell].push_back(cell_mat_ids[offset + icmat]);
          cell_mat_vfracs_[icell].push_back(cell_mat_volfracs[offset + icmat]);
          cell_mat_centroids_[icell].push_back(cell_mat_centroids[offset + icmat]);
        }
      }
      offset += nmats;
    }
  }
  
  /*!
    @brief Used iterative methods tolerances
    @return  Tolerances for iterative methods, 
    here im_tols_.fun_eps is the volume tolerance and
    im_tols_.arg_eps is the tolerance on norm of the change 
    in polar angles defining the normal to the cutting plane
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
    assert(Dim > 1);

    int cellMatID = std::distance(cell_mat_ids_[cellID].begin(),
      std::find(cell_mat_ids_[cellID].begin(), 
                cell_mat_ids_[cellID].end(), matID));

    Point<Dim> cell_cen;
    mesh_.cell_centroid(cellID, &cell_cen);
    Vector<Dim - 1> init_guess = cartesian_to_polar(
      cell_cen - cell_mat_centroids_[cellID][cellMatID]);

    std::function<double(const Vector<Dim - 1>&)> bfgs_obj_fun = 
      [this, &cellID, &cellMatID, &mixed_polys]
      (const Vector<Dim - 1>& cur_arg)->double {
      return centroid_error(cellID, cellMatID, mixed_polys, 
                            polar_to_cartesian(cur_arg));
    };    

    Vector<Dim - 1> ang_min = 
      bfgs<Dim - 1>(bfgs_obj_fun, init_guess, im_tols_.max_num_iter,
                    im_tols_.arg_eps, im_tols_.fun_eps);

    cutting_plane.normal = polar_to_cartesian(ang_min);

    CuttingDistanceSolver<Dim, MatPoly_Clipper> 
      solve_cut_dst(mixed_polys, cutting_plane.normal, im_tols_, planar_faces);

    double target_vol = cell_mat_vfracs_[cellID][cellMatID]*mesh_.cell_volume(cellID);
    solve_cut_dst.set_target_volume(target_vol); 
    std::vector<double> clip_res = solve_cut_dst();
    cutting_plane.dist2origin = clip_res[0];
  }

  /*!
    @brief Given a cell index, calculate the CellMatPoly using the MOF 
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
    NestedDissections<MOF, Dim, MatPoly_Splitter> 
      nested_dissections(*this, cellID, all_convex_);

    // We test all the permutations of the materials order
    nested_dissections.set_cell_materials_order(true);

    int npermutations = nested_dissections.num_materials_orders();
    Tangram::vector<std::shared_ptr<CellMatPoly<Dim>>> 
      permutations_cellmatpoly(npermutations);

    Tangram::transform(make_counting_iterator(0),
                       make_counting_iterator(npermutations),
                       permutations_cellmatpoly.begin(), nested_dissections);

    int nmats = static_cast<int>(cell_mat_ids_[cellID].size());

    int iopt_permutation = -1;
    double min_centroids_error = DBL_MAX;
    for (int iperm = 0; iperm < npermutations; iperm++) {
      std::shared_ptr<CellMatPoly<Dim>> cur_cmp_ptr = 
        permutations_cellmatpoly[iperm];

      double cur_error = 0.0;      
      for (int imat = 0; imat < nmats; imat++) {
        const std::vector<double>& cur_moments = 
          cur_cmp_ptr->material_moments(cell_mat_ids_[cellID][imat]);

        Point<Dim> mat_centroid;
        for (int idim = 0; idim < Dim; idim++)
          mat_centroid[idim] = cur_moments[idim + 1]/cur_moments[0];

        cur_error += (mat_centroid - cell_mat_centroids_[cellID][imat]).norm(false);
      }
      cur_error = sqrt(cur_error);

      if (cur_error < min_centroids_error) {
        iopt_permutation = iperm;
        min_centroids_error = cur_error;
      }
    }                       

    return permutations_cellmatpoly[iopt_permutation];
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
  /*!
    @brief For a given material, normal, and collection of MatPoly's
    finds the position of the plane corresponding to the material's
    volume fraction and computes the distance between the reference
    and the actual aggregated centroids of clipped MatPoly's
    @param[in] cellID Index of the multi-material cell
    @param[in] cellMatID Local index of the clipped material
    @param[in] mixed_polys Vector of material poly's that contain the material
    to clip and (possibly) other materials
    @param[in] plane_normal Direction of the cutting plane
    @return Distance between the reference and the actual centroids
  */  
  double centroid_error(const int cellID, 
                        const int cellMatID,
                        const std::vector< MatPoly<Dim> >& mixed_polys,
                        const Vector<Dim>& plane_normal) const {
    double vol_tol = im_tols_.fun_eps;

    Plane_t<Dim> cutting_plane;
    cutting_plane.normal = plane_normal;

    double target_vol = cell_mat_vfracs_[cellID][cellMatID]*mesh_.cell_volume(cellID);

#ifdef DEBUG
    //Confirm that the clipped volume is smaller than the volume of MatPoly's
    double mixed_polys_vol = 0.0;
    for (int ipoly = 0; ipoly < mixed_polys.size(); ipoly++)
      mixed_polys_vol += mixed_polys[ipoly].moments()[0];
    assert(target_vol < mixed_polys_vol + vol_tol);
#endif

    //Create cutting distance solver
    CuttingDistanceSolver<Dim, MatPoly_Clipper> 
      solve_cut_dst(mixed_polys, cutting_plane.normal, im_tols_, true);

    solve_cut_dst.set_target_volume(target_vol); 
    std::vector<double> clip_res = solve_cut_dst();

#ifdef DEBUG
    // Check if the resulting volume matches the reference value
    double cur_vol_err = std::fabs(clip_res[1] - target_vol);
    if (cur_vol_err > vol_tol) {
      std::cerr << "MOF for cell " << cellID << ", testing normal ( ";
      for (int idim = 0; idim < Dim; idim++)
        std::cerr << plane_normal[idim] << " ";
      std::cerr << "): after " <<
        im_tols_.max_num_iter <<
        " iteration(s) achieved error in volume for material " << 
        cell_mat_ids_[cellID][cellMatID] << " is " << cur_vol_err << 
        ", volume tolerance is " << vol_tol << 
        ", volume of the split MatPoly's is " << mixed_polys_vol << 
        ", target volume is " << target_vol << std::endl;
    }
#endif    

    Point<Dim> mat_centroid;
    for (int idim = 0; idim < Dim; idim++)
      mat_centroid[idim] = clip_res[idim + 2]/clip_res[1];
    
    return (mat_centroid - cell_mat_centroids_[cellID][cellMatID]).norm();
  }

  const Mesh_Wrapper& mesh_;
  const IterativeMethodTolerances_t im_tols_;
  const bool all_convex_;
  std::vector< std::vector<int> > cell_mat_ids_;
  std::vector< std::vector<double> > cell_mat_vfracs_;
  std::vector< std::vector< Point<Dim> > > cell_mat_centroids_;
  std::vector<int> icells_to_reconstruct;
};  // class MOF

} // namespace Tangram

#endif    