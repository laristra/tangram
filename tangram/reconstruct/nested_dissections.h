/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#ifndef NESTED_DISSECTIONS_H
#define NESTED_DISSECTIONS_H

#include <vector>
#include <numeric>
#include <algorithm>
#include "tangram/support/tangram.h"
#include "tangram/support/MatPoly.h"


/*!
 @file nested_dissections.h
  @brief Implements the nested dissections algorithm. 
  Each instance of this class is associated with a specific CellID. 
  If is possible to sequentially use the same instance of this class 
  for different mesh cells because cellID is a reference to the
  variable that can be modified externally. Materials order should
  generally be updated if cellID is modified.
  
  @tparam Reconstructor A reconstructor that implements an algorithm 
  for finding the position of a plane between a given material and the rest
  @tparam Dim The spatial dimension of the problem
  @tparam MatPoly_Splitter An operator for splitting a vector of MatPoly's 
  into half-space sets with a cutting plane
*/

namespace Tangram {

template <class Reconstructor, int Dim, class MatPoly_Splitter>
class NestedDissections {
public:
  /*!
    @brief Constructor for the nested dissections algorithm.
    @param[in] rec A reconstructor capable of finding the position of
    a cutting plane
    @param[in] cellID Index of a cell for which interface reconstruction 
    is performed
    @param[in] convex_cell Flag indicating whether the cell is convex
  */
  explicit NestedDissections(const Reconstructor& rec,
                             const int& cellID,
                             const bool convex_cell) : 
                             reconstructor_(rec), cell_id_(cellID),
                             convex_cell_(convex_cell) {}

  /*!
    @brief Specifies the order in which materials are clipped
    @param[in] cell_materials_order Vector of local material indices, if a cell 
    has three materials, it will have three entries (0, 1, 2) in certain order
  */
  void set_cell_materials_order(const std::vector<int>& cell_materials_order) {
    cell_materials_order_.clear();                              
    cell_materials_order_.push_back(cell_materials_order);
  }

  /*!
    @brief Sets materials clipping order: either all ordering permutations are
    enabled, or materials are clipped in the order they are listed for the cell
    @param[in] enable_permutations Flag indicating whether the order permutations
    are to be used. If enabled, ID of the permutation can be used as the parameter
    of the operator
  */
  void set_cell_materials_order(const bool enable_permutations) {
    int nmats = static_cast<int>(reconstructor_.cell_materials(cell_id_).size());
    cell_materials_order_.resize(1);
    cell_materials_order_[0].resize(nmats);
    std::iota(cell_materials_order_[0].begin(), cell_materials_order_[0].end(), 0);                              

    if (enable_permutations) {
      std::vector<int> cur_mat_order = cell_materials_order_[0];
      while (std::next_permutation(cur_mat_order.begin(), cur_mat_order.end()))
        cell_materials_order_.push_back(cur_mat_order);
    }
  }

  /*!
    @brief The number of materials orders for which the respective CellMatPoly's
    can be constructed. If permuations are enabled, it is the number of possible
    permutations, otherwise it is equal to 1.
    @return  Number of materials orders, their index can be used as a parameter
    for this class's operator.
  */
  int num_materials_orders() { return static_cast<int>(cell_materials_order_.size()); }

  /*!
    @brief Calculates the CellMatPoly for the given interface reconstruction method
    using the nested dissections algorithm. If materials order permutations are enabled,
    permutation ID can be used to choose an alternative ordering. 
  */
  std::shared_ptr< CellMatPoly<Dim> > operator()(const int permutation_ID = 0) const {
    assert(unsigned(permutation_ID) < cell_materials_order_.size());

    // Get material indices for the cell from the reconstructor
    const std::vector<int>& mat_ids = reconstructor_.cell_materials(cell_id_);    
    int nmats = static_cast<int>(mat_ids.size());

    std::shared_ptr< CellMatPoly<Dim> > cmp_ptr = 
      std::make_shared< CellMatPoly<Dim> >(cell_id_);

    //MatPoly corresponding to a non-convex cell can have non-planar faces
    MatPoly<Dim> cell_matpoly;
    if (convex_cell_)
      cell_matpoly = reconstructor_.cell_matpoly(cell_id_);
    else
      reconstructor_.cell_matpoly(cell_id_).faceted_matpoly(&cell_matpoly);
    
    double vol_tol = reconstructor_.iterative_methods_tolerances()[0].fun_eps;
    double dst_tol = reconstructor_.iterative_methods_tolerances()[0].arg_eps;

    Plane_t<Dim> cutting_plane {};
    HalfSpaceSets_t<Dim> hs_sets; // Structure containing vectors of MatPoly's and their
                                  // aggregated moments for the respective half-spaces

    //Original vector of mixed MatPoly's consists of the cell's MatPoly
    hs_sets.upper_halfspace_set.matpolys = { cell_matpoly };

    //Create Splitter instance: assumes split MatPoly's are convex
    MatPoly_Splitter split_matpolys(hs_sets.upper_halfspace_set.matpolys, 
                                    cutting_plane, vol_tol, dst_tol, true);

    for (int imat = 0; imat < nmats; imat++) {
      int matid = mat_ids[cell_materials_order_[permutation_ID][imat]];
      
      MatPolySet_t<Dim>* single_mat_set_ptr;
      //On the last iteration the remaining part is single-material,
      //so we don't need to split it
      if (imat == nmats - 1)
        single_mat_set_ptr = &hs_sets.upper_halfspace_set;
      else {
        // Find cutting plane position: assumes all faces are planar
        // Note that this method generally does NOT require MatPoly_Splitter,
        // but uses MatPoly_Clipper instead. Nested dissections itself
        // has neither MatPoly_Clipper, nor MeshWrapper, those are
        // handled by the reconstructor that created this nested dissections
        // instance.
        reconstructor_.get_plane_position(cell_id_, matid, 
                                          hs_sets.upper_halfspace_set.matpolys, 
                                          cutting_plane, true);

        // On the first cut we have to decompose a non-convex MatPoly
        // This is done inside the loop to make the first step cheaper: 
        // finding the position of the cutting plane normally requires 
        // only computation of moments, which we assume to be possible 
        // without decomposing into tetrahedrons (e.g. if r3d is used). 
        if (not convex_cell_ and imat == 0) {
          hs_sets.upper_halfspace_set.matpolys.clear();
          std::vector<int> face_group_ids = reconstructor_.cell_face_group_ids(cell_id_, true);
          cell_matpoly.decompose(hs_sets.upper_halfspace_set.matpolys, &face_group_ids);
        }

        hs_sets = split_matpolys();
        //Single-material MatPoly's are below the plane
        single_mat_set_ptr = &hs_sets.lower_halfspace_set;

        // Filter out MatPoly's with volumes below tolerance
        int ismp = 0;
        int nb_single_mat_points = single_mat_set_ptr->matpolys.size();
        while (ismp < nb_single_mat_points) {
          if (single_mat_set_ptr->matpolys[ismp].moments()[0] >= vol_tol) {
            ismp++;
          } else {
            for (int im = 0; im < Dim + 1; im++) {
              single_mat_set_ptr->moments[im] -=
                single_mat_set_ptr->matpolys[ismp].moments()[im];
            }

            auto deleted = single_mat_set_ptr->matpolys.begin() + ismp;
            single_mat_set_ptr->matpolys.erase(deleted);
            nb_single_mat_points = single_mat_set_ptr->matpolys.size();
          }
        }
      }

      // Add single-material poly's below the plane to CellMatPoly
      for (auto& cur_mat_poly : single_mat_set_ptr->matpolys) {
        cur_mat_poly.set_mat_id(matid);
        cmp_ptr->add_matpoly(cur_mat_poly);
      }

      if (not single_mat_set_ptr->matpolys.empty())
        cmp_ptr->assign_material_moments(matid, single_mat_set_ptr->moments);
    }

    return cmp_ptr;
  }

private:
  const Reconstructor& reconstructor_;
  const int& cell_id_;
  const bool convex_cell_;
  std::vector< std::vector<int> > cell_materials_order_;
}; // class NestedDissections

} // namespace Tangram

#endif
 
