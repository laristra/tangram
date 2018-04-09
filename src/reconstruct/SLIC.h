/*
 Copyright (c) 2017, Los Alamos National Security, LLC
 All rights reserved.
 
 Copyright 2017. Los Alamos National Security, LLC. This software was produced
 under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National
 Laboratory (LANL), which is operated by Los Alamos National Security, LLC for
 the U.S. Department of Energy. The U.S. Government has rights to use,
 reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS
 NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
 LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
 derivative works, such modified software should be clearly marked, so as not to
 confuse it with the version available from LANL.
 
 Additionally, redistribution and use in source and binary forms, with or
 without modification, are permitted provided that the following conditions are
 met:
 
 1. Redistributions of source code must retain the above copyright notice,
 this list of conditions and the following disclaimer.
 2. Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.
 3. Neither the name of Los Alamos National Security, LLC, Los Alamos
 National Laboratory, LANL, the U.S. Government, nor the names of its
 contributors may be used to endorse or promote products derived from this
 software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND
 CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL
 SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
 IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef TANGRAM_SLIC_H_
#define TANGRAM_SLIC_H_

#include <vector>
#include <algorithm>
#include <numeric>
#include <cmath>
#include <memory>
#include <float.h>

#include "tangram/support/Point.h"
#include "tangram/driver/CellMatPoly.h"
#include "tangram/support/MatPoly.h"
#include "tangram/intersect/split_r3d.h"

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

  template <class Mesh_Wrapper, int Dim, class MatPoly_Splitter>
  class SLIC {
  public:
    /*!
     @brief Constructor performing a SLIC algorithm for interface reconstruction.
     */
    explicit SLIC(const Mesh_Wrapper & Mesh) : mesh_(Mesh) {
      // For now
      assert(Dim == 3);
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
      
      auto nc = mesh_.num_entities(Entity_kind::CELL);
      cell_mat_offsets_.resize(nc);
      cell_mat_offsets_[0] = 0;
      for (int c(1); c < nc; ++c)
        cell_mat_offsets_[c] = cell_mat_offsets_[c-1] + cell_num_mats_[c-1];
    }
    
    void set_cell_indices_to_operate_on(std::vector<int> const& cellIDs_to_op_on) {
      icells_to_reconstruct = cellIDs_to_op_on;
    }
    
    /*!
     @brief Given a cell index, calculate the CellMatPoly for this reconstruction
     */
    std::shared_ptr<CellMatPoly<Dim>> operator()(const int cell_op_ID) const {
      double vfrac_tol = 1.0e-13;

      int cellID = icells_to_reconstruct[cell_op_ID];
      auto numMats = cell_num_mats_[cellID];
      
      CellMatPoly<Dim>* cellpoly = new CellMatPoly<Dim>(cellID);
      
      // We are assuming here - for now - that cells are rectangular prisms for
      // simplicity
      
      auto iStart = cell_mat_offsets_[cellID];
      
      //Sets of MatPoly's on two sides of the cutting plane
      HalfSpaceSets_t<3> hs_sets;

      // For every chopped off material, we split MatPoly's above the plane
      hs_sets.upper_halfspace_set.matpolys.resize(1);
      mesh_.cell_get_matpoly(cellID, &hs_sets.upper_halfspace_set.matpolys[0]);

      // Just going along x-direction - again assuming rectangular prisms only
      Plane_t cutting_plane;
      cutting_plane.normal = Vector3(1.0, 0.0, 0.0);

      //Create Splitter instance
      MatPoly_Splitter split_matpolys(hs_sets.upper_halfspace_set.matpolys, 
                                      cutting_plane, true);

      //Cutting from left to right: find the x range
      double xmin = DBL_MAX, xmax = DBL_MIN;
      for (int immp = 0; immp < hs_sets.upper_halfspace_set.matpolys.size(); immp++) {
        const MatPoly<3>& cur_matpoly = hs_sets.upper_halfspace_set.matpolys[immp];
        for (int ivrt = 0; ivrt < cur_matpoly.num_vertices(); ivrt++) {
          double pt_x = cur_matpoly.vertex_point(ivrt)[0];
          if (pt_x < xmin) xmin = pt_x;
          if (pt_x > xmax) xmax = pt_x;
        }
      }

      double cell_volume = mesh_.cell_volume(cellID);
      double x_cut = xmin;
      for (int iMat(0); iMat < numMats; iMat++) {
        auto vfrac = cell_mat_volfracs_[iStart + iMat];
        // If the mass fraction is too small, skip it
        if (vfrac < vfrac_tol) continue;
          
        const MatPolySet_t<3>* single_mat_set_ptr;
        //On the last iteration the remaining part is single-material,
        //so we don't need to split it
        if (iMat == numMats - 1)
          single_mat_set_ptr = &hs_sets.upper_halfspace_set;
        else {
          // Find distance to origin for the cutting plane
          x_cut += vfrac*(xmax - xmin);
          cutting_plane.dist2origin = -x_cut;

          hs_sets = split_matpolys();

          // Confirm that the volume fraction matches the reference value
          assert(std::fabs(hs_sets.lower_halfspace_set.moments[0]/cell_volume - vfrac) < vfrac_tol);

          //Chopped off single-material MatPoly's are below the plane
          single_mat_set_ptr = &hs_sets.lower_halfspace_set;
        }
        //Add single-material MatPoly's to CellMatPoly
        for (int ismp = 0; ismp < single_mat_set_ptr->matpolys.size(); ismp++) {
          const MatPoly<3>& cur_matpoly = single_mat_set_ptr->matpolys[ismp];
          // Flatten the face vertices for the single-material MatPoly
          int nfaces = cur_matpoly.num_faces();
          std::vector<int> nface_vrts(nfaces);
          std::vector<int> faces_vrts;
          for (int iface = 0; iface < nfaces; iface++) {
            const std::vector<int>& face_ivrts = cur_matpoly.face_vertices(iface);
            nface_vrts[iface] = face_ivrts.size();
            faces_vrts.insert(faces_vrts.end(), face_ivrts.begin(), face_ivrts.end());
          }
          //Add the MatPoly below the cutting plane to CellMatPoly
          (*cellpoly).add_matpoly(cell_mat_ids_[iStart + iMat], 
                                  cur_matpoly.num_vertices(), 
                                  &cur_matpoly.points()[0],
                                  nullptr, nullptr,
                                  nfaces, &nface_vrts[0], &faces_vrts[0],
                                  nullptr, nullptr);
        }
      }
      
      return std::shared_ptr<CellMatPoly<Dim>>(cellpoly);
    }
  private:
    const Mesh_Wrapper & mesh_;
    std::vector<int> cell_num_mats_;
    std::vector<int> cell_mat_ids_;
    std::vector<double> cell_mat_volfracs_;
    std::vector<int> cell_mat_offsets_;
    std::vector<int> icells_to_reconstruct;
  };  // class SLIC
}  // namespace Tangram


#endif  // TANGRAM_SLIC_H_
