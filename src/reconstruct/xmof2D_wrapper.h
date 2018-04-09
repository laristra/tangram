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

#ifndef TANGRAM_xmof2D_h
#define TANGRAM_xmof2D_h

#include <memory>
#include "xmof2D.h"
#include "tangram/support/Point.h"
#include "tangram/driver/CellMatPoly.h"

/*!
 @file xmof2D_wrapper.h
 @brief Wrapper for the X-MOF 2D interface reconstruction method.
 Implements methods required by the standard Tangram driver.
 Requires XMOF2D library to be linked.
 */

namespace Tangram {
  
  /*!
   @class XMOF2D_Wrapper "xmof2D_wrapper.h"
   @brief Calculates material interfaces and constructs CellMatPoly objects
   for specified cells using the eXtended Moments-of-Fluid (X-MOF) method.

   @tparam Mesh_Wrapper A lightweight wrapper to a specific input mesh
   implementation that provides required functionality
   @tparam Dim The spatial dimension of the problem: has to be 2.
   */
  template <class Mesh_Wrapper, int Dim, class MatPoly_Splitter=void>
  class XMOF2D_Wrapper {
  public:
    /*!
     @brief Constructor initializing the X-MOF 2D reconstructor.
     @param[in] ir_tolerances User-specified tolerances for interface reconstruction in
     XMOF2D format
     */
    XMOF2D_Wrapper(const Mesh_Wrapper& Mesh,
                   const XMOF2D::IRTolerances* ir_tolerances = NULL) : mesh_(Mesh) {
      assert(Dim == 2);
      XMOF2D::MeshConfig mesh_cfg;
      int nnodes = mesh_.num_owned_nodes() + mesh_.num_ghost_nodes();
      mesh_cfg.nodes_coords.resize(nnodes);
      for (int inode = 0; inode < nnodes; inode++) {
        Point<Dim> cnode;
        mesh_.node_get_coordinates(inode, &cnode);
        mesh_cfg.nodes_coords[inode] = XMOF2D::Point2D(cnode[0], cnode[1]);
      }
      int nfaces = mesh_.num_owned_faces() + mesh_.num_ghost_faces();
      mesh_cfg.ifaces_nodes.resize(nfaces);
      for (int iface = 0; iface < nfaces; iface++)
        mesh_.face_get_nodes(iface, &mesh_cfg.ifaces_nodes[iface]);
      
      int ncells = mesh_.num_owned_cells() + mesh_.num_ghost_cells();
      mesh_cfg.icells_faces.resize(ncells);
      for (int icell = 0; icell < ncells; icell++) {
        std::vector<int> fdirs;
        mesh_.cell_get_faces_and_dirs(icell, &mesh_cfg.icells_faces[icell], &fdirs);
      }
      mesh_cfg.cells_material.clear();
      mesh_cfg.cells_material.resize(ncells, -1);
      
      XMOF2D::IRTolerances ir_tol;
      if (ir_tolerances)
        ir_tol = *ir_tolerances;
      else {
        ir_tol.dist_eps = 1.0e-15;
        ir_tol.div_eps = 1.0e-6;
        ir_tol.ddot_eps = 1.0e-14;
        ir_tol.vfrac_eps = 1.0e-14;
        ir_tol.ang_eps = 1.0e-13;
        ir_tol.mof_max_iter = 10000;
      }
     
      xmof_ir = std::make_shared<XMOF2D::XMOF_Reconstructor>(mesh_cfg, ir_tol);
    }
    /*!
     @brief Pass in the material data for use in the reconstruction.
     @param[in] cell_num_mats A vector specifying the
     number of materials in each cell of the mesh partition.
     @param[in] cell_mat_ids A vector of length (sum(cell_num_mats)) specifying
     the ID of each material in each cell
     @param[in] cell_mat_volfracs A vector of length(sum(cell_num_mats))
     specifying the volume fraction of each material in each cell.
     @param[in] cell_mat_centroids A vector of length(sum(cell_num_mats))
     specifying the centroids of each material in each cell.
     */
    void set_volume_fractions(std::vector<int> const& cell_num_mats,
                              std::vector<int> const& cell_mat_ids,
                              std::vector<double> const& cell_mat_volfracs,
                              std::vector<Point<Dim>> const& cell_mat_centroids) {
      assert(Dim == 2);
      int ncells = xmof_ir->get_base_mesh().ncells();
      assert(cell_num_mats.size() == ncells);
      XMOF2D::CellsMatData mat_data;
      mat_data.cells_materials.resize(ncells);
      mat_data.cells_vfracs.resize(ncells);
      mat_data.cells_centroids.resize(ncells);
      std::vector<int> offsets(cell_num_mats.size());
      for (int icell = 0; icell < ncells - 1; icell++)
        offsets[icell + 1] = offsets[icell] + cell_num_mats[icell];
      for (int icell = 0; icell < ncells; icell++) {
        int nmats = cell_num_mats[icell];
        mat_data.cells_materials[icell].insert(mat_data.cells_materials[icell].begin(),
          cell_mat_ids.begin() + offsets[icell],
          cell_mat_ids.begin() + offsets[icell] + nmats);
        mat_data.cells_vfracs[icell].insert(mat_data.cells_vfracs[icell].begin(),
          cell_mat_volfracs.begin() + offsets[icell],
          cell_mat_volfracs.begin() + offsets[icell] + nmats);
        for (int imat = 0; imat < nmats; imat++)
          mat_data.cells_centroids[icell].push_back(
            XMOF2D::Point2D(cell_mat_centroids[offsets[icell] + imat][0],
                            cell_mat_centroids[offsets[icell] + imat][1]));
      }
      xmof_ir->set_materials_data(mat_data);
    }
    
    /*!
     @brief Pass in the local partition indices of cells for which CellMatPoly objects 
     are to be constructed. If the index is in the list, a CellMatPoly object will be
     created even for a single-material cell.
     @param[in] cellIDs_to_op_on A vector of length up to (num_local_partition_cells) 
     specifying the local indices of cells for which CellMatPoly objects are requested.
     */
    void set_cell_indices_to_operate_on(std::vector<int> const& cellIDs_to_op_on) {
      icells_to_reconstruct = cellIDs_to_op_on;
    }

    /*!
     @brief Given a cell index, calculate the CellMatPoly based on the material data
     */
    std::shared_ptr<CellMatPoly<Dim>> operator()(const int cell_op_ID) const {
      assert(Dim == 2);
      assert(cell_op_ID < icells_to_reconstruct.size());
      int cellID = icells_to_reconstruct[cell_op_ID];
      
      CellMatPoly<Dim>* cell_mat_poly = new CellMatPoly<Dim>(cellID);
      const XMOF2D::BaseMesh& mesh = xmof_ir->get_base_mesh();
      if (xmof_ir->get_cell_materials(cellID).size() > 1) {
        xmof_ir->construct_minimesh(cellID);
        const XMOF2D::MiniMesh& cminimesh = mesh.get_cell(cellID).get_minimesh();
        for (int isc = 0; isc < cminimesh.ncells(); isc++) {
          const XMOF2D::Cell& subcell = cminimesh.get_cell(isc);
          int nvrts = subcell.nfaces();
          std::vector<Tangram::Point<Dim>> subcell_vrts(nvrts);
          std::vector<Tangram::Entity_kind> vrts_parentkind(nvrts);
          std::vector<int> vrts_iparent(nvrts);
          for (int ivrt = 0; ivrt < nvrts; ivrt++) {
            subcell_vrts[ivrt] = Tangram::Point<Dim>(subcell.get_node_crd(ivrt).x,
                                                     subcell.get_node_crd(ivrt).y);
            vrts_iparent[ivrt] = subcell.get_node(ivrt).iparent();
            vrts_parentkind[ivrt] = (vrts_iparent[ivrt] == -1) ?
            Tangram::Entity_kind::UNKNOWN_KIND : Tangram::Entity_kind::NODE;
          }
          std::vector<Tangram::Entity_kind> sides_parentkind(nvrts);
          std::vector<int> sides_iparent(nvrts);
          for (int iside = 0; iside < nvrts; iside++) {
            sides_iparent[iside] = subcell.get_face(iside).iparent();
            sides_parentkind[iside] = (sides_iparent[iside] == -1) ?
            Tangram::Entity_kind::UNKNOWN_KIND : Tangram::Entity_kind::FACE;
            
          }
          (*cell_mat_poly).add_matpoly(subcell.get_material_index(), nvrts,
                                       &subcell_vrts[0],
                                       &vrts_parentkind[0], &vrts_iparent[0],
                                       &sides_parentkind[0], &sides_iparent[0]);
        }
      }
      else {
        const XMOF2D::Cell& ccell = mesh.get_cell(cellID);
        int nvrts = ccell.nfaces();
        std::vector<Tangram::Point<Dim>> cell_vrts(nvrts);
        for (int ivrt = 0; ivrt < nvrts; ivrt++)
          cell_vrts[ivrt] = Tangram::Point<Dim>(ccell.get_node_crd(ivrt).x,
                                                ccell.get_node_crd(ivrt).y);
        (*cell_mat_poly).add_matpoly(xmof_ir->get_cell_materials(cellID)[0],
                                     nvrts, &cell_vrts[0],
                                     NULL, NULL, NULL, NULL);
      }
    
      return std::shared_ptr<CellMatPoly<Dim>>(cell_mat_poly);
    }

  private:
    const Mesh_Wrapper& mesh_; // Provided base mesh wrapper
    std::shared_ptr<XMOF2D::XMOF_Reconstructor> xmof_ir; // XMOF2D reconstructor object
    std::vector<int> icells_to_reconstruct; // List of cells to create CellMatPoly objects for
  }; // class XMOF2D_Wrapper
}  // namespace Tangram

#endif /* TANGRAM_xmof2D_h */
