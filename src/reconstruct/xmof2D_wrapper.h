//
//  xmof2D.h
//  X-MOF-2D
//
//  Created by Evgeny Kikinzon.
//  Copyright © 2017 Los Alamos National Laboratory. All rights reserved.
//

#ifndef xmof2D_h
#define xmof2D_h

#include "xmof2D.h"

#include "tangram/support/Point.h"
#include "tangram/driver/CellMatPoly.h"

namespace Tangram {
  
  template <class Mesh_Wrapper, int Dim>
  class XMOF2D_Wrapper {
  public:
    XMOF2D_Wrapper(const Mesh_Wrapper& Mesh,
                   const XMOF2D::IRTolerances* ir_tolerances = NULL) : mesh_(Mesh) {
      assert(Dim == 2);
      XMOF2D::MeshConfig mesh_cfg;
      XMOF2D::GlobalIndData gid_data;
      int nnodes = mesh_.num_owned_nodes() + mesh_.num_ghost_nodes();
      mesh_cfg.nodes_coords.resize(nnodes);
      gid_data.nodes_global_ind.resize(nnodes);
      for (int inode = 0; inode < nnodes; inode++) {
        Point<Dim> cnode;
        mesh_.node_get_coordinates(inode, &cnode);
        mesh_cfg.nodes_coords[inode] = XMOF2D::Point2D(cnode[0], cnode[1]);
        gid_data.nodes_global_ind[inode] = mesh_.get_global_id(inode, Tangram::Entity_kind::NODE);
      }
      int nfaces = mesh_.num_owned_faces() + mesh_.num_ghost_faces();
      mesh_cfg.ifaces_nodes.resize(nfaces);
      gid_data.faces_global_ind.resize(nfaces);
      for (int iface = 0; iface < nfaces; iface++) {
        mesh_.face_get_nodes(iface, &mesh_cfg.ifaces_nodes[iface]);
        gid_data.faces_global_ind[iface] = mesh_.get_global_id(iface, Tangram::Entity_kind::FACE);
      }
      
      int ncells = mesh_.num_owned_cells() + mesh_.num_ghost_cells();
      mesh_cfg.icells_faces.resize(ncells);
      mesh_cfg.ifaces_out_cells.resize(nfaces, -1);
      gid_data.cells_global_ind.resize(ncells);
      for (int icell = 0; icell < ncells; icell++) {
        std::vector<int> fdirs;
        mesh_.cell_get_faces_and_dirs(icell, &mesh_cfg.icells_faces[icell], &fdirs);
        for (int iside = 0; iside < fdirs.size(); iside++) {
          int iface = mesh_cfg.icells_faces[icell][iside];
          if (fdirs[iside] == 1)
            mesh_cfg.ifaces_out_cells[iface] = icell;
        }
        gid_data.cells_global_ind[icell] = mesh_.get_global_id(icell, Tangram::Entity_kind::CELL);
      }
      for (int iface = 0; iface < nfaces; iface++) {
        if (mesh_cfg.ifaces_out_cells[iface] == -1) {
          std::vector<int> fcells;
          mesh_.face_get_cells(iface, &fcells);
          assert(fcells.size() == 1);
          mesh_cfg.ifaces_out_cells[iface] = fcells[0];
          std::reverse(mesh_cfg.ifaces_nodes[iface].begin(), mesh_cfg.ifaces_nodes[iface].end());
        }
      }
      mesh_cfg.cells_material.clear();
      mesh_cfg.cells_material.resize(ncells, -1);
      
      XMOF2D::IRTolerances ir_tol;
      if (ir_tolerances)
        ir_tol = *ir_tolerances;
      else {
        ir_tol.dist_eps = 1.0e-15;
        ir_tol.div_eps = 1.0e-6;
        ir_tol.ddot_eps = 1.0e-12;
        ir_tol.vfrac_eps = 1.0e-13;
        ir_tol.ang_eps = 1.0e-12;
        ir_tol.mof_max_iter = 1e3;
      }
     
      xmof_ir = std::make_shared<XMOF2D::XMOF_Reconstructor>(mesh_cfg, ir_tol);
      xmof_ir->set_global_ind(gid_data);
    }
    
    void set_volume_fractions(std::vector<int> const& cell_num_mats,
                              std::vector<int> const& cell_mat_ids,
                              std::vector<double> const& cell_mat_volfracs,
                              std::vector<Point<Dim>> const& cell_mat_centroids) {
      assert(Dim == 2);
      int ncells = xmof_ir->get_base_mesh().ncells();
      assert(cell_num_mats.size() >= ncells);
      XMOF2D::CellsMatData mat_data;
      mat_data.cells_materials.resize(ncells);
      mat_data.cells_vfracs.resize(ncells);
      mat_data.cells_centroids.resize(ncells);
      std::vector<int> offsets(cell_num_mats.size());
      for (int cgi = 0; cgi < cell_num_mats.size() - 1; cgi++)
        offsets[cgi + 1] = offsets[cgi] + cell_num_mats[cgi];
      for (int cellID = 0; cellID < ncells; cellID++) {
        int cell_gid = xmof_ir->cell_global_ind(cellID);
        int nmats = cell_num_mats[cell_gid];
        mat_data.cells_materials[cellID].resize(nmats);
        mat_data.cells_vfracs[cellID].resize(nmats);
        std::copy(cell_mat_ids.begin() + offsets[cell_gid],
                  cell_mat_ids.begin() + offsets[cell_gid] + nmats,
                  mat_data.cells_materials[cellID].begin());
        std::copy(cell_mat_volfracs.begin() + offsets[cell_gid],
                  cell_mat_volfracs.begin() + offsets[cell_gid] + nmats,
                  mat_data.cells_vfracs[cellID].begin());
        for (int imat = 0; imat < nmats; imat++)
          mat_data.cells_centroids[cellID].push_back(
            XMOF2D::Point2D(cell_mat_centroids[offsets[cell_gid] + imat][0],
                            cell_mat_centroids[offsets[cell_gid] + imat][1]));
      }
      xmof_ir->set_materials_data(mat_data);
    }
    
    void set_cell_indices_to_operate_on(std::vector<int> const& cellIDs_to_op_on) {
      icells_to_reconstruct = cellIDs_to_op_on;
    }

    std::shared_ptr<CellMatPoly<Dim>> operator()(const int cell_op_ID) const {
      assert(Dim == 2);
      assert(cell_op_ID < icells_to_reconstruct.size());
      int cellID = icells_to_reconstruct[cell_op_ID];
      
//      if (xmof_ir->get_cell_materials(cellID).size() == 1) {
//        return std::shared_ptr<CellMatPoly<Dim>>(nullptr);
      
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
    const Mesh_Wrapper& mesh_;
    std::shared_ptr<XMOF2D::XMOF_Reconstructor> xmof_ir;
    std::vector<int> icells_to_reconstruct;
  };
}

#endif /* xmof2D_h */
