/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#ifndef ADVANCE_MOMENTS_H_
#define ADVANCE_MOMENTS_H_

#include <stdlib.h>
#include "tangram/support/tangram.h"
#include "tangram/utility/rk4.h"

/*!
 @brief For the given mesh, velocity field, and material data for cells' inverse images
 from the previous time step, get material moments on the current time step.
 @tparam Mesh_Wrapper A lightweight wrapper to a specific input mesh
 implementation that provides required functionality.

 @param[in] mesh Mesh wrapper.
 @param[in] mesh_mat_IDs IDs of materials in the mesh
 @param[in] mesh_mat_vols Total volumes of materials in the mesh
 @param[in] velocity Velocity field
 @param[in] t_prev Time of the previous step
 @param[in] dt Time step
 @param[in] cell_num_mats Number of material in each mesh cell, vector of length cell_num
 @param[in] cell_mat_ids Indices of materials in each mesh cell, a flat vector, requires
 computations of offsets
 @param[in, out] cell_mat_volfracs Volume fractions of materials in each mesh cell, a flat
 vector, requires computations of offsets (in: previous time step, out: current time step)
 @param[in, out] cell_mat_centroids Centroids of materials in each mesh cell, a flat vector,
 requires computations of offsets (in: previous time step, out: current time step)
 @param[in] vol_tol Volume tolerance                           
*/
template <class Mesh_Wrapper, int D>
void advance_moments(const Mesh_Wrapper& mesh,
                     const std::vector<int>& mesh_mat_IDs,
                     const std::vector<double>& mesh_mat_vols,
                     const std::function<Tangram::Vector<D>(double, 
                       const Tangram::Point<D>&)>& velocity,
                     double t_prev,
                     double dt,
                     const std::vector<int>& cell_num_mats,
                     const std::vector<int>& cell_mat_ids,
                     std::vector<double>& cell_mat_volfracs,
                     std::vector< Tangram::Point<D> >& cell_mat_centroids,
                     double vol_tol) {

  int ncells = mesh.num_owned_cells() + mesh.num_ghost_cells();
  std::vector<int> offsets(ncells, 0);
  for (int icell = 0; icell < ncells - 1; icell++)
    offsets[icell + 1] = offsets[icell] + cell_num_mats[icell];

  // Advance the positions of centroids
  for (int icell = 0; icell < ncells; icell++)
    for (int icmat = 0; icmat < cell_num_mats[icell]; icmat++)
      cell_mat_centroids[offsets[icell] + icmat] = Tangram::runge_kutta_4(
        velocity, cell_mat_centroids[offsets[icell] + icmat], t_prev, dt);

  // Calculate total volume fraction of material each cell was assigned
  std::vector<double> total_mat_volfracs(ncells, 0.0);
  for (int icell = 0; icell < ncells; icell++)
    for (int icmat = 0; icmat < cell_num_mats[icell]; icmat++)
      total_mat_volfracs[icell] += cell_mat_volfracs[offsets[icell] + icmat];

  // Determine cells that have too much material and push it out
  for (int icell = 0; icell < ncells; icell++) {
    double cell_size = mesh.cell_volume(icell);
    if (total_mat_volfracs[icell] < 1.0 + vol_tol/cell_size)
      continue;

    // Sort material indices in the order of increasing volume fractions
    int ncmats = cell_num_mats[icell];
    std::vector<int> sorted_icmats = Tangram::sorted_indices(std::vector<double>(
      cell_mat_volfracs.begin() + offsets[icell], 
      cell_mat_volfracs.begin() + offsets[icell] + ncmats), true);

    std::vector<int> iadj_cells;
    mesh.cell_get_node_adj_cells(icell, Tangram::Entity_type::PARALLEL_OWNED, &iadj_cells);
    for (int iscm = 0; iscm < ncmats; iscm++) {
      int icmat = sorted_icmats[iscm];
      int mat_id = cell_mat_ids[offsets[icell] + icmat];
      double extra_mat_vol = cell_size*cell_mat_volfracs[offsets[icell] + icmat]*
        (1.0 - 1.0/total_mat_volfracs[icell]);

      std::vector<int> ireceivers = iadj_cells;
      while (extra_mat_vol >= vol_tol/ncmats) {
        std::vector<int> ireceivers_upd;
        double slice_vol = extra_mat_vol/ireceivers.size();
        for (int irc = 0; irc < ireceivers.size(); irc++) {
          int iadj_cell = ireceivers[irc];
          double adj_cell_size = mesh.cell_volume(iadj_cell);
          double adjc_capacity = adj_cell_size*(1.0 - total_mat_volfracs[iadj_cell]);
          if (adjc_capacity < vol_tol)
            continue;
          
          int adjc_matid = std::distance(cell_mat_ids.begin() + offsets[iadj_cell],
            std::find(cell_mat_ids.begin() + offsets[iadj_cell],
                      cell_mat_ids.begin() + offsets[iadj_cell] + cell_num_mats[iadj_cell],
                      mat_id));
          if (adjc_matid == cell_num_mats[iadj_cell])
            continue;

          double vol_delta = (adjc_capacity < slice_vol) ? adjc_capacity : slice_vol;
          adjc_capacity -= vol_delta;

          double adjc_vf_delta = vol_delta/adj_cell_size;
          cell_mat_volfracs[offsets[iadj_cell] + adjc_matid] += adjc_vf_delta;
          total_mat_volfracs[iadj_cell] += adjc_vf_delta;
          if (adjc_capacity >= vol_tol)
            ireceivers_upd.push_back(iadj_cell);

          extra_mat_vol -= vol_delta;
          double vf_delta = vol_delta/cell_size;
          cell_mat_volfracs[offsets[icell] + icmat] -= vf_delta;
          total_mat_volfracs[icell] -= vf_delta;
        }

        if (!ireceivers_upd.empty())
          ireceivers = ireceivers_upd;
        else {
          double extra_vf_delta = extra_mat_vol/cell_size;
          cell_mat_volfracs[offsets[icell] + icmat] -= extra_vf_delta;
          total_mat_volfracs[icell] -= extra_vf_delta;
          break;
        }
      }
    }
  }

  // Determine how much material is left to be distributed
  std::vector<double> extra_mesh_mat_vols = mesh_mat_vols;
  for (int icell = 0; icell < ncells; icell++) {
    double cell_size = mesh.cell_volume(icell);
    for (int icmat = 0; icmat < cell_num_mats[icell]; icmat++) {
      int imatid = std::distance(mesh_mat_IDs.begin(),
        std::find(mesh_mat_IDs.begin(), mesh_mat_IDs.end(), 
                  cell_mat_ids[offsets[icell] + icmat]));      
      extra_mesh_mat_vols[imatid] -= cell_mat_volfracs[offsets[icell] + icmat]*cell_size;
    }
  }   

  // Determine cells that have too little material and fill them up
  for (int icell = 0; icell < ncells; icell++) {
    double cell_size = mesh.cell_volume(icell);
    if (total_mat_volfracs[icell] > 1.0 - vol_tol/cell_size)
      continue;  

    int ncmats = cell_num_mats[icell];

    // Find global leftover volumes for contained materials
    std::vector<int> imatids(ncmats);
    double total_extra_mat_vol = 0.0;
    for (int icmat = 0; icmat < ncmats; icmat++) {
      imatids[icmat] = std::distance(mesh_mat_IDs.begin(),
        std::find(mesh_mat_IDs.begin(), mesh_mat_IDs.end(), 
                  cell_mat_ids[offsets[icell] + icmat]));
      if (extra_mesh_mat_vols[imatids[icmat]] > 0.0)
        total_extra_mat_vol += extra_mesh_mat_vols[imatids[icmat]];
    }

    // Weight the material based on the leftover volumes
    std::vector<double> cmat_weights(ncmats);
    if (total_extra_mat_vol < vol_tol)
      cmat_weights.assign(ncmats, 1.0/ncmats);
    else {
      for (int icmat = 0; icmat < ncmats; icmat++)
        cmat_weights[icmat] = extra_mesh_mat_vols[imatids[icmat]]/total_extra_mat_vol;
    }

    //Fill the cell up to full size
    double empty_vf = 1.0 - total_mat_volfracs[icell];
    for (int icmat = 0; icmat < ncmats; icmat++) {
      if (cmat_weights[icmat] == 0.0)
        continue;

      double vf_delta = empty_vf*cmat_weights[icmat];
      cell_mat_volfracs[offsets[icell] + icmat] += vf_delta;
      total_mat_volfracs[icell] += vf_delta;
      extra_mesh_mat_vols[imatids[icmat]] -= vf_delta*cell_size;
    }
  }   

  for (int icell = 0; icell < ncells; icell++) {
    double mat_volfrac = 0.0;
    for (int icmat = 0; icmat < cell_num_mats[icell]; icmat++)
      mat_volfrac += cell_mat_volfracs[offsets[icell] + icmat];

    assert(std::fabs(mat_volfrac - 1.0) < vol_tol/mesh.cell_volume(icell));
  }
}


#endif