/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#ifndef GET_MATERIAL_MOMENTS_H_
#define GET_MATERIAL_MOMENTS_H_

#include <stdlib.h>
#include "tangram/support/tangram.h"
#include "tangram/intersect/split_r2d.h"
#include "tangram/intersect/split_r3d.h"
#include "tangram/utility/rpgtools/primitives.h"
#include "tangram/utility/rpgtools/cuts.h"

/*!
 @brief For a given mesh and a sequence of planes computes volume fractions
 and centroids. r3d is used to compute intersection moments. For every plane
 in the sequence, the remaining parts of the domain in the lower half-space 
 are assigned the corresponding material and excluded from further consideration.
 @tparam Mesh_Wrapper A lightweight wrapper to a specific input mesh
                      implementation that provides required functionality.

 @param[in] mesh Mesh wrapper.
 @param[in] planar_interfaces Sequence of planar material interfaces
 @param[in] material_IDs IDs of materials corresponding to the lower half-space of each
 plane; the last ID is for the material in the upper half-space of the last plane in the sequence
 @param[out] cell_num_mats Number of material in each mesh cell, vector of length cell_num
 @param[out] cell_mat_ids Indices of materials in each mesh cell, a flat vector, requires
                          computations of offsets
 @param[out] cell_mat_volfracs Volume fractions of materials in each mesh cell, a flat
                               vector, requires computations of offsets
 @param[out] cell_mat_centroids Centroids of materials in each mesh cell, a flat vector,
                                requires computations of offsets
 @param[out] reference_mat_polys For every cell and every material inside that cell, 
 the collection of single-material polyhedra containing that material
 @param[in] vol_tol Volume tolerance 
 @param[in] dst_tol Distance tolerance 
 @param[in] decompose_cells If mesh has non-convex cells, this flag should be set to true
 in order to decompose cells into tetrahedrons                            
*/
template <class Mesh_Wrapper>
void get_material_moments(const Mesh_Wrapper& mesh,
                          const std::vector< Tangram::Plane_t<3> >& planar_interfaces,
                          const std::vector<int>& material_IDs,
                          std::vector<int>& cell_num_mats,
                          std::vector<int>& cell_mat_ids,
                          std::vector<double>& cell_mat_volfracs,
                          std::vector< Tangram::Point<3> >& cell_mat_centroids,
                          const double vol_tol,
                          const double dst_tol,
                          const bool decompose_cells,
                          std::vector< std::vector< std::vector<r3d_poly> > >*
                            reference_mat_polys = nullptr) {
  int nplanes = static_cast<int>(planar_interfaces.size());
  assert(material_IDs.size() == nplanes + 1);
  
  int ncells = mesh.num_owned_cells() + mesh.num_ghost_cells();
  
  //Get convex MatPoly's for every cell
  std::vector< std::vector< Tangram::MatPoly<3> > > cells_polys(ncells);
  std::vector<int> cur_polys_cellID;
  for (int icell = 0; icell < ncells; icell++) {
    Tangram::MatPoly<3> mat_poly;
    Tangram::cell_get_matpoly(mesh, icell, &mat_poly, dst_tol);

    if (decompose_cells) {
      mat_poly.facetize_decompose(cells_polys[icell]);
      cur_polys_cellID.resize(
        cur_polys_cellID.size() + cells_polys[icell].size(), icell);
    }
    else {
      cells_polys[icell].push_back(mat_poly);
      cur_polys_cellID.push_back(icell);
    }
  }

  int ncur_polys = static_cast<int>(cur_polys_cellID.size());
  r3d_poly* cur_r3d_polys = new r3d_poly [ncur_polys];

  const int POLY_ORDER = 1;
  int nmoments = R3D_NUM_MOMENTS(POLY_ORDER);
  r3d_real r3d_moments[R3D_NUM_MOMENTS(POLY_ORDER)];

  int ir3d_poly = 0;
  std::vector< std::vector<double> > cur_polys_moments(ncur_polys);
  for (int icell = 0; icell < ncells; icell++)
    for (int ipoly = 0; ipoly < cells_polys[icell].size(); ipoly++) {
      Tangram::matpoly_to_r3dpoly(cells_polys[icell][ipoly], cur_r3d_polys[ir3d_poly]);

      r3d_reduce(&cur_r3d_polys[ir3d_poly], r3d_moments, POLY_ORDER);
      cur_polys_moments[ir3d_poly].resize(nmoments);
      for (int im = 0; im < nmoments; im++)
        cur_polys_moments[ir3d_poly][im] += r3d_moments[im];

      ir3d_poly++;
    }

  std::vector< std::vector<int> > cells_mat_ids(ncells);
  if (reference_mat_polys != nullptr) {
    reference_mat_polys->clear();
    reference_mat_polys->resize(ncells);
  }
  std::vector< std::vector< std::vector<double> > > cells_mat_moments(ncells);
  for (int iplane = 0; iplane < nplanes; iplane++) {
    int cur_mat_id = material_IDs[iplane];
    r3d_plane cur_r3d_plane;
    for (int ixyz = 0; ixyz < 3; ixyz++)
      cur_r3d_plane.n.xyz[ixyz] = planar_interfaces[iplane].normal[ixyz];
    cur_r3d_plane.d = planar_interfaces[iplane].dist2origin;

    r3d_poly* lower_hs_polys = new r3d_poly [ncur_polys];
    r3d_poly* upper_hs_polys = new r3d_poly [ncur_polys];
    r3d_split(cur_r3d_polys, (r3d_int) ncur_polys, cur_r3d_plane, 
              upper_hs_polys, lower_hs_polys);

    delete [] cur_r3d_polys;          

    std::vector<int> iremaining_polys;
    std::vector< std::vector<double> > remaining_polys_moments;
    for (int ipoly = 0; ipoly < ncur_polys; ipoly++) {
      bool nnz_cutoff = false;
      if (lower_hs_polys[ipoly].nverts > 0) {
        r3d_reduce(&lower_hs_polys[ipoly], r3d_moments, POLY_ORDER);
        if (r3d_moments[0] >= vol_tol) {
          nnz_cutoff = true;
          // Poly below the plane is cut off by the plane, add it to 
          // the cell's list of single-material poly's
          int icell = cur_polys_cellID[ipoly];
          int cell_mat_id = std::distance(cells_mat_ids[icell].begin(), 
            std::find(cells_mat_ids[icell].begin(), 
                      cells_mat_ids[icell].end(), cur_mat_id));

          if (cell_mat_id == cells_mat_ids[icell].size()) {
            cells_mat_ids[icell].resize(cell_mat_id + 1);
            cells_mat_moments[icell].resize(cell_mat_id + 1);
            cells_mat_moments[icell][cell_mat_id].resize(nmoments, 0.0);
            cells_mat_ids[icell][cell_mat_id] = cur_mat_id;
            if (reference_mat_polys != nullptr) 
              (*reference_mat_polys)[icell].resize(cell_mat_id + 1);
          }

          if (reference_mat_polys != nullptr)
            (*reference_mat_polys)[icell][cell_mat_id].push_back(lower_hs_polys[ipoly]);
          for (int im = 0; im < nmoments; im++)
            cells_mat_moments[icell][cell_mat_id][im] += r3d_moments[im];
        }
      }

      if (upper_hs_polys[ipoly].nverts > 0) {
        if (nnz_cutoff) {
          for (int im = 0; im < nmoments; im++)
            r3d_moments[im] = cur_polys_moments[ipoly][im] - r3d_moments[im];
        } 
        else {
          for (int im = 0; im < nmoments; im++)
            r3d_moments[im] = cur_polys_moments[ipoly][im];
        }
        if (r3d_moments[0] >= vol_tol) {
          int irpoly = iremaining_polys.size();
          iremaining_polys.push_back(ipoly);
          remaining_polys_moments.resize(irpoly + 1);
          remaining_polys_moments[irpoly].resize(nmoments);
          for (int im = 0; im < nmoments; im++)
            remaining_polys_moments[irpoly][im] = r3d_moments[im];
        }
      }
    }

    delete [] lower_hs_polys;

    ncur_polys = static_cast<int>(iremaining_polys.size());
    r3d_poly* remaining_polys = new r3d_poly [ncur_polys];
    std::vector<int> remaining_polys_cellID(ncur_polys);
    for (int irpoly = 0; irpoly < ncur_polys; irpoly++) {
      remaining_polys[irpoly] = upper_hs_polys[iremaining_polys[irpoly]];
      remaining_polys_cellID[irpoly] = cur_polys_cellID[iremaining_polys[irpoly]];
    }

    delete [] upper_hs_polys;
    cur_r3d_polys = remaining_polys;
    cur_polys_cellID = remaining_polys_cellID;
    cur_polys_moments = remaining_polys_moments;
  }

  //Poly's that are left after all planes were processed are also single-material
  int cur_mat_id = material_IDs[nplanes];
  for (int ipoly = 0; ipoly < ncur_polys; ipoly++) {
    int icell = cur_polys_cellID[ipoly];
    int cell_mat_id = std::distance(cells_mat_ids[icell].begin(),
      std::find(cells_mat_ids[icell].begin(), 
                cells_mat_ids[icell].end(), cur_mat_id));

    if (cell_mat_id == cells_mat_ids[icell].size()) {
      cells_mat_ids[icell].resize(cell_mat_id + 1);
      cells_mat_moments[icell].resize(cell_mat_id + 1);
      cells_mat_moments[icell][cell_mat_id].resize(nmoments, 0.0);
      cells_mat_ids[icell][cell_mat_id] = cur_mat_id;
      if (reference_mat_polys != nullptr)
        (*reference_mat_polys)[icell].resize(cell_mat_id + 1);
    }

    if (reference_mat_polys != nullptr)
      (*reference_mat_polys)[icell][cell_mat_id].push_back(cur_r3d_polys[ipoly]);
    for (int im = 0; im < nmoments; im++)
      cells_mat_moments[icell][cell_mat_id][im] += cur_polys_moments[ipoly][im];
  }

  delete [] cur_r3d_polys;

  cell_num_mats.resize(ncells);
  cell_mat_ids.clear();
  cell_mat_volfracs.clear();
  cell_mat_centroids.clear();

  for (int icell = 0; icell < ncells; icell++) {
    cell_num_mats[icell] = cells_mat_ids[icell].size();
    cell_mat_ids.insert(cell_mat_ids.end(), cells_mat_ids[icell].begin(), 
                        cells_mat_ids[icell].end());
    int offset = cell_mat_volfracs.size();
    cell_mat_volfracs.resize(offset + cell_num_mats[icell]);
    cell_mat_centroids.resize(offset + cell_num_mats[icell]);
    double cell_volume = mesh.cell_volume(icell);
    for (int imat = 0; imat < cell_num_mats[icell]; imat++) {
      cell_mat_volfracs[offset + imat] = 
        cells_mat_moments[icell][imat][0]/cell_volume;
      for (int idim = 0; idim < 3; idim++)
        cell_mat_centroids[offset + imat][idim] = 
          cells_mat_moments[icell][imat][idim + 1]/cells_mat_moments[icell][imat][0];
    }
  }
}

/*!
 @brief For a given mesh and a sphere computes volume fractions and centroids of 
 two material: one inside and the other outside of the sphere.
 r3d is used to compute intersection moments. 
 @tparam Mesh_Wrapper A lightweight wrapper to a specific input mesh
                      implementation that provides required functionality.

 @param[in] mesh Mesh wrapper.
 @param[in] material_IDs First entry is assigned to the material outside the sphere;
 the second entry is assigned to the material inside the sphere
 @param[in] center Center of the sphere
 @param[in] radius Radius of the sphere
 @param[in] nquadrant_samples Discrete representation of the sphere is used,
 this number determines how many points are sampled per a coordinate plane (xy, xz) 
 quadrant
 @param[out] cell_num_mats Number of material in each mesh cell, vector of length cell_num
 @param[out] cell_mat_ids Indices of materials in each mesh cell, a flat vector, requires
                          computations of offsets
 @param[out] cell_mat_volfracs Volume fractions of materials in each mesh cell, a flat
                               vector, requires computations of offsets
 @param[out] cell_mat_centroids Centroids of materials in each mesh cell, a flat vector,
                                requires computations of offsets
 @param[in] vol_tol Volume tolerance
 @param[in] dst_tol Distance tolerance
 @param[in] decompose_cells If mesh has non-convex cells, this flag should be set to true
 in order to decompose cells into tetrahedrons 
 @param[out] reference_mat_polys For every mesh cell and every material inside that cell, 
 a pointer to the collection of single-material polyhedra containing that material
                           
*/
template <class Mesh_Wrapper>
void get_material_moments(const Mesh_Wrapper& mesh,
                          const std::vector<int>& material_IDs,
                          const Tangram::Point<3>& center,
                          const double radius,
                          const int nquadrant_samples,
                          std::vector<int>& cell_num_mats,
                          std::vector<int>& cell_mat_ids,
                          std::vector<double>& cell_mat_volfracs,
                          std::vector< Tangram::Point<3> >& cell_mat_centroids,
                          const double vol_tol,
                          const double dst_tol,
                          const bool decompose_cells,
                          std::vector< std::vector< std::vector<r3d_poly> > >*
                            reference_mat_polys = nullptr) {
  assert(material_IDs.size() == 2);

  int ncells = mesh.num_owned_cells() + mesh.num_ghost_cells();
  
  Tangram::MatPoly<3> sphere_poly = sphere(center, radius, nquadrant_samples, dst_tol);

  std::vector< std::shared_ptr<RefPolyData_t> > mesh_polys;
  mesh_to_r3d_polys<Mesh_Wrapper>(mesh, mesh_polys, dst_tol, decompose_cells);

  std::vector< std::vector< std::shared_ptr<RefPolyData_t> > > ref_poly_sets(2);
  apply_poly(mesh_polys, sphere_poly, ref_poly_sets[1], ref_poly_sets[0],
             vol_tol, dst_tol, true);

  finalize_ref_data(mesh_polys, ref_poly_sets, material_IDs, cell_num_mats, cell_mat_ids, 
    cell_mat_volfracs, cell_mat_centroids, reference_mat_polys);
}

/*!
 @brief For a given mesh and concentric spheres computes volume fractions and centroids 
 of materials inside and outside of spheres. For each sphere in a sequence, the corresponding
 material is assigned to the volume inside the sphere that is not contained in any of the
 spheres preceding it in the sequence.
 r3d is used to compute intersection moments. 
 @tparam Mesh_Wrapper A lightweight wrapper to a specific input mesh
                      implementation that provides required functionality.

 @param[in] mesh Mesh wrapper.
 @param[in] material_IDs IDs of materials corresponding to the the interior of each
 sphere; the last ID is for the material outside of all the spheres
 @param[in] center Center of the spheres
 @param[in] radius Radius' of the spheres
 @param[in] nquadrant_samples Discrete representation of the spheres is used,
 this number determines how many points are sampled per a coordinate plane (xy, xz) 
 quadrant
 @param[out] cell_num_mats Number of material in each mesh cell, vector of length cell_num
 @param[out] cell_mat_ids Indices of materials in each mesh cell, a flat vector, requires
                          computations of offsets
 @param[out] cell_mat_volfracs Volume fractions of materials in each mesh cell, a flat
                               vector, requires computations of offsets
 @param[out] cell_mat_centroids Centroids of materials in each mesh cell, a flat vector,
                                requires computations of offsets
 @param[in] vol_tol Volume tolerance
 @param[in] dst_tol Distance tolerance 
 @param[in] decompose_cells If mesh has non-convex cells, this flag should be set to true
 in order to decompose cells into tetrahedrons 
 @param[out] reference_mat_polys For every mesh cell and every material inside that cell, 
 a pointer to the collection of single-material polyhedra containing that material
                           
*/
template <class Mesh_Wrapper>
void get_material_moments(const Mesh_Wrapper& mesh,
                          const std::vector<int>& material_IDs,
                          const Tangram::Point<3>& center,
                          const std::vector<double>& radius,
                          const int nquadrant_samples,
                          std::vector<int>& cell_num_mats,
                          std::vector<int>& cell_mat_ids,
                          std::vector<double>& cell_mat_volfracs,
                          std::vector< Tangram::Point<3> >& cell_mat_centroids,
                          const double vol_tol,
                          const double dst_tol,
                          const bool decompose_cells,
                          std::vector< std::vector< std::vector<r3d_poly> > >*
                            reference_mat_polys = nullptr) {
  int nspheres = static_cast<int>(radius.size());
  assert(material_IDs.size() == nspheres + 1);

  int ncells = mesh.num_owned_cells() + mesh.num_ghost_cells();
  
  std::vector< std::shared_ptr<RefPolyData_t> > mesh_polys, cur_polys_data, rem_polys_data;
  mesh_to_r3d_polys<Mesh_Wrapper>(mesh, mesh_polys, dst_tol, decompose_cells);

  std::vector< std::vector< std::shared_ptr<RefPolyData_t> > > ref_poly_sets(nspheres + 1);
  cur_polys_data = mesh_polys;
  for (int isphere = 0; isphere < nspheres; isphere++) {
    Tangram::MatPoly<3> sphere_poly = sphere(center, radius[isphere], 
                                             nquadrant_samples, dst_tol);
    apply_poly(cur_polys_data, sphere_poly, ref_poly_sets[isphere], rem_polys_data, 
               vol_tol, dst_tol, true);
    cur_polys_data = rem_polys_data;
    rem_polys_data.clear();
  }
  ref_poly_sets[nspheres] = cur_polys_data;
  finalize_ref_data(mesh_polys,ref_poly_sets, material_IDs, cell_num_mats, cell_mat_ids, 
    cell_mat_volfracs, cell_mat_centroids, reference_mat_polys);
}

/*!
 @brief For a given mesh and a sequence of lines computes volume fractions
 and centroids. r2d is used to compute intersection moments. For every line
 in the sequence, the remaining parts of the domain in the lower half-plane 
 are assigned the corresponding material and excluded from further consideration.
 @tparam Mesh_Wrapper A lightweight wrapper to a specific input mesh
                      implementation that provides required functionality.

 @param[in] mesh Mesh wrapper.
 @param[in] linear_interfaces Sequence of linear material interfaces
 @param[in] material_IDs IDs of materials corresponding to the lower half-plane of each
 line; the last ID is for the material in the upper half-plane of the last line in the sequence
 @param[out] cell_num_mats Number of material in each mesh cell, vector of length cell_num
 @param[out] cell_mat_ids Indices of materials in each mesh cell, a flat vector, requires
                          computations of offsets
 @param[out] cell_mat_volfracs Volume fractions of materials in each mesh cell, a flat
                               vector, requires computations of offsets
 @param[out] cell_mat_centroids Centroids of materials in each mesh cell, a flat vector,
                                requires computations of offsets
 @param[out] reference_mat_polys For every cell and every material inside that cell, 
 the collection of single-material polygons containing that material
 @param[in] vol_tol Volume tolerance
 @param[in] dst_tol Distance tolerance
 @param[in] decompose_cells If mesh has non-convex cells, this flag should be set to true
 in order to decompose cells into triangles                            
*/
template <class Mesh_Wrapper>
void get_material_moments(const Mesh_Wrapper& mesh,
                          const std::vector< Tangram::Plane_t<2> >& linear_interfaces,
                          const std::vector<int>& material_IDs,
                          std::vector<int>& cell_num_mats,
                          std::vector<int>& cell_mat_ids,
                          std::vector<double>& cell_mat_volfracs,
                          std::vector< Tangram::Point<2> >& cell_mat_centroids,
                          const double vol_tol,
                          const double dst_tol,
                          const bool decompose_cells,
                          std::vector< std::vector< std::vector<r2d_poly> > >*
                            reference_mat_polys = nullptr) {
  int nlines = static_cast<int>(linear_interfaces.size());
  assert(material_IDs.size() == nlines + 1);
  
  int ncells = mesh.num_owned_cells() + mesh.num_ghost_cells();
  
  //Get convex MatPoly's for every cell
  std::vector< std::vector< Tangram::MatPoly<2> > > cells_polys(ncells);
  std::vector<int> cur_polys_cellID;
  for (int icell = 0; icell < ncells; icell++) {
    Tangram::MatPoly<2> mat_poly;
    Tangram::cell_get_matpoly(mesh, icell, &mat_poly, dst_tol);

    if (decompose_cells) {
      mat_poly.decompose(cells_polys[icell]);
      cur_polys_cellID.resize(
        cur_polys_cellID.size() + cells_polys[icell].size(), icell);
    }
    else {
      cells_polys[icell].push_back(mat_poly);
      cur_polys_cellID.push_back(icell);
    }
  }

  int ncur_polys = static_cast<int>(cur_polys_cellID.size());
  r2d_poly* cur_r2d_polys = new r2d_poly [ncur_polys];

  const int POLY_ORDER = 1;
  int nmoments = R2D_NUM_MOMENTS(POLY_ORDER);
  r2d_real r2d_moments[R2D_NUM_MOMENTS(POLY_ORDER)];

  int ir2d_poly = 0;
  std::vector< std::vector<double> > cur_polys_moments(ncur_polys);
  for (int icell = 0; icell < ncells; icell++)
    for (int ipoly = 0; ipoly < cells_polys[icell].size(); ipoly++) {
      Tangram::matpoly_to_r2dpoly(cells_polys[icell][ipoly], cur_r2d_polys[ir2d_poly]);

      r2d_reduce(&cur_r2d_polys[ir2d_poly], r2d_moments, POLY_ORDER);
      cur_polys_moments[ir2d_poly].resize(nmoments);
      for (int im = 0; im < nmoments; im++)
        cur_polys_moments[ir2d_poly][im] += r2d_moments[im];

      ir2d_poly++;
    }

  std::vector< std::vector<int> > cells_mat_ids(ncells);
  if (reference_mat_polys != nullptr) {
    reference_mat_polys->clear();
    reference_mat_polys->resize(ncells);
  }
  std::vector< std::vector< std::vector<double> > > cells_mat_moments(ncells);
  for (int iline = 0; iline < nlines; iline++) {
    int cur_mat_id = material_IDs[iline];
    r2d_plane cur_r2d_plane;
    for (int ixy = 0; ixy < 2; ixy++)
      cur_r2d_plane.n.xy[ixy] = -linear_interfaces[iline].normal[ixy];
    cur_r2d_plane.d = -linear_interfaces[iline].dist2origin;

    r2d_poly* lower_hp_polys = new r2d_poly [ncur_polys];
    r2d_poly* upper_hp_polys = new r2d_poly [ncur_polys];
    r2d_split(cur_r2d_polys, (r2d_int) ncur_polys, cur_r2d_plane, 
              upper_hp_polys, lower_hp_polys);

    delete [] cur_r2d_polys;          

    std::vector<int> iremaining_polys;
    std::vector< std::vector<double> > remaining_polys_moments;
    for (int ipoly = 0; ipoly < ncur_polys; ipoly++) {
      bool nnz_cutoff = false;
      if (lower_hp_polys[ipoly].nverts > 0) {
        r2d_reduce(&lower_hp_polys[ipoly], r2d_moments, POLY_ORDER);
        if (r2d_moments[0] >= vol_tol) {
          nnz_cutoff = true;
          // Poly below the line is cut off by the line, add it to 
          // the cell's list of single-material poly's
          int icell = cur_polys_cellID[ipoly];
          int cell_mat_id = std::distance(cells_mat_ids[icell].begin(),
            std::find(cells_mat_ids[icell].begin(), 
                      cells_mat_ids[icell].end(), cur_mat_id));

          if (cell_mat_id == cells_mat_ids[icell].size()) {
            cells_mat_ids[icell].resize(cell_mat_id + 1);
            cells_mat_moments[icell].resize(cell_mat_id + 1);
            cells_mat_moments[icell][cell_mat_id].resize(nmoments, 0.0);
            cells_mat_ids[icell][cell_mat_id] = cur_mat_id;
            if (reference_mat_polys != nullptr)
              (*reference_mat_polys)[icell].resize(cell_mat_id + 1);
          }

          if (reference_mat_polys != nullptr)
            (*reference_mat_polys)[icell][cell_mat_id].push_back(lower_hp_polys[ipoly]);
          for (int im = 0; im < nmoments; im++)
            cells_mat_moments[icell][cell_mat_id][im] += r2d_moments[im];
        }
      }

      if (upper_hp_polys[ipoly].nverts > 0) {
        if (nnz_cutoff) {
          for (int im = 0; im < nmoments; im++)
            r2d_moments[im] = cur_polys_moments[ipoly][im] - r2d_moments[im];
        } 
        else {
          for (int im = 0; im < nmoments; im++)
            r2d_moments[im] = cur_polys_moments[ipoly][im];
        }
        if (r2d_moments[0] >= vol_tol) {
          int irpoly = iremaining_polys.size();
          iremaining_polys.push_back(ipoly);
          remaining_polys_moments.resize(irpoly + 1);
          remaining_polys_moments[irpoly].resize(nmoments);
          for (int im = 0; im < nmoments; im++)
            remaining_polys_moments[irpoly][im] = r2d_moments[im];
        }
      }
    }

    delete [] lower_hp_polys;

    ncur_polys = static_cast<int>(iremaining_polys.size());
    r2d_poly* remaining_polys = new r2d_poly [ncur_polys];
    std::vector<int> remaining_polys_cellID(ncur_polys);
    for (int irpoly = 0; irpoly < ncur_polys; irpoly++) {
      remaining_polys[irpoly] = upper_hp_polys[iremaining_polys[irpoly]];
      remaining_polys_cellID[irpoly] = cur_polys_cellID[iremaining_polys[irpoly]];
    }

    delete [] upper_hp_polys;
    cur_r2d_polys = remaining_polys;
    cur_polys_cellID = remaining_polys_cellID;
    cur_polys_moments = remaining_polys_moments;
  }

  //Poly's that are left after all lines were processed are also single-material
  int cur_mat_id = material_IDs[nlines];
  for (int ipoly = 0; ipoly < ncur_polys; ipoly++) {
    int icell = cur_polys_cellID[ipoly];
    int cell_mat_id = std::distance(cells_mat_ids[icell].begin(),
      std::find(cells_mat_ids[icell].begin(), 
                cells_mat_ids[icell].end(), cur_mat_id));

    if (cell_mat_id == cells_mat_ids[icell].size()) {
      cells_mat_ids[icell].resize(cell_mat_id + 1);
      cells_mat_moments[icell].resize(cell_mat_id + 1);
      cells_mat_moments[icell][cell_mat_id].resize(nmoments, 0.0);
      cells_mat_ids[icell][cell_mat_id] = cur_mat_id;
      if (reference_mat_polys != nullptr)
        (*reference_mat_polys)[icell].resize(cell_mat_id + 1);
    }

    if (reference_mat_polys != nullptr)
      (*reference_mat_polys)[icell][cell_mat_id].push_back(cur_r2d_polys[ipoly]);
    for (int im = 0; im < nmoments; im++)
      cells_mat_moments[icell][cell_mat_id][im] += cur_polys_moments[ipoly][im];
  }

  delete [] cur_r2d_polys;

  cell_num_mats.resize(ncells);
  cell_mat_ids.clear();
  cell_mat_volfracs.clear();
  cell_mat_centroids.clear();

  for (int icell = 0; icell < ncells; icell++) {
    cell_num_mats[icell] = cells_mat_ids[icell].size();
    cell_mat_ids.insert(cell_mat_ids.end(), cells_mat_ids[icell].begin(), 
                        cells_mat_ids[icell].end());
    int offset = cell_mat_volfracs.size();
    cell_mat_volfracs.resize(offset + cell_num_mats[icell]);
    cell_mat_centroids.resize(offset + cell_num_mats[icell]);
    double cell_volume = mesh.cell_volume(icell);
    for (int imat = 0; imat < cell_num_mats[icell]; imat++) {
      cell_mat_volfracs[offset + imat] = 
        cells_mat_moments[icell][imat][0]/cell_volume;
      for (int idim = 0; idim < 2; idim++)
        cell_mat_centroids[offset + imat][idim] = 
          cells_mat_moments[icell][imat][idim + 1]/cells_mat_moments[icell][imat][0];
    }
  }
}

/*!
 @brief For a given mesh and a circle computes volume fractions and centroids of 
 two material: one inside and the other outside of the circle.
 r2d is used to compute intersection moments. 
 @tparam Mesh_Wrapper A lightweight wrapper to a specific input mesh
                      implementation that provides required functionality.

 @param[in] mesh Mesh wrapper.
 @param[in] material_IDs First entry is assigned to the material outside the circle;
 the second entry is assigned to the material inside the circle
 @param[in] center Center of the circle
 @param[in] radius Radius of the circle
 @param[in] nquadrant_samples Discrete representation of the circle is used,
 this number determines how many points are sampled per the xy-coordinate plane quadrant
 @param[out] cell_num_mats Number of material in each mesh cell, vector of length cell_num
 @param[out] cell_mat_ids Indices of materials in each mesh cell, a flat vector, requires
                          computations of offsets
 @param[out] cell_mat_volfracs Volume fractions of materials in each mesh cell, a flat
                               vector, requires computations of offsets
 @param[out] cell_mat_centroids Centroids of materials in each mesh cell, a flat vector,
                                requires computations of offsets
 @param[in] vol_tol Volume tolerance
 @param[in] dst_tol Distance tolerance 
 @param[in] decompose_cells If mesh has non-convex cells, this flag should be set to true
 in order to decompose cells into tetrahedrons 
 @param[out] reference_mat_polys For every mesh cell and every material inside that cell, 
 a pointer to the collection of single-material polyhedra containing that material
                           
*/
template <class Mesh_Wrapper>
void get_material_moments(const Mesh_Wrapper& mesh,
                          const std::vector<int>& material_IDs,
                          const Tangram::Point<2>& center,
                          const double radius,
                          const int nquadrant_samples,
                          std::vector<int>& cell_num_mats,
                          std::vector<int>& cell_mat_ids,
                          std::vector<double>& cell_mat_volfracs,
                          std::vector< Tangram::Point<2> >& cell_mat_centroids,
                          const double vol_tol,
                          const double dst_tol,
                          const bool decompose_cells,
                          std::vector< std::vector< std::vector<r2d_poly> > >*
                            reference_mat_polys = nullptr) {
  assert(material_IDs.size() == 2);

  int ncells = mesh.num_owned_cells() + mesh.num_ghost_cells();
  
  // Create a MatPoly that represents the circle
  int ncircle_samples = 4*nquadrant_samples;
  std::vector<Tangram::Point<2>> circle_pts;
  circle_pts.reserve(ncircle_samples);

  for (int ipt = 0; ipt < ncircle_samples; ipt++) {
    double angle = 2*ipt*PI/ncircle_samples;
    circle_pts.push_back(Tangram::Point<2>(
      center[0] + radius*cos(angle),
      center[1] + radius*sin(angle)));
  }

  Tangram::MatPoly<2> circle_poly;
  circle_poly.initialize(circle_pts, dst_tol);

  // Get the lines containing the circle's faces
  std::vector< Tangram::Plane_t<2> > face_lines;
  circle_poly.face_planes(face_lines);

  //Get convex MatPoly's for every cell
  std::vector< std::vector< Tangram::MatPoly<2> > > cells_polys(ncells);
  std::vector<int> cur_polys_cellID;
  for (int icell = 0; icell < ncells; icell++) {
    Tangram::MatPoly<2> mat_poly;
    Tangram::cell_get_matpoly(mesh, icell, &mat_poly);

    if (decompose_cells) {
      mat_poly.decompose(cells_polys[icell]);
      cur_polys_cellID.resize(
        cur_polys_cellID.size() + cells_polys[icell].size(), icell);
    }
    else {
      cells_polys[icell].push_back(mat_poly);
      cur_polys_cellID.push_back(icell);
    }
  }

  int ncur_polys = static_cast<int>(cur_polys_cellID.size());
  r2d_poly* cur_r2d_polys = new r2d_poly [ncur_polys];

  const int POLY_ORDER = 1;
  int nmoments = R2D_NUM_MOMENTS(POLY_ORDER);
  r2d_real r2d_moments[R2D_NUM_MOMENTS(POLY_ORDER)];

  int ir2d_poly = 0;
  std::vector< std::vector<double> > cur_polys_moments(ncur_polys);
  for (int icell = 0; icell < ncells; icell++)
    for (int ipoly = 0; ipoly < cells_polys[icell].size(); ipoly++) {
      Tangram::matpoly_to_r2dpoly(cells_polys[icell][ipoly], cur_r2d_polys[ir2d_poly]);

      r2d_reduce(&cur_r2d_polys[ir2d_poly], r2d_moments, POLY_ORDER);
      cur_polys_moments[ir2d_poly].resize(nmoments);
      for (int im = 0; im < nmoments; im++)
        cur_polys_moments[ir2d_poly][im] += r2d_moments[im];

      ir2d_poly++;
    }

  std::vector< std::vector<int> > cells_mat_ids(ncells);
  if (reference_mat_polys != nullptr) {
    reference_mat_polys->clear();
    reference_mat_polys->resize(ncells);
  }
  std::vector< std::vector< std::vector<double> > > cells_mat_moments(ncells);
  int nlines = static_cast<int>(face_lines.size());
  for (int iline = 0; iline < nlines; iline++) {
    r2d_plane cur_r2d_plane;
    for (int ixy = 0; ixy < 2; ixy++)
      cur_r2d_plane.n.xy[ixy] = -face_lines[iline].normal[ixy];
    cur_r2d_plane.d = -face_lines[iline].dist2origin;

    r2d_poly* lower_hp_polys = new r2d_poly [ncur_polys];
    r2d_poly* upper_hp_polys = new r2d_poly [ncur_polys];
    r2d_split(cur_r2d_polys, (r2d_int) ncur_polys, cur_r2d_plane, 
              upper_hp_polys, lower_hp_polys);

    delete [] cur_r2d_polys;          

    std::vector<int> iremaining_polys;
    std::vector< std::vector<double> > remaining_polys_moments;
    for (int ipoly = 0; ipoly < ncur_polys; ipoly++) {
      bool nnz_cutoff = false;
      if (upper_hp_polys[ipoly].nverts > 0) {
        r2d_reduce(&upper_hp_polys[ipoly], r2d_moments, POLY_ORDER);
        if (r2d_moments[0] >= vol_tol) {
          nnz_cutoff = true;
          // Poly above the line is exterior with respect to the circle,
          // add it to the list of cell's poly's, which currently can only
          // contain exterior material with ID material_IDs[0]
          int icell = cur_polys_cellID[ipoly];

          if (cells_mat_ids[icell].empty()) {
            cells_mat_ids[icell].push_back(material_IDs[0]);
            cells_mat_moments[icell].resize(1);
            cells_mat_moments[icell][0].resize(nmoments, 0.0);

            if (reference_mat_polys != nullptr)
              (*reference_mat_polys)[icell].resize(1);
          }

          if (reference_mat_polys != nullptr)
            (*reference_mat_polys)[icell][0].push_back(upper_hp_polys[ipoly]);

          for (int im = 0; im < nmoments; im++)
            cells_mat_moments[icell][0][im] += r2d_moments[im];
        }
      }

      if (lower_hp_polys[ipoly].nverts > 0) {
        if (nnz_cutoff) {
          for (int im = 0; im < nmoments; im++)
            r2d_moments[im] = cur_polys_moments[ipoly][im] - r2d_moments[im];
        } 
        else {
          for (int im = 0; im < nmoments; im++)
            r2d_moments[im] = cur_polys_moments[ipoly][im];
        }
        if (r2d_moments[0] >= vol_tol) {
          int irpoly = iremaining_polys.size();
          iremaining_polys.push_back(ipoly);
          remaining_polys_moments.resize(irpoly + 1);
          remaining_polys_moments[irpoly].resize(nmoments);
          for (int im = 0; im < nmoments; im++)
            remaining_polys_moments[irpoly][im] = r2d_moments[im];
        }
      }
    }

    delete [] upper_hp_polys;

    ncur_polys = static_cast<int>(iremaining_polys.size());
    r2d_poly* remaining_polys = new r2d_poly [ncur_polys];
    std::vector<int> remaining_polys_cellID(ncur_polys);
    for (int irpoly = 0; irpoly < ncur_polys; irpoly++) {
      remaining_polys[irpoly] = lower_hp_polys[iremaining_polys[irpoly]];
      remaining_polys_cellID[irpoly] = cur_polys_cellID[iremaining_polys[irpoly]];
    }

    delete [] lower_hp_polys;
    cur_r2d_polys = remaining_polys;
    cur_polys_cellID = remaining_polys_cellID;
    cur_polys_moments = remaining_polys_moments;
  }

  //Poly's that are left after all lines were processed make up the circle
  for (int ipoly = 0; ipoly < ncur_polys; ipoly++) {
    int icell = cur_polys_cellID[ipoly];
    int cell_mat_id = std::distance(cells_mat_ids[icell].begin(),
      std::find(cells_mat_ids[icell].begin(), 
                cells_mat_ids[icell].end(), material_IDs[1]));

    if (cell_mat_id == cells_mat_ids[icell].size()) {
      cells_mat_ids[icell].push_back(material_IDs[1]);
      cells_mat_moments[icell].resize(cell_mat_id + 1);
      cells_mat_moments[icell][cell_mat_id].resize(nmoments, 0.0);

      if (reference_mat_polys != nullptr)
        (*reference_mat_polys)[icell].resize(cell_mat_id + 1);
    }

    if (reference_mat_polys != nullptr)
      (*reference_mat_polys)[icell][cell_mat_id].push_back(cur_r2d_polys[ipoly]);

    for (int im = 0; im < nmoments; im++)
      cells_mat_moments[icell][cell_mat_id][im] += cur_polys_moments[ipoly][im];
  }

  delete [] cur_r2d_polys;

  cell_num_mats.resize(ncells);
  cell_mat_ids.clear();
  cell_mat_volfracs.clear();
  cell_mat_centroids.clear();

  for (int icell = 0; icell < ncells; icell++) {
    cell_num_mats[icell] = cells_mat_ids[icell].size();
    cell_mat_ids.insert(cell_mat_ids.end(), cells_mat_ids[icell].begin(), 
                        cells_mat_ids[icell].end());
    int offset = cell_mat_volfracs.size();
    cell_mat_volfracs.resize(offset + cell_num_mats[icell]);
    cell_mat_centroids.resize(offset + cell_num_mats[icell]);
    double cell_volume = mesh.cell_volume(icell);
    for (int imat = 0; imat < cell_num_mats[icell]; imat++) {
      cell_mat_volfracs[offset + imat] = 
        cells_mat_moments[icell][imat][0]/cell_volume;
      for (int idim = 0; idim < 2; idim++)
        cell_mat_centroids[offset + imat][idim] = 
          cells_mat_moments[icell][imat][idim + 1]/cells_mat_moments[icell][imat][0];
    }
  }
}

#endif
