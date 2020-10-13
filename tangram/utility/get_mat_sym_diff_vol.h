/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#ifndef GET_MAT_SYM_DIFF_VOL_H_
#define GET_MAT_SYM_DIFF_VOL_H_

#include <stdlib.h>
#include "tangram/support/tangram.h"
#include "tangram/driver/CellMatPoly.h"
#include "tangram/intersect/split_rnd.h"

/*!
 @brief For a given cell, computes per-material volume of symmetric difference
 between reference polyhedra with that material and reconstructed MatPoly's from the
 CellMatPoly object.
 @param[in] reference_mat_polys for every cell material, the collection of reference
 single-material polyhedra containing that material
 @param[in] ref_mat_ids IDs of cell materials
 @param[in] ref_mat_vol reference volumes of cell materials
 @param[in] result_ptr pointer to CellMatPoly object containing the results of interface
 reconstruction for this cell
 @param[out] mat_sym_diff_vol per-material volumes of symmetric difference 
 @param[in] convex_matpolys If material polyhedra in CellMatPoly are non-convex, 
 this flag should be set to true in order to decompose them into tetrahedrons
*/
void get_mat_sym_diff_vol(const std::vector< std::vector<r3d_poly> >& reference_mat_polys,
                          const std::vector<int>& ref_mat_ids,
                          const std::vector<double>& ref_mat_vol,
                          const std::shared_ptr< Tangram::CellMatPoly<3> >& result_ptr,
                          std::vector<double>& mat_sym_diff_vol,
                          bool convex_matpolys) {
  int nref_mats = ref_mat_ids.size();
  mat_sym_diff_vol.resize(nref_mats);

  for (int icmat = 0; icmat < nref_mats; icmat++) {
    int material_id = ref_mat_ids[icmat];
    auto res_mat_polys = result_ptr->get_matpolys(material_id);

    if (res_mat_polys.empty())
      mat_sym_diff_vol[icmat] = ref_mat_vol[icmat];
    else {
      mat_sym_diff_vol[icmat] =
        result_ptr->material_moments(material_id)[0] + ref_mat_vol[icmat];

      for (auto&& ref_poly : reference_mat_polys[icmat]) {
        for (auto&& res_mat_poly : res_mat_polys) {
          std::vector<double> intersection_moments;
          Tangram::get_intersection_moments(res_mat_poly, ref_poly,
                                            intersection_moments, convex_matpolys);
          mat_sym_diff_vol[icmat] -= 2*intersection_moments[0];                               
        }
      }
      if (mat_sym_diff_vol[icmat] < 0.0)
        mat_sym_diff_vol[icmat] = 0.0;
    }
  }
}

/*!
 @brief For a given cell, computes per-material area of symmetric difference
 between reference polygons with that material and reconstructed MatPoly's from the
 CellMatPoly object.
 @param[in] reference_mat_polys for every cell material, the collection of reference
 single-material polygons containing that material
 @param[in] ref_mat_ids IDs of cell materials
 @param[in] ref_mat_vol reference areas of cell materials
 @param[in] result_ptr pointer to CellMatPoly object containing the results of interface
 reconstruction for this cell
 @param[out] mat_sym_diff_vol per-material areas of symmetric difference 
 @param[in] convex_matpolys If material polygons in CellMatPoly are non-convex, 
 this flag should be set to true in order to decompose them into triangles
*/
void get_mat_sym_diff_vol(const std::vector< std::vector<r2d_poly> >& reference_mat_polys,
                          const std::vector<int>& ref_mat_ids,
                          const std::vector<double>& ref_mat_vol,
                          const std::shared_ptr< Tangram::CellMatPoly<2> >& result_ptr,
                          std::vector<double>& mat_sym_diff_vol,
                          bool convex_matpolys) {
  int nref_mats = ref_mat_ids.size();
  mat_sym_diff_vol.resize(nref_mats);

  for (int icmat = 0; icmat < nref_mats; icmat++) {
    int material_id = ref_mat_ids[icmat];
    auto res_mat_polys = result_ptr->get_matpolys(material_id);

    if (res_mat_polys.empty())
      mat_sym_diff_vol[icmat] = ref_mat_vol[icmat];
    else {
      mat_sym_diff_vol[icmat] =
        result_ptr->material_moments(material_id)[0] + ref_mat_vol[icmat];

      for (auto&& ref_poly : reference_mat_polys[icmat]) {
        for (auto&& res_mat_poly : res_mat_polys) {
          std::vector<double> intersection_moments;
          Tangram::get_intersection_moments(res_mat_poly, ref_poly,
                                            intersection_moments, convex_matpolys);
          mat_sym_diff_vol[icmat] -= 2*intersection_moments[0];                               
        }
      }
      if (mat_sym_diff_vol[icmat] < 0.0)
        mat_sym_diff_vol[icmat] = 0.0;
    }
  }
}
#endif
