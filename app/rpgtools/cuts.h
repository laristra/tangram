/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#ifndef RPG_TOOLS_CUTS_H_
#define RPG_TOOLS_CUTS_H_

#include <stdlib.h>
#include "tangram/support/tangram.h"
#include "tangram/intersect/split_r3d.h"

#define R3D_POLY_ORDER 1

struct RefPolySet_t {
  Tangram::vector<r3d_poly> r3d_polys;
  std::vector<int> ipoly2cellID;
  std::vector< std::vector<double> > moments;

  void clear() { r3d_polys.clear(); ipoly2cellID.clear(); moments.clear(); }
};

class r3d_poly_intersect_check {
public:
  r3d_poly_intersect_check(const Tangram::vector<r3d_poly>& r3dpolys,
                           const std::vector<int>& IDs_to_check,
                           const Tangram::MatPoly<3>& convex_matpoly) :
                           r3dpolys_(r3dpolys), IDs_to_check_(IDs_to_check),
                           convex_matpoly_(convex_matpoly) {
    std::vector< Tangram::Plane_t<3> > fplanes;
    convex_matpoly_.face_planes(fplanes);

    int nplanes = (int) fplanes.size();
    face_planes_.resize(nplanes);
    //Face planes, used to find the intersection with convex_matpoly
    for (int iplane = 0; iplane < nplanes; iplane++) {
      for (int ixyz = 0; ixyz < 3; ixyz++)
        face_planes_[iplane].n.xyz[ixyz] = -fplanes[iplane].normal[ixyz];
      face_planes_[iplane].d = -fplanes[iplane].dist2origin;
    }
  }

  //return value: -1 if outside, 1 if inside, 0 if intersects
  int operator()(const int checkID) const {
    int polyID = IDs_to_check_[checkID];
    bool intersects = false;
    for (int ivrt = 0; ivrt < r3dpolys_[polyID].nverts; ivrt++) {
      Tangram::Point3 cur_vrt;
      for (int ixyz = 0; ixyz < 3; ixyz++)
        cur_vrt[ixyz] = r3dpolys_[polyID].verts[ivrt].pos.xyz[ixyz];

      if (Tangram::point_inside_matpoly(convex_matpoly_, cur_vrt, true)) {
        if (!intersects && (ivrt > 0)) return 0;
        intersects = true;
      }
      else if (intersects) return 0;  
    }

    if (intersects) return 1;

    r3d_poly intersection = r3dpolys_[polyID];
    for (int iplane = 0; iplane < face_planes_.size(); iplane++) {
      r3d_clip(&intersection, &face_planes_[iplane], 1);
      if (intersection.nverts == 0) return -1;
    }
    
    return 0;
  }

private:
  const Tangram::vector<r3d_poly>& r3dpolys_;
  const std::vector<int>& IDs_to_check_;
  const Tangram::MatPoly<3>& convex_matpoly_;
  std::vector<r3d_plane> face_planes_;
};

template <int POLY_ORDER>
struct r3d_hs_polys_t {
  r3d_poly lower_hs_poly;
  r3d_poly upper_hs_poly;
  r3d_real lower_hs_moments[R3D_NUM_MOMENTS(POLY_ORDER)];
  r3d_real upper_hs_moments[R3D_NUM_MOMENTS(POLY_ORDER)];
};

template <int POLY_ORDER>
class r3d_split_operator {
public:
  r3d_split_operator(const Tangram::vector<r3d_poly>& r3dpolys,
                     const r3d_plane& cutting_plane,
                     const std::vector< std::vector<double> >& poly_moments) :
                     r3dpolys_(r3dpolys), cutting_plane_(cutting_plane),
                     poly_moments_(poly_moments) {}

  r3d_hs_polys_t<POLY_ORDER> operator()(const int polyID) const {
    r3d_hs_polys_t<POLY_ORDER> hs_polys;

    r3d_split(&r3dpolys_[polyID], 1, cutting_plane_, 
              &hs_polys.upper_hs_poly, &hs_polys.lower_hs_poly);

    bool nnz_cutoff = false;
    if (hs_polys.lower_hs_poly.nverts > 0) {
      r3d_reduce(&hs_polys.lower_hs_poly, hs_polys.lower_hs_moments, POLY_ORDER);
      if (hs_polys.lower_hs_moments[0] > std::numeric_limits<double>::epsilon())
        nnz_cutoff = true;
      else
        hs_polys.lower_hs_poly.nverts = 0;
    }

    if (hs_polys.upper_hs_poly.nverts > 0) {
      if (nnz_cutoff) {
        for (int im = 0; im < R3D_NUM_MOMENTS(POLY_ORDER); im++)
          hs_polys.upper_hs_moments[im] = 
            poly_moments_[polyID][im] - hs_polys.lower_hs_moments[im];
      } 
      else {
        for (int im = 0; im < R3D_NUM_MOMENTS(POLY_ORDER); im++)
          hs_polys.upper_hs_moments[im] = poly_moments_[polyID][im];
      }
      if (hs_polys.upper_hs_moments[0] <= std::numeric_limits<double>::epsilon())
        hs_polys.upper_hs_poly.nverts = 0;
    }

    return hs_polys;
  }

private:
  const Tangram::vector<r3d_poly>& r3dpolys_;
  const r3d_plane& cutting_plane_;
  const std::vector< std::vector<double> >& poly_moments_;
};

template <class Mesh_Wrapper>
void mesh_to_r3d_polys(const Mesh_Wrapper& mesh,
                       RefPolySet_t& all_polys,
                       bool decompose_cells = true) {
  int ncells = mesh.num_owned_cells() + mesh.num_ghost_cells();

  all_polys.r3d_polys.clear();
  all_polys.ipoly2cellID.clear();
  all_polys.moments.clear();

  int nmoments = R3D_NUM_MOMENTS(R3D_POLY_ORDER);
  r3d_real r3d_moments[R3D_NUM_MOMENTS(R3D_POLY_ORDER)];

  for (int icell = 0; icell < ncells; icell++) {
    Tangram::MatPoly<3> mat_poly;
    Tangram::cell_get_matpoly(mesh, icell, &mat_poly);

    std::vector< Tangram::MatPoly<3> > cell_polys;
    if (decompose_cells) {
      mat_poly.facetize_decompose(cell_polys);
      all_polys.ipoly2cellID.resize(
        all_polys.ipoly2cellID.size() + cell_polys.size(), icell);
    }
    else {
      cell_polys.push_back(mat_poly);
      all_polys.ipoly2cellID.push_back(icell);
    }

    int num_polys = (int) cell_polys.size();
    int offset = (int) all_polys.r3d_polys.size();
    all_polys.r3d_polys.resize(offset + num_polys);
    all_polys.moments.resize(offset + num_polys);

    for (int icp = 0; icp < num_polys; icp++) {
      Tangram::matpoly_to_r3dpoly(cell_polys[icp], 
                                  all_polys.r3d_polys[offset + icp]);

      r3d_reduce(&all_polys.r3d_polys[offset + icp], r3d_moments, R3D_POLY_ORDER);
      all_polys.moments[offset + icp].resize(nmoments);
      for (int im = 0; im < nmoments; im++)
        all_polys.moments[offset + icp][im] += r3d_moments[im];
    }
  }
}

void apply_plane(const Tangram::Plane_t<3>& planar_interface,
                 RefPolySet_t& multimat_polys,
                 RefPolySet_t& singlemat_polys) {
  assert(!multimat_polys.r3d_polys.empty());

  singlemat_polys.r3d_polys.clear();
  singlemat_polys.ipoly2cellID.clear();
  singlemat_polys.moments.clear();
  
  int num_mmpolys = (int) multimat_polys.r3d_polys.size();
  int nmoments = (int) multimat_polys.moments[0].size();

  r3d_plane r3d_planar_iface;
  for (int ixyz = 0; ixyz < 3; ixyz++)
    r3d_planar_iface.n.xyz[ixyz] = planar_interface.normal[ixyz];
  r3d_planar_iface.d = planar_interface.dist2origin;

  r3d_split_operator<R3D_POLY_ORDER> 
    r3d_split_op(multimat_polys.r3d_polys, r3d_planar_iface, multimat_polys.moments);
  Tangram::vector< r3d_hs_polys_t<R3D_POLY_ORDER> > hs_polys(num_mmpolys);

  Tangram::transform(Tangram::make_counting_iterator(0), 
                     Tangram::make_counting_iterator(num_mmpolys),
                     hs_polys.begin(), r3d_split_op);

  std::vector<int> ism_polys, imm_polys;
  for (int ipoly = 0; ipoly < num_mmpolys; ipoly++) {
    if (hs_polys[ipoly].lower_hs_poly.nverts > 0) {
      // Poly below the plane is added to the list of single-material poly's
      int ism_poly = ism_polys.size();
      ism_polys.push_back(ipoly);
    }

    if (hs_polys[ipoly].upper_hs_poly.nverts > 0) {
      // Poly above the plane is added to the list of multi-material poly's
      int imm_poly = imm_polys.size();
      imm_polys.push_back(ipoly);
    }
  }

  int num_smpolys = (int) ism_polys.size();
  singlemat_polys.r3d_polys.resize(num_smpolys);
  singlemat_polys.ipoly2cellID.resize(num_smpolys);
  singlemat_polys.moments.resize(num_smpolys);
  for (int ismp = 0; ismp < num_smpolys; ismp++) {
    int isplit_poly = ism_polys[ismp];
    singlemat_polys.r3d_polys[ismp] = hs_polys[isplit_poly].lower_hs_poly;
    singlemat_polys.ipoly2cellID[ismp] = 
      multimat_polys.ipoly2cellID[isplit_poly];
    singlemat_polys.moments[ismp].resize(nmoments);
    for (int im = 0; im < nmoments; im++)
      singlemat_polys.moments[ismp][im] = 
        hs_polys[isplit_poly].lower_hs_moments[im];
  }
      
  num_mmpolys = (int) imm_polys.size();
  multimat_polys.r3d_polys.resize(num_mmpolys);
  std::vector<int> irem_poly2cellID(num_mmpolys);
  multimat_polys.moments.resize(num_mmpolys);
  for (int immp = 0; immp < num_mmpolys; immp++) {
    int isplit_poly = imm_polys[immp];
    multimat_polys.r3d_polys[immp] = hs_polys[isplit_poly].upper_hs_poly;
    irem_poly2cellID[immp] = 
      multimat_polys.ipoly2cellID[isplit_poly];
    multimat_polys.moments[immp].resize(nmoments);
    for (int im = 0; im < nmoments; im++)
      multimat_polys.moments[immp][im] = 
        hs_polys[isplit_poly].upper_hs_moments[im];
  }
  
  multimat_polys.ipoly2cellID = irem_poly2cellID;
}

void sort_wrt_convex_poly(const Tangram::MatPoly<3>& convex_matpoly,
                          const Tangram::vector<r3d_poly>& all_polys,
                          std::vector<int>& iexterior_polys,
                          std::vector<int>& iinstersecting_polys,
                          std::vector<int>& iinterior_polys ) {
  assert (!all_polys.empty());

  std::vector< Tangram::Plane_t<3> > face_planes;
  convex_matpoly.face_planes(face_planes);
  Tangram::BoundingBox_t<3> poly_bbox = convex_matpoly.bounding_box();

  int nfaces = (int) face_planes.size();
  int npolys = (int) all_polys.size(); 

  iexterior_polys.clear();
  std::vector<int> iin_box_polys;
  for (int ipoly = 0; ipoly < npolys; ipoly++) {
    Tangram::BoundingBox_t<3> cur_bbox = 
      Tangram::r3d_poly_bounding_box(all_polys[ipoly]);
    if(overlapping_boxes(cur_bbox, poly_bbox))
      iin_box_polys.push_back(ipoly);
    else
      iexterior_polys.push_back(ipoly);
  }

  //Sort out the polys inside the box
  int num_in_box = (int) iin_box_polys.size();

  r3d_poly_intersect_check 
    r3d_isect_check(all_polys, iin_box_polys, convex_matpoly);
  Tangram::vector<int> check_result(num_in_box);

  Tangram::transform(Tangram::make_counting_iterator(0), 
                     Tangram::make_counting_iterator(num_in_box),
                     check_result.begin(), r3d_isect_check);

  iinterior_polys.clear();
  iinstersecting_polys.clear();
  for (int iip = 0; iip < num_in_box; iip++)
    switch (check_result[iip]) {
      case -1: iexterior_polys.push_back(iin_box_polys[iip]);
               break;
      case  0: iinstersecting_polys.push_back(iin_box_polys[iip]);
               break;
      case  1: iinterior_polys.push_back(iin_box_polys[iip]);
               break;
      default: throw std::runtime_error("Undefined position");
    }
}

void append_data(RefPolySet_t& orig_polys,
                 const RefPolySet_t& new_polys) {
  int new_size = (int) (orig_polys.r3d_polys.size() + new_polys.r3d_polys.size());

  orig_polys.r3d_polys.reserve(new_size);
  orig_polys.ipoly2cellID.reserve(new_size);
  orig_polys.moments.reserve(new_size);

  orig_polys.r3d_polys.insert(orig_polys.r3d_polys.end(), 
    new_polys.r3d_polys.begin(), new_polys.r3d_polys.end());
  orig_polys.ipoly2cellID.insert(orig_polys.ipoly2cellID.end(),
    new_polys.ipoly2cellID.begin(), new_polys.ipoly2cellID.end());
  orig_polys.moments.insert(orig_polys.moments.end(), 
    new_polys.moments.begin(), new_polys.moments.end());   
}

void extract_polys(const RefPolySet_t& all_polys,
                   const std::vector<int>& set_ipolys,
                   RefPolySet_t& set_polys ) {

  int num_set_polys = set_ipolys.size();

  set_polys.r3d_polys.resize(num_set_polys);
  set_polys.ipoly2cellID.resize(num_set_polys);
  set_polys.moments.resize(num_set_polys);

  for (int isp = 0; isp < num_set_polys; isp++) {
    set_polys.r3d_polys[isp] = all_polys.r3d_polys[set_ipolys[isp]];
    set_polys.ipoly2cellID[isp] = all_polys.ipoly2cellID[set_ipolys[isp]];
    set_polys.moments[isp] = all_polys.moments[set_ipolys[isp]];
  }
}

void apply_poly(const Tangram::MatPoly<3>& mat_poly,
                RefPolySet_t& multimat_polys,
                RefPolySet_t& interior_polys,
                bool convex_matpoly = false) {
  assert(!multimat_polys.r3d_polys.empty());

  interior_polys.clear();

  if (convex_matpoly) {
    //We need to make sure that we only split r3d_polys that actually
    //intersect the MatPoly
    std::vector<int> iexterior_polys, iinstersecting_polys, iinterior_polys;
    sort_wrt_convex_poly(mat_poly, multimat_polys.r3d_polys, iexterior_polys, 
                         iinstersecting_polys, iinterior_polys);

    RefPolySet_t exterior_polys, intersecting_polys;

    extract_polys(multimat_polys, iexterior_polys, exterior_polys);
    extract_polys(multimat_polys, iinstersecting_polys, intersecting_polys);
    extract_polys(multimat_polys, iinterior_polys, interior_polys);

    multimat_polys.clear();

    std::vector< Tangram::Plane_t<3> > face_planes;
    mat_poly.face_planes(face_planes);

    int nfaces = (int) face_planes.size();

    //We cut with face planes slicing off exterior polys on each step
    for (int iface = 0; iface < nfaces; iface++) {
      RefPolySet_t new_exterior_polys;

      //Face normal is pointing outwards, but we need to keep slicing polys
      //below the plane
      Tangram::Plane_t<3> cutting_plane = -face_planes[iface];

      apply_plane(cutting_plane, intersecting_polys, new_exterior_polys);
      append_data(exterior_polys, new_exterior_polys);
    }

    append_data(interior_polys, intersecting_polys);    
    multimat_polys = exterior_polys;
  }
  else {
    std::vector< Tangram::MatPoly<3> > matpoly_tets;
    mat_poly.facetize_decompose(matpoly_tets);

    int ntets = (int) matpoly_tets.size();

    for (int itet = 0; itet < ntets; itet++) {
      //We need to make sure that we only split r3d_polys that actually
      //intersect the tets of MatPoly
      std::vector<int> iexterior_polys, iinstersecting_polys, iinterior_polys;
      sort_wrt_convex_poly(matpoly_tets[itet], multimat_polys.r3d_polys, 
                           iexterior_polys, iinstersecting_polys, iinterior_polys);
      
      RefPolySet_t tet_interior_polys, exterior_polys, intersecting_polys;
      extract_polys(multimat_polys, iexterior_polys, exterior_polys);
      extract_polys(multimat_polys, iinstersecting_polys, intersecting_polys);
      extract_polys(multimat_polys, iinterior_polys, tet_interior_polys);

      multimat_polys.clear();

      std::vector< Tangram::Plane_t<3> > face_planes;
      matpoly_tets[itet].face_planes(face_planes);

      //We cut with face planes slicing off exterior polys on each step
      for (int iface = 0; iface < 4; iface++) {
        RefPolySet_t new_exterior_polys;

        //Face normal is pointing outwards, but we need to keep slicing polys
        //below the plane
        Tangram::Plane_t<3> cutting_plane = -face_planes[iface];

        apply_plane(cutting_plane, intersecting_polys, new_exterior_polys);
        append_data(exterior_polys, new_exterior_polys);
      }

      append_data(tet_interior_polys, intersecting_polys); 

      //Poly is interior if it's interior for any one tets
      append_data(interior_polys, tet_interior_polys);
      
      //Poly is exterior if it's exterior for all the tets, so we keep
      //cutting out tets from the current exterior list
      multimat_polys = exterior_polys;
    }
  }
}

void finalize_ref_data(const std::vector<RefPolySet_t>& ref_poly_sets,
                       const std::vector<int>& sets_material_IDs,
                       std::vector<int>& cell_num_mats,
                       std::vector<int>& cell_mat_ids,
                       std::vector<double>& cell_mat_volfracs,
                       std::vector< Tangram::Point<3> >& cell_mat_centroids,
                       std::vector< std::vector< std::vector<r3d_poly> > >&
                         reference_mat_polys) {                     
  int ncells = -1, nsets = (int) ref_poly_sets.size();
  for (int iset = 0; iset < nsets; iset++) {
    int set_max_icell = *std::max_element(ref_poly_sets[iset].ipoly2cellID.begin(), 
                                          ref_poly_sets[iset].ipoly2cellID.end());
    if (set_max_icell > ncells) ncells = set_max_icell;  
  }
  ncells++;

  int nmoments = (int) ref_poly_sets[0].moments.size();

  std::vector< std::vector<int> > cells_mat_ids(ncells);
  std::vector< std::vector< std::vector<double> > > cells_mat_moments(ncells);
  reference_mat_polys.clear();
  reference_mat_polys.resize(ncells);

  for (int iset = 0; iset < nsets; iset++) {
    int cur_mat_id = sets_material_IDs[iset];

    for (int ipoly = 0; ipoly < ref_poly_sets[iset].r3d_polys.size(); ipoly++) {
      int icell = ref_poly_sets[iset].ipoly2cellID[ipoly];

      int cell_mat_id = std::distance(cells_mat_ids[icell].begin(), 
        std::find(cells_mat_ids[icell].begin(), cells_mat_ids[icell].end(), cur_mat_id));

      if (cell_mat_id == cells_mat_ids[icell].size()) {
        cells_mat_ids[icell].resize(cell_mat_id + 1);
        reference_mat_polys[icell].resize(cell_mat_id + 1);
        cells_mat_moments[icell].resize(cell_mat_id + 1);
        cells_mat_moments[icell][cell_mat_id].resize(nmoments, 0.0);
        cells_mat_ids[icell][cell_mat_id] = cur_mat_id;
      }

      reference_mat_polys[icell][cell_mat_id].push_back(
        ref_poly_sets[iset].r3d_polys[ipoly]);
      for (int im = 0; im < nmoments; im++)
        cells_mat_moments[icell][cell_mat_id][im] += 
          ref_poly_sets[iset].moments[ipoly][im];
    }
  }

  cell_num_mats.resize(ncells);
  cell_mat_ids.clear();
  cell_mat_volfracs.clear();
  cell_mat_centroids.clear();

  for (int icell = 0; icell < ncells; icell++) {
    cell_num_mats[icell] = cells_mat_ids[icell].size();
    cell_mat_ids.insert(cell_mat_ids.end(), cells_mat_ids[icell].begin(), 
                        cells_mat_ids[icell].end());

    double cell_mats_vol = 0.0;
    for (int icmat = 0; icmat < cell_num_mats[icell]; icmat++)
      cell_mats_vol += cells_mat_moments[icell][icmat][0];

    int offset = cell_mat_volfracs.size();
    cell_mat_volfracs.resize(offset + cell_num_mats[icell]);
    cell_mat_centroids.resize(offset + cell_num_mats[icell]);

    for (int icmat = 0; icmat < cell_num_mats[icell]; icmat++) {
      cell_mat_volfracs[offset + icmat] = 
        cells_mat_moments[icell][icmat][0]/cell_mats_vol;
      for (int idim = 0; idim < 3; idim++)
        cell_mat_centroids[offset + icmat][idim] = 
          cells_mat_moments[icell][icmat][idim + 1]/cells_mat_moments[icell][icmat][0];
    }
  }
}

#endif