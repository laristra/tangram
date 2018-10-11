/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#ifndef RPG_TOOLS_CUTS_H_
#define RPG_TOOLS_CUTS_H_

#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <memory>

extern "C" {
#include "tangram/intersect/r3d.h"
}

#include "tangram/support/tangram.h"
#include "tangram/intersect/split_r3d.h"

#define BUGS_IN_R3D 1

#define R3D_POLY_ORDER 1

/* Functions for intersection-based generation of reference
   material polyhedra. Given data on a collection of polyhedra
   generates sub-collections by splitting the original polyhedra
   with a plane or cutting out a polyhedral shape.
   Uses r3d.
   Collection of single-material polyhedra sets can be 
   converted to the format used by the Tangram's driver */

namespace R3DPOLY {
  enum Position {
    INSIDE, INTERSECTS, OUTSIDE
  };
}

struct RefPolyData_t {
  r3d_poly r3dpoly;   // geometry
  int cellID;         // ID of the containing cell
  r3d_real moments[R3D_NUM_MOMENTS(R3D_POLY_ORDER)]; // moments
};

/*!
 @brief Breaks an r3d_poly into a vector of disjoint r3d_poly's.
*/
void r3d_poly_components(const r3d_poly& r3dpoly, std::vector<r3d_poly>& poly_components) {
  poly_components.clear();

  int total_nvrts = static_cast<int>(r3dpoly.nverts);
  std::vector<int> vrt_marks(total_nvrts, 0);
  std::vector< std::vector<bool> > walked_edge(total_nvrts);
  for (int ivrt = 0; ivrt < total_nvrts; ivrt++)
    walked_edge[ivrt].assign(3, false);

  std::vector<bool> vrt_in_cmp(total_nvrts, false);
  std::vector<int> poly_irvt2cmp_ivrt(total_nvrts, -1);

  int icmp_first_vrt = 0;
  r3d_real cmp_moments[R3D_NUM_MOMENTS(R3D_POLY_ORDER)];

  while(icmp_first_vrt < total_nvrts) {
    r3d_poly cur_cmp;
    cur_cmp.nverts = 0;

    //Walk all faces of the current components
    int iface_first_vrt = icmp_first_vrt;
    int iedge = 0;
    while (iface_first_vrt < total_nvrts) {
      int icur_vrt = iface_first_vrt;

      for (iedge = 0; iedge < 3; iedge++)
        if (!walked_edge[iface_first_vrt][iedge])
          break;

      //Walk the current face
      do {
        vrt_marks[icur_vrt]++;
        
        if (!vrt_in_cmp[icur_vrt]) {
          poly_irvt2cmp_ivrt[icur_vrt] = cur_cmp.nverts;
          cur_cmp.nverts++;
          cur_cmp.verts[poly_irvt2cmp_ivrt[icur_vrt]] = r3dpoly.verts[icur_vrt];
          vrt_in_cmp[icur_vrt] = true;
        }

        int inext_vrt = r3dpoly.verts[icur_vrt].pnbrs[iedge];
        walked_edge[icur_vrt][iedge] = true;

        int ireturn_edge;
        for (ireturn_edge = 0; ireturn_edge < 3; ireturn_edge++)
          if (r3dpoly.verts[inext_vrt].pnbrs[ireturn_edge] == icur_vrt)
            break;

        iedge = (ireturn_edge + 2)%3;
        icur_vrt = inext_vrt;

      } while (icur_vrt != iface_first_vrt);

      for (iface_first_vrt = 0; iface_first_vrt < total_nvrts; iface_first_vrt++)
        if( (vrt_marks[iface_first_vrt] > 0) && (vrt_marks[iface_first_vrt] < 3) )
          break;
    }

    //Neighbors data still uses r3dpoly indices
    for (int icmp_vrt = 0; icmp_vrt < cur_cmp.nverts; icmp_vrt++)
      for (int iedge = 0; iedge < 3; iedge++) {
        int old_id = cur_cmp.verts[icmp_vrt].pnbrs[iedge];
        cur_cmp.verts[icmp_vrt].pnbrs[iedge] = poly_irvt2cmp_ivrt[old_id];
      }

    //Confirm that the component is not empty
    r3d_reduce(&cur_cmp, cmp_moments, R3D_POLY_ORDER);
    if (cmp_moments[0] > std::numeric_limits<double>::epsilon())
      poly_components.push_back(cur_cmp);

    //Find the first vertex of the next component
    for (icmp_first_vrt = 0; icmp_first_vrt < total_nvrts; icmp_first_vrt++)
      if (vrt_marks[icmp_first_vrt] == 0)
        break;
  }
}

/*!
 @brief For a given data on a collection of polyhedra find their positions with
 respect to a convex shape given by a MatPoly object.
*/
class r3d_poly_intersect_check {
public:
  /*!
  @brief Constructor.

  @param[in] polys_data Data on polyhedra
  @param[in] IDs_to_check IDs of polyhedra for which their positions 
  with respect to MatPoly is to be determined
  @param[in] convex_matpoly Convex shape with respect to which polyhedra are tested
  */
  r3d_poly_intersect_check(const std::vector< std::shared_ptr<RefPolyData_t> >& polys_data,
                           const std::vector<int>& IDs_to_check,
                           const Tangram::MatPoly<3>& convex_matpoly) :
                           polys_data_(polys_data), IDs_to_check_(IDs_to_check),
                           convex_matpoly_(convex_matpoly) {
    std::vector< Tangram::Plane_t<3> > fplanes;
    convex_matpoly_.face_planes(fplanes);

    int nplanes = static_cast<int>(fplanes.size());
    face_planes_.resize(nplanes);
    //Face planes, used to find the intersection with convex_matpoly
    for (int iplane = 0; iplane < nplanes; iplane++) {
      for (int ixyz = 0; ixyz < 3; ixyz++)
        face_planes_[iplane].n.xyz[ixyz] = -fplanes[iplane].normal[ixyz];
      face_planes_[iplane].d = -fplanes[iplane].dist2origin;
    }
  }

  /*!
  @brief Finds the position of a polyhedron with respect to the shape
  given by convex_matpoly_.

  @param[in] checkID ID of the polyhedra to test in the list of IDs of polyhedra
  for which their positions is to be determined
  @return position with respect to the shape given by convex_matpoly_
  */
  R3DPOLY::Position operator()(const int checkID) const {
    int polyID = IDs_to_check_[checkID];
    const r3d_poly& r3dpoly = polys_data_[polyID]->r3dpoly;

    bool intersects = false;
    for (int ivrt = 0; ivrt < r3dpoly.nverts; ivrt++) {
      Tangram::Point3 cur_vrt;
      for (int ixyz = 0; ixyz < 3; ixyz++)
        cur_vrt[ixyz] = r3dpoly.verts[ivrt].pos.xyz[ixyz];

      if (Tangram::point_inside_matpoly(convex_matpoly_, cur_vrt, true)) {
        if (!intersects && (ivrt > 0)) return R3DPOLY::Position::INTERSECTS;
        intersects = true;
      }
      else if (intersects) return R3DPOLY::Position::INTERSECTS;  
    }

    if (intersects) return R3DPOLY::Position::INSIDE;

    //Instead of checking if MatPoly vertices are inside the r3d_poly,
    //we look at the volume of the actual intersection
    r3d_poly intersection = r3dpoly;
    for (int iplane = 0; iplane < face_planes_.size(); iplane++) {
      r3d_clip(&intersection, &face_planes_[iplane], 1);
      if (intersection.nverts == 0) return R3DPOLY::Position::OUTSIDE;
    }
    
    return R3DPOLY::Position::INTERSECTS;
  }

private:
  const std::vector< std::shared_ptr<RefPolyData_t> >& polys_data_;
  const std::vector<int>& IDs_to_check_;
  const Tangram::MatPoly<3>& convex_matpoly_;
  std::vector<r3d_plane> face_planes_;
};

struct PolyHalfspaceData_t {
  std::shared_ptr<RefPolyData_t> lower;
  std::shared_ptr<RefPolyData_t> upper;
};

/*!
 @brief For a given data on a collection of polyhedra and a cutting plane,
 splits polyhedra with the plane and returns the collections for the lower 
 and the upper half-spaces.
*/
class r3d_split_operator {
public:
  /*!
  @brief Constructor.

  @param[in] polys_data Data on polyhedra to split
  @param[in] cutting_plane Cutting plane with respect to which half-spaces are defined
  */
  r3d_split_operator(const std::vector< std::shared_ptr<RefPolyData_t> >& polys_data,
                     const r3d_plane& cutting_plane) :
                     polys_data_(polys_data), cutting_plane_(cutting_plane) {}

  /*!
  @brief Splits a polyhedron with the cutting plane and returns data on its parts
  in the lower and the upper half-spaces.

  @param[in] polyID ID of the polyhedron in the vector of polyhedra to split
  @return data on respective parts in the lower and the upper half-spaces
  */
  PolyHalfspaceData_t operator()(const int polyID) const {
    PolyHalfspaceData_t hs_data = {.lower = std::make_shared<RefPolyData_t>(), 
                                   .upper = std::make_shared<RefPolyData_t>() };

    const r3d_poly& r3dpoly = polys_data_[polyID]->r3dpoly;
    const r3d_real* poly_moments = polys_data_[polyID]->moments;
    int cellID = polys_data_[polyID]->cellID;

    r3d_poly split_poly = r3dpoly;
    r3d_split(&split_poly, 1, cutting_plane_, &hs_data.upper->r3dpoly, &hs_data.lower->r3dpoly);

#if BUGS_IN_R3D
    std::shared_ptr<RefPolyData_t> hs_data_ptrs[2] = {hs_data.lower, hs_data.upper};
    for (int ihs = 0; ihs < 2; ihs++) {
      if (hs_data_ptrs[ihs]->r3dpoly.nverts > 0) {
        std::vector<r3d_poly> poly_components;
        r3d_poly_components(hs_data_ptrs[ihs]->r3dpoly, poly_components);
        if (poly_components.empty())
          hs_data_ptrs[ihs]->r3dpoly.nverts = 0;
        else {
          assert(poly_components.size() == 1);
          if (hs_data_ptrs[ihs]->r3dpoly.nverts != poly_components[0].nverts)
            hs_data_ptrs[ihs]->r3dpoly = poly_components[0];
          hs_data_ptrs[ihs]->cellID = cellID;
          r3d_reduce(&hs_data_ptrs[ihs]->r3dpoly, hs_data_ptrs[ihs]->moments, R3D_POLY_ORDER);
        }
      }
    }
#else
    bool nnz_cutoff = false;
    if (hs_data.lower->r3dpoly.nverts > 0) {
      r3d_reduce(&hs_data.lower->r3dpoly, hs_data.lower->moments, R3D_POLY_ORDER);
      if (hs_data.lower->moments[0] > std::numeric_limits<double>::epsilon()) {
        hs_data.lower->cellID = cellID;
        nnz_cutoff = true;
      }
      else
        hs_data.lower->r3dpoly.nverts = 0;
    }

    if (hs_data.upper->r3dpoly.nverts > 0) {
      if (nnz_cutoff) {
        for (int im = 0; im < R3D_NUM_MOMENTS(R3D_POLY_ORDER); im++)
          hs_data.upper->moments[im] = 
            poly_moments[im] - hs_data.lower->moments[im];
      } 
      else {
        for (int im = 0; im < R3D_NUM_MOMENTS(R3D_POLY_ORDER); im++)
          hs_data.upper->moments[im] = poly_moments[im];
      }
      if (hs_data.upper->moments[0] > std::numeric_limits<double>::epsilon())
        hs_data.upper->cellID = cellID;
      else  
        hs_data.upper->r3dpoly.nverts = 0;
    }
#endif
    return hs_data;
  }

private:
  const std::vector< std::shared_ptr<RefPolyData_t> >& polys_data_;
  const r3d_plane& cutting_plane_;
};

/*!
 @brief Generates data for polyhedra corresponding to given mesh cells.
 @tparam Mesh_Wrapper A lightweight wrapper to a specific input mesh
                      implementation that provides required functionality

 @param[in] mesh Mesh wrapper
 @param[out] polys_data Data on polyhedra corresponding to mesh cells
 @param[in] decompose_cells If set to true, mesh cells will be split into
 tetrahedra using their centroids, data on each tetrahedron will be added 
 to the collection
*/
template <class Mesh_Wrapper>
void mesh_to_r3d_polys(const Mesh_Wrapper& mesh,
                       std::vector< std::shared_ptr<RefPolyData_t> >& polys_data,
                       bool decompose_cells = true) {
  int ncells = mesh.num_owned_cells() + mesh.num_ghost_cells();

  polys_data.clear();

  for (int icell = 0; icell < ncells; icell++) {
    Tangram::MatPoly<3> mat_poly;
    Tangram::cell_get_matpoly(mesh, icell, &mat_poly);

    std::vector< Tangram::MatPoly<3> > cell_polys;
    if (decompose_cells)
      mat_poly.facetize_decompose(cell_polys);
    else
      cell_polys.push_back(mat_poly);

    int num_polys = static_cast<int>(cell_polys.size());
    int offset = static_cast<int>(polys_data.size());
    polys_data.resize(offset + num_polys);

    for (int icp = 0; icp < num_polys; icp++) {
      polys_data[offset + icp] = std::make_shared<RefPolyData_t>();
      polys_data[offset + icp]->cellID = icell;

      Tangram::matpoly_to_r3dpoly(cell_polys[icp], 
                                  polys_data[offset + icp]->r3dpoly);

      r3d_reduce(&polys_data[offset + icp]->r3dpoly, 
                 polys_data[offset + icp]->moments, R3D_POLY_ORDER);
    }
  }
}

/*!
 @brief Generates data for polyhedra in the lower and the upper half-spaces
 of a given plane. Can use multiple threads.

 @param[in] polys_data Data for the polyhedra to split
 @param[in] planar_interface Cutting plane defining the half-spaces
 @param[out] lower_hs_data Data on polyhedra in the lower half-space
 @param[out] upper_hs_data Data on polyhedra in the upper half-space
*/
void apply_plane(const std::vector< std::shared_ptr<RefPolyData_t> >& polys_data,
                 const Tangram::Plane_t<3>& planar_interface,
                 std::vector< std::shared_ptr<RefPolyData_t> >& lower_hs_data,
                 std::vector< std::shared_ptr<RefPolyData_t> >& upper_hs_data) {
  lower_hs_data.clear();
  upper_hs_data.clear();
  
  int num_polys = static_cast<int>(polys_data.size());
  int nmoments = R3D_NUM_MOMENTS(R3D_POLY_ORDER);

  r3d_plane r3d_planar_iface;
  for (int ixyz = 0; ixyz < 3; ixyz++)
    r3d_planar_iface.n.xyz[ixyz] = planar_interface.normal[ixyz];
  r3d_planar_iface.d = planar_interface.dist2origin;

  r3d_split_operator r3d_split_op(polys_data, r3d_planar_iface);
  Tangram::vector< PolyHalfspaceData_t > hs_data(num_polys);

  Tangram::transform(Tangram::make_counting_iterator(0), 
                     Tangram::make_counting_iterator(num_polys),
                     hs_data.begin(), r3d_split_op);

  for (int ipoly = 0; ipoly < num_polys; ipoly++) {
    const PolyHalfspaceData_t& poly_hs_data = hs_data[ipoly];
    if (poly_hs_data.lower->r3dpoly.nverts > 0) 
      lower_hs_data.push_back(poly_hs_data.lower);

    if (poly_hs_data.upper->r3dpoly.nverts > 0) 
      upper_hs_data.push_back(poly_hs_data.upper);    
  }
}

/*!
 @brief Checks the positions of polyhedra with respect to a given convex shape.
 Can use multiple threads.

 @param[in] polys_data Data on polyhedra for which their position is to be checked
 @param[in] convex_matpoly Shape with respect to which polyhedra are tested
 @param[out] iinterior_polys IDs of polyhedra inside the shape
 @param[out] iinstersecting_polys IDs of polyhedra intersected by the boundary of the shape
 @param[out] iexterior_polys IDs of polyhedra outside of the shape
*/
void sort_wrt_convex_poly(const std::vector< std::shared_ptr<RefPolyData_t> >& polys_data,
                          const Tangram::MatPoly<3>& convex_matpoly,
                          std::vector<int>& iinterior_polys,
                          std::vector<int>& iinstersecting_polys,
                          std::vector<int>& iexterior_polys ) {
  std::vector< Tangram::Plane_t<3> > face_planes;
  convex_matpoly.face_planes(face_planes);
  Tangram::BoundingBox_t<3> poly_bbox = convex_matpoly.bounding_box();

  int nfaces = static_cast<int>(face_planes.size());
  int npolys = static_cast<int>(polys_data.size()); 

  iexterior_polys.clear();
  std::vector<int> iin_box_polys;
  for (int ipoly = 0; ipoly < npolys; ipoly++) {
    Tangram::BoundingBox_t<3> cur_bbox = 
      Tangram::r3d_poly_bounding_box(polys_data[ipoly]->r3dpoly);
    if(overlapping_boxes(cur_bbox, poly_bbox))
      iin_box_polys.push_back(ipoly);
    else
      iexterior_polys.push_back(ipoly);
  }

  //Sort out the polys inside the box
  int num_in_box = static_cast<int>(iin_box_polys.size());

  r3d_poly_intersect_check r3d_isect_check(polys_data, iin_box_polys, convex_matpoly);
  Tangram::vector<R3DPOLY::Position> check_result(num_in_box);

  Tangram::transform(Tangram::make_counting_iterator(0), 
                     Tangram::make_counting_iterator(num_in_box),
                     check_result.begin(), r3d_isect_check);

  iinterior_polys.clear();
  iinstersecting_polys.clear();
  for (int iip = 0; iip < num_in_box; iip++)
    switch (check_result[iip]) {
      case R3DPOLY::Position::OUTSIDE: 
               iexterior_polys.push_back(iin_box_polys[iip]);
               break;
      case R3DPOLY::Position::INTERSECTS: 
               iinstersecting_polys.push_back(iin_box_polys[iip]);
               break;
      case R3DPOLY::Position::INSIDE: 
               iinterior_polys.push_back(iin_box_polys[iip]);
               break;
      default: throw std::runtime_error("Undefined position");
    }
}

/*!
 @brief Given data on polyhedra, generates its subset based on given IDs.
 
 @param[in] polys_data Data on polyhedra
 @param[in] set_ipolys IDs of polyhedra for which their data is to be extracted
 @param[out] set_data data subset for polyhedra with given IDs
*/
void extract_data(const std::vector< std::shared_ptr<RefPolyData_t> >& polys_data,
                  const std::vector<int>& set_ipolys,
                  std::vector< std::shared_ptr<RefPolyData_t> >& set_data ) {

  int num_set_polys = set_ipolys.size();
  set_data.clear();
  set_data.reserve(num_set_polys);

  for (int isp = 0; isp < num_set_polys; isp++)
    set_data.push_back(polys_data[set_ipolys[isp]]);
}

/*!
 @brief Generates data for polyhedra inside and outside
 of a given shape. Can use multiple threads.

 @param[in] polys_data Data for the polyhedra to sort and split
 @param[in] mat_poly Shape which partitions polyhedra
 @param[out] interior_data Data on polyhedra inside the shape
 @param[out] exterior_data Data on polyhedra outside of the shape
 @param[in] convex_matpoly Specifies if the shape is convex: if set
 to false, the shape will be partitioned into tetrahedra using its centroid
*/
void apply_poly(const std::vector< std::shared_ptr<RefPolyData_t> >& polys_data,
                const Tangram::MatPoly<3>& mat_poly,
                std::vector< std::shared_ptr<RefPolyData_t> >& interior_data,
                std::vector< std::shared_ptr<RefPolyData_t> >& exterior_data,
                bool convex_matpoly = false) {
  if (convex_matpoly) {
    //We need to make sure that we only split r3d_polys that actually
    //intersect the MatPoly
    std::vector<int> iexterior_polys, iinstersecting_polys, iinterior_polys;
    sort_wrt_convex_poly(polys_data, mat_poly, iinterior_polys, 
                         iinstersecting_polys, iexterior_polys);

    extract_data(polys_data, iinterior_polys, interior_data);
    extract_data(polys_data, iexterior_polys, exterior_data);

    std::vector< std::shared_ptr<RefPolyData_t> > intersection_data;
    extract_data(polys_data, iinstersecting_polys, intersection_data);
    
    std::vector< Tangram::Plane_t<3> > face_planes;
    mat_poly.face_planes(face_planes);

    int nfaces = static_cast<int>(face_planes.size());

    //We cut with face planes slicing off exterior polys on each step
    for (int iface = 0; iface < nfaces; iface++) {
      std::vector< std::shared_ptr<RefPolyData_t> > upd_intersection_data, 
                                                    new_exterior_data;
                                                  
      //Face normal is pointing outward
      apply_plane(intersection_data, face_planes[iface], 
                  upd_intersection_data, new_exterior_data);

      intersection_data = upd_intersection_data;
      exterior_data.reserve(exterior_data.size() + new_exterior_data.size());
      exterior_data.insert(exterior_data.end(), 
                           new_exterior_data.begin(), new_exterior_data.end());
    }

    interior_data.reserve(interior_data.size() + intersection_data.size());
    interior_data.insert(interior_data.end(), 
                         intersection_data.begin(), intersection_data.end());
  }
  else {
    interior_data.clear();

    std::vector< Tangram::MatPoly<3> > matpoly_tets;
    mat_poly.facetize_decompose(matpoly_tets);

    int ntets = static_cast<int>(matpoly_tets.size());

    std::vector< std::shared_ptr<RefPolyData_t> > remaining_data = polys_data;
    for (int itet = 0; itet < ntets; itet++) {
      //We need to make sure that we only split r3d_polys that actually
      //intersect the tets of MatPoly
      std::vector<int> iexterior_polys, iinstersecting_polys, iinterior_polys;
      sort_wrt_convex_poly(remaining_data, matpoly_tets[itet], iinterior_polys,
                           iinstersecting_polys, iexterior_polys);
      
      std::vector< std::shared_ptr<RefPolyData_t> > tet_interior_data, 
                                                    intersection_data;

      extract_data(remaining_data, iexterior_polys, exterior_data);
      extract_data(remaining_data, iinstersecting_polys, intersection_data);

      extract_data(remaining_data, iinterior_polys, tet_interior_data);
      interior_data.reserve(interior_data.size() + tet_interior_data.size());
      interior_data.insert(interior_data.end(), 
                           tet_interior_data.begin(), tet_interior_data.end());

      std::vector< Tangram::Plane_t<3> > face_planes;
      matpoly_tets[itet].face_planes(face_planes);

      //We cut with face planes slicing off exterior polys on each step
      for (int iface = 0; iface < 4; iface++) {
        std::vector< std::shared_ptr<RefPolyData_t> > upd_intersection_data, 
                                                      new_exterior_data;
        //Face normal is pointing outward
        apply_plane(intersection_data, face_planes[iface], 
                    upd_intersection_data, new_exterior_data);

        intersection_data = upd_intersection_data;
        exterior_data.reserve(exterior_data.size() + new_exterior_data.size());
        exterior_data.insert(exterior_data.end(), 
                             new_exterior_data.begin(), new_exterior_data.end());
      }

      //Poly is interior if it's interior for any one tets
      interior_data.reserve(interior_data.size() + intersection_data.size());
      interior_data.insert(interior_data.end(), 
                           intersection_data.begin(), intersection_data.end());
          
      //Poly is exterior if it's exterior for all the tets, so we keep
      //cutting out tets from the current exterior list
      remaining_data = exterior_data;
    }
  }
}

unsigned int factorial(unsigned int n)
{
  unsigned int res = 1;
  for(unsigned int i = 0; i < n; ++i)
      res *= i + 1;
  return res;
}

/*!
 @brief For a given collections of single-material reference polyhedra sets
 generates data compatible with Tangram's driver. 

 @param[in] ref_sets_data Collection of data for single-material polyhedra sets
 @param[in] sets_material_IDs Material IDs for every set in the collection
 @param[out] cell_num_mats Number of material in each mesh cell, vector of length cell_num
 @param[out] cell_mat_ids Indices of materials in each mesh cell, a flat vector, requires
                          computations of offsets
 @param[out] cell_mat_volfracs Volume fractions of materials in each mesh cell, a flat
                               vector, requires computations of offsets
 @param[out] cell_mat_centroids Centroids of materials in each mesh cell, a flat vector,
                                requires computations of offsets
 @param[out] reference_mat_polys For every mesh cell and every material inside that cell, 
 the collection of single-material polyhedra containing that material          
 @param[in] permute_order If true, the order of materials in every cell will be 
 additionally permuted a random number of times                     
*/
void finalize_ref_data(const std::vector< std::vector< std::shared_ptr<RefPolyData_t> > >& 
                         ref_sets_data,
                       const std::vector<int>& sets_material_IDs,
                       std::vector<int>& cell_num_mats,
                       std::vector<int>& cell_mat_ids,
                       std::vector<double>& cell_mat_volfracs,
                       std::vector< Tangram::Point<3> >& cell_mat_centroids,
                       std::vector< std::vector< std::vector<r3d_poly> > >&
                         reference_mat_polys,
                       bool permute_order=false) {                     
  int ncells = -1, nsets = static_cast<int>(ref_sets_data.size());
  for (int iset = 0; iset < nsets; iset++) {
    int set_max_cellID = -1;
    for (int ipoly = 0; ipoly < ref_sets_data[iset].size(); ipoly++)
      if (ref_sets_data[iset][ipoly]->cellID > set_max_cellID)
        set_max_cellID = ref_sets_data[iset][ipoly]->cellID;
    
    if (set_max_cellID > ncells) ncells = set_max_cellID;  
  }
  ncells++;

  int nmoments = R3D_NUM_MOMENTS(R3D_POLY_ORDER);

  std::vector< std::vector<int> > cells_mat_ids(ncells);
  std::vector< std::vector< std::vector<double> > > cells_mat_moments(ncells);
  reference_mat_polys.clear();
  reference_mat_polys.resize(ncells);

  for (int iset = 0; iset < nsets; iset++) {
    int cur_mat_id = sets_material_IDs[iset];

    for (int ipoly = 0; ipoly < ref_sets_data[iset].size(); ipoly++) {
      int icell = ref_sets_data[iset][ipoly]->cellID;

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
        ref_sets_data[iset][ipoly]->r3dpoly);
      for (int im = 0; im < nmoments; im++)
        cells_mat_moments[icell][cell_mat_id][im] += 
          ref_sets_data[iset][ipoly]->moments[im];
    }
  }

  cell_num_mats.resize(ncells);
  cell_mat_ids.clear();
  cell_mat_volfracs.clear();
  cell_mat_centroids.clear();

  if (permute_order)
    srand(20150420);

  for (int icell = 0; icell < ncells; icell++) {
    int ncmats = static_cast<int>(cells_mat_ids[icell].size());
    cell_num_mats[icell] = ncmats;

    // If requested, we permuate the order of materials a random number of times
    // This can be used for testing accuracy of material order dependent methods
    if (permute_order) {
      int max_npermutations = factorial(ncmats);
      int cur_npermutations = rand()%max_npermutations;

      if (cur_npermutations != 0) {
        std::vector<int> new_mat_order(ncmats);
        std::iota(new_mat_order.begin(), new_mat_order.end(), 0);
        for (int ip = 0; ip < cur_npermutations; ip++)
          std::next_permutation(new_mat_order.begin(), new_mat_order.end());

        std::vector<int> new_cell_mat_ids(ncmats);
        std::vector< std::vector<double> > new_cell_mat_moments(ncmats);
        std::vector< std::vector<r3d_poly> > new_cell_ref_mat_polys(ncmats);
        for (int icmat = 0; icmat <ncmats; icmat++) {
          new_cell_mat_ids[icmat] = cells_mat_ids[icell][new_mat_order[icmat]];
          new_cell_mat_moments[icmat] = cells_mat_moments[icell][new_mat_order[icmat]];
          new_cell_ref_mat_polys[icmat] = reference_mat_polys[icell][new_mat_order[icmat]];
        }
        cells_mat_ids[icell] = new_cell_mat_ids;
        cells_mat_moments[icell] = new_cell_mat_moments;
        reference_mat_polys[icell] = new_cell_ref_mat_polys;
      }
    }

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