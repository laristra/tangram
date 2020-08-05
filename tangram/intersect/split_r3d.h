/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#ifndef TANGRAM_INTERSECT_SPLIT_R3D_H_
#define TANGRAM_INTERSECT_SPLIT_R3D_H_

#include <vector>
#include <algorithm>
#include <memory>

// tangram includes
extern "C" {
#include "wonton/intersect/r3d/r3d.h"
}
#include "tangram/support/tangram.h"
#include "tangram/support/MatPoly.h"

/*!
 @file split_r3d.h
 @brief Routines for splitting convex or non-convex MatPoly's
 by a cutting plane. Uses R3D.
 A convex MatPoly will be split into two MatPoly's.
 A non-convex MatPoly will be decomposed into convex polyhedra using
 its centroid, several MatPoly's could be returned for each half-plane.
 All the faces of the input MatPoly are assumed to be planar and no
 additional facetization is performed, so when sequentially splitting of a
 faceted MatPoly, no additional facets will be introduced.
 */

namespace Tangram {

/*!
  @brief Converts a MatPoly to a polyhedron in R3D format.
  @param[in] mat_poly MatPoly object to convert
  @param[out] r3dpoly Corresponding R3D polyhedron
*/
inline
void
matpoly_to_r3dpoly(const MatPoly<3>& mat_poly, r3d_poly& r3dpoly) {

  //Translate coordinates of vertices to R3D format
  const std::vector<Point3>& matpoly_vrts = mat_poly.points();
  int nvrts = mat_poly.num_vertices();
  r3d_rvec3* r3dized_poly_vrts;
  r3dized_poly_vrts = new r3d_rvec3 [nvrts];
  for (int ivrt = 0; ivrt < nvrts; ivrt++)
    for (int ixyz = 0; ixyz < 3; ixyz++)
      r3dized_poly_vrts[ivrt].xyz[ixyz] = matpoly_vrts[ivrt][ixyz];
  
  //Get faces of the polyhedron in terms of local indices of its vertices
  int nfaces = mat_poly.num_faces();
  r3d_int* nface_vrts;
  nface_vrts = new r3d_int [nfaces];
  r3d_int** r3dized_poly_faces;
  r3dized_poly_faces = new r3d_int* [nfaces];
  for (int iface = 0; iface < nfaces; iface++) {
    const std::vector<int>& matpoly_faces = mat_poly.face_vertices(iface);
    int face_nvrts = static_cast<int>(matpoly_faces.size());
    nface_vrts[iface] = face_nvrts;
    r3dized_poly_faces[iface] = new r3d_int [face_nvrts];
    for (int ivrt = 0; ivrt < face_nvrts; ivrt++)
      r3dized_poly_faces[iface][ivrt] = matpoly_faces[ivrt];
  }
  
  r3d_init_poly(&r3dpoly, r3dized_poly_vrts, nvrts, r3dized_poly_faces, nface_vrts, nfaces);
  
  delete [] r3dized_poly_vrts;
  delete [] nface_vrts;
  for (int iface = 0; iface < nfaces; iface++)
    delete [] r3dized_poly_faces[iface];
  delete [] r3dized_poly_faces;
}

/*!
 @brief Breaks an r3d_poly into a vector of disjoint components.
 @param[in] r3dpoly Original r3dpoly that can have disjoint components
 @param[out] poly_components Vector of r3dpoly's, one for each components
 @param[in] vol_tol Volume tolerance, only components with volumes above it
 will be included in the output vector
*/
inline
void r3d_poly_components(const r3d_poly& r3dpoly, std::vector<r3d_poly>& poly_components,
                         const double vol_tol) {
  const int R3D_POLY_ORDER = 1;
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
    for (int icmp_vrt = 0; icmp_vrt < cur_cmp.nverts; icmp_vrt++) {
      for (int k = 0; k < 3; k++) {
        int old_id = cur_cmp.verts[icmp_vrt].pnbrs[k];
        cur_cmp.verts[icmp_vrt].pnbrs[k] = poly_irvt2cmp_ivrt[old_id];
      }
    }

    //Confirm that the component is not empty
    r3d_reduce(&cur_cmp, cmp_moments, R3D_POLY_ORDER);
    if (cmp_moments[0] >= vol_tol)
      poly_components.push_back(cur_cmp);

    //Find the first vertex of the next component
    for (icmp_first_vrt = 0; icmp_first_vrt < total_nvrts; icmp_first_vrt++)
      if (vrt_marks[icmp_first_vrt] == 0)
        break;
  }
}

/*!
  @brief Converts a polyhedron in R3D format to MatPoly's.
  @param[in] r3dpoly R3D polyhedron to convert
  @param[out] mat_polys Corresponding vector of MatPoly's: 
  will contain as many MatPoly's as there are components 
  in the R3D polyhedron
  @param[in] dst_tol Distance tolerance
  @param[in] reference_pts Preferred points to snap vertices to
*/
inline
void
r3dpoly_to_matpolys(const r3d_poly& r3dpoly,
                    std::vector< MatPoly<3> >& mat_polys,
                    const double vol_tol,
                    const double dst_tol,
                    const std::vector< Point<3> >* reference_pts = nullptr) {
  mat_polys.clear();
  std::vector<r3d_poly> poly_components; 
  r3d_poly_components(r3dpoly, poly_components, vol_tol);

  if (poly_components.empty())
    return;

  int nb_poly_components = poly_components.size();

  mat_polys.reserve(nb_poly_components);
  for (int ir3dpoly = 0; ir3dpoly < nb_poly_components; ir3dpoly++) {
    r3d_brep* poly_brep;
    r3d_int ncomponents;
    r3d_poly poly_copy = poly_components[ir3dpoly];
    r3d_init_brep(&poly_copy, &poly_brep, &ncomponents);

    if (ncomponents == 0) {
      r3d_free_brep(&poly_brep, 0);
      continue;
    }

    assert(ncomponents == 1);
    
    int nvrts = poly_brep[0].numvertices;
    std::vector< Point<3> > curpoly_vrts(nvrts);
    for (int ivrt = 0; ivrt < nvrts; ivrt++)
      for (int ixyz = 0; ixyz < 3; ixyz++)
        curpoly_vrts[ivrt][ixyz] = poly_brep[0].vertices[ivrt].xyz[ixyz];

    int nfaces = poly_brep[0].numfaces;
    std::vector< std::vector<int> > curpoly_faces(nfaces);
    for (int iface = 0; iface < nfaces; iface++) {
      int face_nverts = poly_brep[0].numvertsperface[iface];
      curpoly_faces[iface].resize(face_nverts);
      for (int ifv = 0; ifv < face_nverts; ifv++)
        curpoly_faces[iface][ifv] = poly_brep[0].faceinds[iface][ifv];
    }

    // Eliminate degeneracies that r3d has potentially introduced: vertices
    // within dst_tol from each other are replaced by a single vertex,
    // hanging nodes are identified and removed, degenerate faces with
    // less than three vertices are eliminated.
    MatPoly<3> fit_poly = natural_selection(curpoly_vrts, curpoly_faces, 
                                            dst_tol, reference_pts);
    // We do not store polyhedra with less than four faces, 
    // as they are clearly degenerate
    if (not fit_poly.points().empty()) {
      mat_polys.push_back(fit_poly);
    }

    mat_polys.shrink_to_fit();
    r3d_free_brep(&poly_brep, ncomponents);
  }
}

/*!
  @brief Splits a convex MatPoly into two (convex) MatPoly's
  with a cutting plane.
  @param[in] mat_poly Convex MatPoly to split
  @param[in] cutting_plane Cutting plane to split with
  @param[out] lower_halfspace_poly MatPoly below the plane
  @param[out] upper_halfspace_poly MatPoly above the plane
  @param[out] lower_halfspace_moments Moments of MatPoly below the plane
  @param[out] upper_halfspace_moments Moments of MatPoly above the plane
*/
inline
void
split_convex_matpoly_r3d(const MatPoly<3>& mat_poly,
                         const Plane_t<3>& cutting_plane,
                         MatPoly<3>& lower_halfspace_poly,
                         MatPoly<3>& upper_halfspace_poly,
                         std::vector<double>& lower_halfspace_moments,
                         std::vector<double>& upper_halfspace_moments,
                         double vol_tol, double dst_tol) {
  //Translate the cutting plane to R3D format
  r3d_plane r3d_cut_plane;
  for (int ixyz = 0; ixyz < 3; ixyz++)
    r3d_cut_plane.n.xyz[ixyz] = cutting_plane.normal[ixyz];
  r3d_cut_plane.d = cutting_plane.dist2origin;

  r3d_poly r3dized_poly;
  r3d_poly r3d_subpolys[2];
  matpoly_to_r3dpoly(mat_poly, r3dized_poly);
  r3d_split(&r3dized_poly, 1, r3d_cut_plane, &r3d_subpolys[1], &r3d_subpolys[0]);

  MatPoly<3>* subpoly_ptrs[2] = {&lower_halfspace_poly, &upper_halfspace_poly};
  std::vector<double>* subpoly_moments_ptrs[2] = {&lower_halfspace_moments, 
                                                  &upper_halfspace_moments};

  const int POLY_ORDER = 1;
  r3d_real r3d_moments[R3D_NUM_MOMENTS(POLY_ORDER)];
  int iempty_subpoly = -1;
  for (int isp = 0; isp < 2; isp++) {
    //Check if the subpoly is empty
    //Note that if one of the subpoly's is empty, then the other is the
    //original MatPoly, so we just need the ID (0 or 1) of the empty subpoly
    //to get the result    
    if (r3d_subpolys[isp].nverts == 0) {
      subpoly_ptrs[isp]->clear();
      subpoly_moments_ptrs[isp]->assign(4, 0.0);
      iempty_subpoly = isp;
      break;
    }

    //Find the moments for a subpoly
    r3d_reduce(&r3d_subpolys[isp], r3d_moments, POLY_ORDER);

    if (r3d_moments[0] < vol_tol) {
      subpoly_ptrs[isp]->clear();
      subpoly_moments_ptrs[isp]->assign(4, 0.0);
      iempty_subpoly = isp;
      break;
    }    

    subpoly_moments_ptrs[isp]->assign(r3d_moments, r3d_moments + 4);
  }

  //If there are no empty subpolys, we return two resulting pieces.
  //If one of the subpolys is empty, we return the original MatPoly to
  //ensure consistency and prevent unnecessary changes
  if (iempty_subpoly == -1) {
    for (int isp = 0; isp < 2; isp++) {
      //Get a MatPoly for a subpoly
      std::vector< MatPoly<3> > sub_matpoly;
      r3dpoly_to_matpolys(r3d_subpolys[isp], sub_matpoly, vol_tol, dst_tol, 
                          &mat_poly.points());
      int ncomponents = static_cast<int>(sub_matpoly.size());
      if (ncomponents > 1) {
        // Filter out degenerate components
        int subpoly_id = 0;
        while (static_cast<int>(sub_matpoly.size()) > subpoly_id) {
          if (sub_matpoly[subpoly_id].moments()[0] >= vol_tol) subpoly_id++;
          else sub_matpoly.erase(sub_matpoly.begin() + subpoly_id);
        }
      }
        
      if (ncomponents == 1) 
        *subpoly_ptrs[isp] = sub_matpoly[0];
      else
        throw std::runtime_error("Non-convex MatPoly is split using the method for convex MatPoly's!");
    }
  }
  else {
    int ifull_poly = (iempty_subpoly + 1)%2;
    matpoly_to_r3dpoly(mat_poly, r3d_subpolys[ifull_poly]);
    r3d_reduce(&r3d_subpolys[ifull_poly], r3d_moments, POLY_ORDER);
    subpoly_moments_ptrs[ifull_poly]->assign(r3d_moments, r3d_moments + 4);

    *subpoly_ptrs[ifull_poly] = mat_poly;
  }
}

/*!
 * \class SplitR3D  3-D splitting algorithm
 *
 * Splits a vector of MatPoly's with a cutting plane. The result is HalfSpaceSets_t
 * structure that contains two sets of convex MatPoly's: one corresponds to components below,
 * and the other to components above the plane. Each set has associated moments, which
 * are aggregated moments of all components in the set.
 * 
 * Important: input MatPoly's are assumed to have planar faces.
 * By default, input MatPoly's are assumed to be non-covex and are decomposed into
 * convex components. If all_convex is set to true in the constructor, all input
 * MatPoly's are assumed to be convex and will NOT be decomposed.
 */

class SplitR3D {
 public:
  SplitR3D(const std::vector< MatPoly<3> >& matpolys,
           const Plane_t<3>& cutting_plane, 
           const double vol_tol,
           const double dst_tol, 
           const bool all_convex) : 
           matpolys_(matpolys), 
           cutting_plane_(cutting_plane),
           vol_tol_(vol_tol),
           dst_tol_(dst_tol),
           all_convex_(all_convex) {}

  /*! 
    @brief Splits a MatPoly into two sets of convex MatPoly's
    @return Structure containing two sets of MatPoly's, one for each half-space,
    with respective aggregated moments.
  */
  HalfSpaceSets_t<3> operator() () const {
    HalfSpaceSets_t<3> hs_sets;

    std::vector< MatPoly<3> >* hs_subpolys_ptrs[2] = {
      &hs_sets.lower_halfspace_set.matpolys,
      &hs_sets.upper_halfspace_set.matpolys
    };

    std::vector<double>* hs_moments_ptrs[2] = {
      &hs_sets.lower_halfspace_set.moments,
      &hs_sets.upper_halfspace_set.moments
    };

    std::vector< MatPoly<3> > convex_components;
    const std::vector< MatPoly<3> >* convex_polys;
    if (all_convex_)
      convex_polys = &matpolys_;
    else {
      int nncpolys = static_cast<int>(matpolys_.size());
      for (int incp = 0; incp < nncpolys; incp++)
        matpolys_[incp].decompose(convex_components);

      convex_polys = &convex_components;
    }

    for (auto&& moment : hs_moments_ptrs)
      moment->assign(4, 0.0);

    int hs_poly_count[2] = {0, 0};
    int npolys = static_cast<int>(convex_polys->size());

    for (int icp = 0; icp < npolys; icp++) {
      MatPoly<3> cur_subpolys[2];
      std::vector<double> cur_moments[2];
      split_convex_matpoly_r3d((*convex_polys)[icp], cutting_plane_,
                               cur_subpolys[0], cur_subpolys[1],
                               cur_moments[0], cur_moments[1], vol_tol_, dst_tol_);

      int cur_face_group_id = (*convex_polys)[icp].face_group_id();
      for (int ihs = 0; ihs < 2; ihs++) {
        if (cur_subpolys[ihs].num_vertices() != 0) {
          if (cur_face_group_id >= 0)
            cur_subpolys[ihs].set_face_group_id(cur_face_group_id);          
          hs_subpolys_ptrs[ihs]->emplace_back(cur_subpolys[ihs]);
          for (int im = 0; im < 4; im++)
            (*hs_moments_ptrs[ihs])[im] += cur_moments[ihs][im];
          hs_poly_count[ihs]++;
        }
      }
    }

    for (int ihs = 0; ihs < 2; ihs++) {
      if (hs_poly_count[ihs] == 0)
        hs_moments_ptrs[ihs]->assign(4, 0.0);
    }

    return hs_sets;
  }

  //! Default constructor (disabled)
  SplitR3D() = delete;

  //! Assignment operator (disabled)
  SplitR3D& operator = (const SplitR3D&) = delete;

 private:
  const std::vector< MatPoly<3> >& matpolys_;
  const Plane_t<3>& cutting_plane_;
  double vol_tol_;
  double dst_tol_;  
  bool all_convex_;
};

/*!
 * \class ClipR3D  3-D clipping algorithm
 *
 * For a vector of MatPoly's, computes the aggregated moments of their 
 * parts below the cutting plane.
 * This can be interpreted as moments of chopped off parts in the nested 
 * dissections algorithm.
 * 
 * In this routine, we utilize the fact that R3D can clip a non-convex 
 * polyhedron with a plane.
 * 
 * By default, the faces of MatPoly's are assumed to be non-planar
 * and are facetized. If planar_faces is set to true,
 * faces of all MatPoly's are assumed to be planar and are NOT facetized.
 *
 */

class ClipR3D {
 public:
  ClipR3D(double vol_tol) : vol_tol_(vol_tol) {}

  /*! 
    @brief Set the cutting plane used by the functor. 
    @param[in] cutting_plane Cutting plane to clip with. 
  */
  void set_plane(const Plane_t<3>& cutting_plane) {
    for (int ixyz = 0; ixyz < 3; ixyz++)
      r3d_cut_plane_.n.xyz[ixyz] = -cutting_plane.normal[ixyz];
    r3d_cut_plane_.d = -cutting_plane.dist2origin;
  }

  /*! 
    @brief Set the material polyhedra to be used by the functor.
    @param[in] matpolys Material polyhedra to be clipped.
    @param[in] planar_faces Flag indicating whether the faces of
    material polyhedra are planar: if not, they are facetized.
  */
  void set_matpolys(const std::vector< MatPoly<3> >& matpolys,
                    const bool planar_faces) {
    int npolys = static_cast<int>(matpolys.size());
    r3d_polys_.resize(npolys);
    if (planar_faces) {
      for(int ipoly = 0; ipoly < npolys; ipoly++)
        matpoly_to_r3dpoly(matpolys[ipoly], r3d_polys_[ipoly]);
    }
    else {
      for(int ipoly = 0; ipoly < npolys; ipoly++) {
        MatPoly<3> faceted_poly;
        matpolys[ipoly].faceted_matpoly(&faceted_poly);
        matpoly_to_r3dpoly(faceted_poly, r3d_polys_[ipoly]);
      }
    }
  }

  /*! 
    @brief Computes moments of provided material polyhedra 
    @return Vector containing aggregated moments of all
    provided material polyhedra
  */
  std::vector<double> aggregated_moments() {
    std::vector<double> agg_moments;
    if (not r3d_polys_.empty()) {
      agg_moments.assign(4, 0.0);

      const int POLY_ORDER = 1;
      r3d_real r3d_moments[R3D_NUM_MOMENTS(POLY_ORDER)];

      for (auto&& poly : r3d_polys_) {
        if (poly.nverts != 0) {
          r3d_reduce(&poly, r3d_moments, POLY_ORDER);
          for (int im = 0; im < 4; im++)
            agg_moments[im] += r3d_moments[im];
        }
      }
    }

    if (agg_moments[0] < vol_tol_)
      agg_moments.assign(4, 0.0);

    return agg_moments;
  }

  /*! 
    @brief Computes moments of chopped off components: 
    the parts of MatPoly's below the plane
    @return Vector containing aggregated moments 
    of the chopped off part of all MatPoly's
  */
  std::vector<double> operator() () const {
    std::vector<double> below_plane_moments(4, 0.0);

    const int POLY_ORDER = 1;
    r3d_real r3d_moments[R3D_NUM_MOMENTS(POLY_ORDER)];

    for (auto&& poly : r3d_polys_) {
      if (poly.nverts != 0) {
        //r3d does in-place clipping and does not have the const modifier
        //for the plane: we need to make copies
        r3d_poly clipped_poly = poly;
        r3d_plane plane_copy = r3d_cut_plane_;
        r3d_clip(&clipped_poly, &plane_copy, 1);
        if (clipped_poly.nverts != 0) {
          //Find the moments for the part in the lower halfspace
          r3d_reduce(&clipped_poly, r3d_moments, POLY_ORDER);
          for (int im = 0; im < 4; im++)
            below_plane_moments[im] += r3d_moments[im];
        }
      }
    }
    
    if (below_plane_moments[0] < vol_tol_)
      below_plane_moments.assign(4, 0.0);
      
    return below_plane_moments;
  }

  //! Assignment operator (disabled)
  ClipR3D& operator = (const ClipR3D&) = delete;

 private:
  std::vector<r3d_poly> r3d_polys_;
  r3d_plane r3d_cut_plane_;
  double vol_tol_;  
};

/*!
  @brief Computes moments of the intersection of a MatPoly and an r3d_poly. 
  If MatPoly is convex, r3d_poly will be clipped with planes containing the
  faces of MatPoly. Otherwise, MatPoly will be decomposed into tetrahedral
  MatPolys, each of which will be intersected with r3d_poly. 
  Note that r3d_poly does not need to be convex.
  @param[in] mat_poly MatPoly object to intersect with
  @param[in] r3dpoly r3d_poly that is intersected with MatPoly
  @param[out] intersection_moments Moments of the intersection
  @param[in] convex_matpoly flag indicating if MatPoly is convex: if not
  it will be decomposed into tetrahedra
*//*!
  @brief Computes moments of the intersection of a MatPoly and an r3d_poly. 
  If MatPoly is convex, r3d_poly will be clipped with planes containing the
  faces of MatPoly. Otherwise, MatPoly will be decomposed into tetrahedral
  MatPolys, each of which will be intersected with r3d_poly. 
  Note that r3d_poly does not need to be convex.
  @param[in] mat_poly MatPoly object to intersect with
  @param[in] r3dpoly r3d_poly that is intersected with MatPoly
  @param[out] intersection_moments Moments of the intersection
  @param[in] convex_matpoly flag indicating if MatPoly is convex: if not
  it will be decomposed into tetrahedra
*/
inline
void get_intersection_moments(const MatPoly<3>& mat_poly,
                              const r3d_poly& r3dpoly,
                              std::vector<double>& intersection_moments,
                              bool convex_matpoly) {   
  const int POLY_ORDER = 1;
  r3d_real r3d_moments[R3D_NUM_MOMENTS(POLY_ORDER)];

  if (convex_matpoly) {                                
    std::vector< Plane_t<3> > face_planes;
    mat_poly.face_planes(face_planes);
    r3d_poly intersection = r3dpoly;

    int nplanes = static_cast<int>(face_planes.size());
    auto* r3d_face_planes = new r3d_plane[nplanes];

    for (int iplane = 0; iplane < nplanes; iplane++) {
      for (int ixyz = 0; ixyz < 3; ixyz++) {
        r3d_face_planes[iplane].n.xyz[ixyz] = -face_planes[iplane].normal[ixyz];
      }
      r3d_face_planes[iplane].d = -face_planes[iplane].dist2origin;
    }

    r3d_clip(&intersection, r3d_face_planes, (r3d_int) nplanes);
    
    delete [] r3d_face_planes;

    r3d_reduce(&intersection, r3d_moments, POLY_ORDER);
    intersection_moments.resize(4);

    for (int im = 0; im < 4; im++)
      intersection_moments[im] = r3d_moments[im];
  }
  else {
    std::vector< MatPoly<3> > mat_poly_tets;
    mat_poly.facetize_decompose(mat_poly_tets);

    int ntets = static_cast<int>(mat_poly_tets.size());
    intersection_moments.assign(4, 0.0);

    for (int itet = 0; itet < ntets; itet++) {
      std::vector< Plane_t<3> > face_planes;
      mat_poly_tets[itet].face_planes(face_planes);

      if(face_planes.size() < 4) continue;

      r3d_poly intersection = r3dpoly; 
      r3d_plane r3d_face_planes[4];

      for (int iplane = 0; iplane < 4; iplane++) {
        for (int ixyz = 0; ixyz < 3; ixyz++) {
          r3d_face_planes[iplane].n.xyz[ixyz] = -face_planes[iplane].normal[ixyz];
        }
        r3d_face_planes[iplane].d = -face_planes[iplane].dist2origin;
      }

      r3d_clip(&intersection, r3d_face_planes, 4);

      r3d_reduce(&intersection, r3d_moments, POLY_ORDER);
      for (int im = 0; im < 4; im++)
        intersection_moments[im] += r3d_moments[im];
    }
  }
}

inline
BoundingBox_t<3> r3d_poly_bounding_box(const r3d_poly& r3dpoly) {
  BoundingBox_t<3> bbox;
  for (int ivrt = 0; ivrt < r3dpoly.nverts; ivrt++) {
    for (int idim = 0; idim < 3; idim++) {
      if (r3dpoly.verts[ivrt].pos.xyz[idim] < bbox.min[idim])
        bbox.min[idim] = r3dpoly.verts[ivrt].pos.xyz[idim];

      if (r3dpoly.verts[ivrt].pos.xyz[idim] > bbox.max[idim])
        bbox.max[idim] = r3dpoly.verts[ivrt].pos.xyz[idim];
    }
  }

  return bbox;
}

}  // namespace Tangram

#endif  // TANGRAM_INTERSECT_SPLIT_R3D_H_
