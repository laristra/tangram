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
void
matpoly_to_r3dpoly(const MatPoly<3>& mat_poly,
                   r3d_poly& r3dpoly) {
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
  @brief Converts a polyhedron in R3D format to MatPoly's.
  @param[in] r3dpoly R3D polyhedron to convert
  @param[out] mat_polys Corresponding vector of MatPoly's: 
  will contain as many MatPoly's as there are components 
  in the R3D polyhedron
*/
void
r3dpoly_to_matpolys(const r3d_poly& r3dpoly,
                    std::vector< MatPoly<3> >& mat_polys) {
  r3d_brep* poly_brep;
  r3d_int ncomponents;
  r3d_poly poly_copy = r3dpoly;
  r3d_init_brep(&poly_copy, &poly_brep, &ncomponents);

  mat_polys.clear();
  if (ncomponents == 0) {
    r3d_free_brep(&poly_brep, 0);
    return;
  }

  mat_polys.reserve(ncomponents);
  for (int ipoly = 0; ipoly < ncomponents; ipoly++) {
    int nvrts = poly_brep[ipoly].numvertices;
    std::vector<Point3> curpoly_vrts;
    curpoly_vrts.reserve(nvrts);
    // We only store unique vertices, so we need a map from r3d node indices
    // to MatPoly node indices for when we process faces
    std::vector<int> r3d2matpoly_vrt_ids(nvrts, -1);
    for (int ivrt = 0; ivrt < nvrts; ivrt++) {
      Point3 cur_vrt;
      for (int ixyz = 0; ixyz < 3; ixyz++)
        cur_vrt[ixyz] = poly_brep[ipoly].vertices[ivrt].xyz[ixyz];
      // Check if this point is already stored
      for (int i = 0; i < curpoly_vrts.size(); i++)
        //Points are considered equivalent if all their respective coordinates 
        //are within machine epsilon from each other
        if (cur_vrt == curpoly_vrts[i]) {
          r3d2matpoly_vrt_ids[ivrt] = i;
          break;
        }
      // If the point is unique, we store it
      if (r3d2matpoly_vrt_ids[ivrt] == -1) {
        r3d2matpoly_vrt_ids[ivrt] = curpoly_vrts.size();
        curpoly_vrts.push_back(cur_vrt);
      }
    }
    curpoly_vrts.shrink_to_fit();

    if (curpoly_vrts.size() < 4) continue;

    int nfaces = poly_brep[ipoly].numfaces;
    std::vector< std::vector<int> > curpoly_faces(nfaces);
    for (int iface = 0; iface < nfaces; iface++) {
      int face_nverts = poly_brep[ipoly].numvertsperface[iface];
      for (int ifv = 0; ifv < face_nverts; ifv++) {
        int cur_vrt_id = r3d2matpoly_vrt_ids[poly_brep[ipoly].faceinds[iface][ifv]];
        // We only add unique node indices to the list of face's nodes
        if (std::find(curpoly_faces[iface].begin(), curpoly_faces[iface].end(), 
                      cur_vrt_id) == curpoly_faces[iface].end())
          curpoly_faces[iface].push_back(cur_vrt_id);
      }
    }

    std::vector<int> hanging_nodes;
    bool removed_nodes = true;
    while (removed_nodes) {
      //Removing nodes can result in new degenerate faces, so keep removing
      //such faces until no hanging nodes are left
      removed_nodes = false;
      // Filter out degenerate faces
      int ind_face = 0;
      while (ind_face < curpoly_faces.size())
        if (curpoly_faces[ind_face].size() > 2) ind_face++;
        else {
          //A degenerate face can collapse into a hanging point on an edge, i.e.
          //a node that is connected to less than three other nodes.
          //We eliminate such nodes from faces: converting a polyhedron that has
          //hanging nodes using r3d_init_poly can result in an r3d_poly
          //with no vertices
          for (int ifn = 0; ifn < curpoly_faces[ind_face].size(); ifn++) {
            int inode = curpoly_faces[ind_face][ifn];
            std::vector<int> iconnected_nodes;
            std::vector<std::pair<int,int>> node_in_faces;
            for (int iface = 0; iface < curpoly_faces.size(); iface++) {
              int in_face_id = std::distance(curpoly_faces[iface].begin(),
                                             std::find(curpoly_faces[iface].begin(),
                                                       curpoly_faces[iface].end(), inode));
              int nface_vrts = curpoly_faces[iface].size();
              if (in_face_id != nface_vrts) {
                //Node of a degenerate face is also a node of the current face
                node_in_faces.push_back(std::pair<int,int>(iface, in_face_id));
                //We look for unique adjacent nodes in the current face
                for (int p0n1 = 0; p0n1 < 2; p0n1++) {
                  int adj_node_face_id = (nface_vrts + in_face_id + 2*p0n1 - 1)%nface_vrts;
                  if (std::find(iconnected_nodes.begin(), iconnected_nodes.end(),
                                curpoly_faces[iface][adj_node_face_id]) == iconnected_nodes.end())
                    iconnected_nodes.push_back(curpoly_faces[iface][adj_node_face_id]);
                }
              }
            }  
            //Check if the current face node is a hanging node, if so, remove it
            if (iconnected_nodes.size() < 3) {
              removed_nodes = true;
              hanging_nodes.push_back(inode);
              for (int icf = 0; icf < node_in_faces.size(); icf++) {
                int iface = node_in_faces[icf].first;
                curpoly_faces[iface].erase(curpoly_faces[iface].begin() + 
                                           node_in_faces[icf].second);
              }
            }
          }
          //Now, remove the degenerate face
          curpoly_faces.erase(curpoly_faces.begin() + ind_face);  
        }
    }

    //We remove nodes no longer present in the poly: passing to r3d_init_poly
    //vertices that are not used by any faces can result in an empty r3d_poly
    if (!hanging_nodes.empty()) {
      nvrts = static_cast<int>(curpoly_vrts.size());
      std::vector<int> old2new_vrt_ids(nvrts, -1);
      std::vector<Point3> shrunk_vrts;
      shrunk_vrts.reserve(nvrts - hanging_nodes.size());
      for (int ivrt = 0; ivrt < nvrts; ivrt++) {
        // Check if this point was removed
        if (std::find(hanging_nodes.begin(), hanging_nodes.end(), ivrt) ==
            hanging_nodes.end()) {
            old2new_vrt_ids[ivrt] = shrunk_vrts.size();
            shrunk_vrts.push_back(curpoly_vrts[ivrt]);
        }
      }
      curpoly_vrts = shrunk_vrts;  

      nfaces = static_cast<int>(curpoly_faces.size());
      std::vector< std::vector<int> > updated_faces(nfaces);
      for (int iface = 0; iface < nfaces; iface++) {
        int face_nverts = static_cast<int>(curpoly_faces[iface].size());
        updated_faces[iface].resize(face_nverts);
        for (int ifv = 0; ifv < face_nverts; ifv++) 
          updated_faces[iface][ifv] = old2new_vrt_ids[curpoly_faces[iface][ifv]];
      }
      curpoly_faces = updated_faces;
    }

    // We do not store polyhedra with less than four faces, 
    // as they are clearly degenerate
    if (curpoly_faces.size() > 3) {
      int inew_poly = static_cast<int>(mat_polys.size());
      mat_polys.push_back(MatPoly<3>());
      mat_polys[inew_poly].initialize(curpoly_vrts, curpoly_faces);
    }
  }
  mat_polys.shrink_to_fit();

  r3d_free_brep(&poly_brep, ncomponents);
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
void
split_convex_matpoly_r3d(const MatPoly<3>& mat_poly,
                         const Plane_t<3>& cutting_plane,
                         MatPoly<3>& lower_halfspace_poly,
                         MatPoly<3>& upper_halfspace_poly,
                         std::vector<double>& lower_halfspace_moments,
                         std::vector<double>& upper_halfspace_moments) {
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
  for (int isp = 0; isp < 2; isp++) {
    //Check if the subpoly is empty
    if (r3d_subpolys[isp].nverts == 0) {
      subpoly_ptrs[isp]->clear();
      subpoly_moments_ptrs[isp]->clear();
      continue;
    }

    //Find the moments for a subpoly
    r3d_reduce(&r3d_subpolys[isp], r3d_moments, POLY_ORDER);

    if (r3d_moments[0] <= std::numeric_limits<double>::epsilon()) {
      subpoly_ptrs[isp]->clear();
      subpoly_moments_ptrs[isp]->clear();
      continue;
    }    

    subpoly_moments_ptrs[isp]->assign(r3d_moments, r3d_moments + 4);

    //Get a MatPoly for a subpoly
    std::vector< MatPoly<3> > sub_matpoly;
    r3dpoly_to_matpolys(r3d_subpolys[isp], sub_matpoly);
    int ncomponents = static_cast<int>(sub_matpoly.size());
    if (ncomponents > 1) {
      // Filter out degenerate components
      int ind_subpoly = 0;
      while (ind_subpoly < sub_matpoly.size())
        if (sub_matpoly[ind_subpoly].moments()[0] > 
            std::numeric_limits<double>::epsilon()) ind_subpoly++;
        else sub_matpoly.erase(sub_matpoly.begin() + ind_subpoly);  
    }
      
    if (ncomponents == 1) 
      *subpoly_ptrs[isp] = sub_matpoly[0];
    else
      throw std::runtime_error("Non-convex MatPoly is split using the method for convex MatPoly's!");
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
           const bool all_convex = false) : 
           matpolys_(matpolys), 
           cutting_plane_(cutting_plane),
           all_convex_(all_convex) {}

  /*! 
    @brief Splits a MatPoly into two sets of convex MatPoly's
    @return Structure containing two sets of MatPoly's, one for each half-space,
    with respective aggregated moments.
  */
  HalfSpaceSets_t<3> operator() () const {
    HalfSpaceSets_t<3> hs_sets;

    std::vector< MatPoly<3> >* hs_subpolys_ptrs[2] = {
      &hs_sets.lower_halfspace_set.matpolys, &hs_sets.upper_halfspace_set.matpolys};
    std::vector<double>* hs_moments_ptrs[2] = {
      &hs_sets.lower_halfspace_set.moments, &hs_sets.upper_halfspace_set.moments}; 

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

    for (int ihs = 0; ihs < 2; ihs++)
      hs_moments_ptrs[ihs]->assign(4, 0.0);

    int hs_poly_count[2] = {0, 0};
    int npolys = static_cast<int>(convex_polys->size());

    for (int icp = 0; icp < npolys; icp++) {
      MatPoly<3> cur_subpolys[2];
      std::vector<double> cur_moments[2];
      split_convex_matpoly_r3d((*convex_polys)[icp], cutting_plane_,
                               cur_subpolys[0], cur_subpolys[1],
                               cur_moments[0], cur_moments[1]);

      for (int ihs = 0; ihs < 2; ihs++)
        if (!cur_moments[ihs].empty()) {
          hs_subpolys_ptrs[ihs]->emplace_back(cur_subpolys[ihs]);
          for (int im = 0; im < 4; im++)
            (*hs_moments_ptrs[ihs])[im] += cur_moments[ihs][im];
          hs_poly_count[ihs]++;
        }
    }

    for (int ihs = 0; ihs < 2; ihs++)
      if (hs_poly_count[ihs] == 0)
        hs_moments_ptrs[ihs]->clear();

    return hs_sets;
  }

  //! Default constructor (disabled)
  SplitR3D() = delete;

  //! Assignment operator (disabled)
  SplitR3D& operator = (const SplitR3D&) = delete;

 private:
  const std::vector< MatPoly<3> >& matpolys_;
  const Plane_t<3>& cutting_plane_;
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
  ClipR3D() {}

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
                    const bool planar_faces = false) {
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
    if (!r3d_polys_.empty()) {
      agg_moments.assign(4, 0.0);

      const int POLY_ORDER = 1;
      r3d_real r3d_moments[R3D_NUM_MOMENTS(POLY_ORDER)];
      for(int ipoly = 0; ipoly < r3d_polys_.size(); ipoly++) {
        if (r3d_polys_[ipoly].nverts != 0) {
          r3d_reduce(&r3d_polys_[ipoly], r3d_moments, POLY_ORDER);
          if (r3d_moments[0] > std::numeric_limits<double>::epsilon())
            for (int im = 0; im < 4; im++)
              agg_moments[im] += r3d_moments[im];
        }
      }
    }

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

    int poly_counter = 0;
    const int POLY_ORDER = 1;
    r3d_real r3d_moments[R3D_NUM_MOMENTS(POLY_ORDER)];

    for (int ipoly = 0; ipoly < r3d_polys_.size(); ipoly++) {
      if (r3d_polys_[ipoly].nverts != 0) {
        //r3d does in-place clipping and does not have the const modifier
        //for the plane: we need to make copies
        r3d_poly clipped_poly = r3d_polys_[ipoly];
        r3d_plane plane_copy = r3d_cut_plane_;
        r3d_clip(&clipped_poly, &plane_copy, 1);
        if (clipped_poly.nverts != 0) {
          //Find the moments for the part in the lower halfspace
          r3d_reduce(&clipped_poly, r3d_moments, POLY_ORDER);
          if (r3d_moments[0] > std::numeric_limits<double>::epsilon()) {
            for (int im = 0; im < 4; im++)
              below_plane_moments[im] += r3d_moments[im];
            poly_counter++;
          }
        }
      }
    }
    
    if (poly_counter == 0)
      below_plane_moments.clear();
      
    return below_plane_moments;
  }

  //! Assignment operator (disabled)
  ClipR3D& operator = (const ClipR3D&) = delete;

 private:
  std::vector<r3d_poly> r3d_polys_;
  r3d_plane r3d_cut_plane_;
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
void get_intersection_moments(const MatPoly<3>& mat_poly,
                              const r3d_poly& r3dpoly,
                              std::vector<double>& intersection_moments,
                              bool convex_matpoly = false) {   
  const int POLY_ORDER = 1;
  r3d_real r3d_moments[R3D_NUM_MOMENTS(POLY_ORDER)];

  if (convex_matpoly) {                                
    std::vector< Plane_t<3> > face_planes;
    mat_poly.face_planes(face_planes);
    r3d_poly intersection = r3dpoly;

    int nplanes = static_cast<int>(face_planes.size());
    r3d_plane* r3d_face_planes = new r3d_plane [nplanes];
    for (int iplane = 0; iplane < nplanes; iplane++) {
      for (int ixyz = 0; ixyz < 3; ixyz++)
        r3d_face_planes[iplane].n.xyz[ixyz] = -face_planes[iplane].normal[ixyz];
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
        for (int ixyz = 0; ixyz < 3; ixyz++)
          r3d_face_planes[iplane].n.xyz[ixyz] = -face_planes[iplane].normal[ixyz];
        r3d_face_planes[iplane].d = -face_planes[iplane].dist2origin;
      }

      r3d_clip(&intersection, r3d_face_planes, 4);

      r3d_reduce(&intersection, r3d_moments, POLY_ORDER);
      for (int im = 0; im < 4; im++)
        intersection_moments[im] += r3d_moments[im];
    }
  }
}

BoundingBox_t<3> r3d_poly_bounding_box(const r3d_poly& r3dpoly) {
  BoundingBox_t<3> bbox;
  for (int ivrt = 0; ivrt < r3dpoly.nverts; ivrt++)
    for (int idim = 0; idim < 3; idim++) {
      if (r3dpoly.verts[ivrt].pos.xyz[idim] < bbox.min[idim])
        bbox.min[idim] = r3dpoly.verts[ivrt].pos.xyz[idim];

      if (r3dpoly.verts[ivrt].pos.xyz[idim] > bbox.max[idim])
        bbox.max[idim] = r3dpoly.verts[ivrt].pos.xyz[idim];
    }
  return bbox;
}

}  // namespace Tangram

#endif  // TANGRAM_INTERSECT_SPLIT_R3D_H_
