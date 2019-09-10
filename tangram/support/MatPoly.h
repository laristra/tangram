/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#ifndef TANGRAM_MATPOLY_H_
#define TANGRAM_MATPOLY_H_

#include <vector>
#include <array>
#include <algorithm>
#include <numeric>

#include "tangram/support/tangram.h"

namespace Tangram {

double polygon3d_area(const std::vector<Point3>& points,
                      const std::vector<int>& poly_vertices,
                      double dst_tol) {
  double poly_area = 0.0;

  int nvrts = static_cast<int>(poly_vertices.size());
  int nnz_edges = 0;
  if (nvrts == 3) {
    Vector<3> side_vec[2];
    for (int iside = 0; iside < 2; iside++) {
      side_vec[iside] = points[poly_vertices[iside + 1]] - points[poly_vertices[0]];
      if (side_vec[iside].norm() >= dst_tol) nnz_edges++;
    }
    if ((nnz_edges == 2) && 
        ((points[poly_vertices[2]] - points[poly_vertices[1]]).norm() >= dst_tol))
      poly_area = 0.5*cross(side_vec[0], side_vec[1]).norm();
      
    return poly_area;
  }

  Point3 gcenter;
  for (int ivrt = 0; ivrt < nvrts; ivrt++)
    gcenter += points[poly_vertices[ivrt]];
  gcenter /= nvrts;
  
  for (int ivrt = 0; ivrt < nvrts; ivrt++) {
    int ifv = poly_vertices[ivrt];
    int isv = poly_vertices[(ivrt + 1)%nvrts];
    Vector<3> side_vec = points[isv] - points[ifv];
    if (side_vec.norm() >= dst_tol) nnz_edges++;

    double tri_size = 
      0.5*cross(side_vec, gcenter - points[ifv]).norm();

    poly_area += tri_size;
  }

  if (nnz_edges < 3)
    poly_area = 0.0;

  return poly_area;
}

void polygon3d_moments(const std::vector<Point3>& points,
                       const std::vector<int>& poly_vertices,
                       std::vector<double>& poly_moments,
                       double dst_tol) {
  poly_moments.assign(4, 0.0);

  int nvrts = static_cast<int>(poly_vertices.size());
  int nnz_edges = 0;
  if (nvrts == 3) {
    Vector<3> side_vec[2];
    for (int iside = 0; iside < 2; iside++) {
      side_vec[iside] = points[poly_vertices[iside + 1]] - points[poly_vertices[0]];
      if (side_vec[iside].norm() >= dst_tol) nnz_edges++;
    }

    if ((nnz_edges == 2) && 
        ((points[poly_vertices[2]] - points[poly_vertices[1]]).norm() >= dst_tol)) {
      poly_moments[0] = 0.5*cross(side_vec[0], side_vec[1]).norm();

      for (int ivrt = 0; ivrt < 3; ivrt++)
        for (int idim = 0; idim < 3; idim++)
          poly_moments[idim + 1] += points[poly_vertices[ivrt]][idim];
      for (int idim = 0; idim < 3; idim++)    
        poly_moments[idim + 1] *= poly_moments[0]/3.0;
    }

    return;
  }

  Point3 gcenter;
  for (int ivrt = 0; ivrt < nvrts; ivrt++)
    gcenter += points[poly_vertices[ivrt]];
  gcenter /= nvrts;
  
  for (int ivrt = 0; ivrt < nvrts; ivrt++) {
    int ifv = poly_vertices[ivrt];
    int isv = poly_vertices[(ivrt + 1)%nvrts];
    Vector<3> side_vec = points[isv] - points[ifv];
    if (side_vec.norm() >= dst_tol) nnz_edges++;

    double tri_size = 
      0.5*cross(side_vec, gcenter - points[ifv]).norm();

    poly_moments[0] += tri_size;
    for (int idim = 0; idim < 3; idim++)
      poly_moments[idim + 1] += 
        tri_size*(points[ifv][idim] + points[isv][idim] + gcenter[idim])/3.0;
  }

  if (nnz_edges < 3)
    poly_moments.assign(4, 0.0);
}

Vector<3> polygon3d_normal(const std::vector<Point3>& points,
                           const std::vector<int>& poly_vertices,
                           double dst_tol) {
  Vector3 poly_normal(0.0, 0.0, 0.0);

  int nvrts = static_cast<int>(poly_vertices.size());
  if (nvrts == 3) {
    int nnz_edges = 0;
    Vector<3> side_vec[2];
    for (int iside = 0; iside < 2; iside++) {
      side_vec[iside] = points[poly_vertices[iside + 1]] - points[poly_vertices[0]];
      if (side_vec[iside].norm() >= dst_tol) nnz_edges++;
    }

    if ((nnz_edges == 2) && 
        ((points[poly_vertices[2]] - points[poly_vertices[1]]).norm() >= dst_tol)) {
      poly_normal = cross(side_vec[0], side_vec[1]);
      poly_normal.normalize();
    }
    
    return poly_normal;
  }

  std::vector<double> poly_moments;
  polygon3d_moments(points, poly_vertices, poly_moments, dst_tol);

  if (poly_moments[0] == 0.0)
    return poly_normal;

  Point3 poly_centroid;
  for (int idim = 0; idim < 3; idim++)
    poly_centroid[idim] = poly_moments[idim + 1]/poly_moments[0];

  for (int ivrt = 0; ivrt < nvrts; ivrt++) {
    Vector<3> side_vec = points[poly_vertices[(ivrt + 1)%nvrts]] - 
                         points[poly_vertices[ivrt]];
    Vector<3> tri_normal = cross(side_vec, poly_centroid - points[poly_vertices[ivrt]]);
    poly_normal += tri_normal;
  }
  poly_normal.normalize();

  return poly_normal;  
}

template <int D>
class MatPoly {
 public:
  /*!
    @brief Constructor, undefined material ID corresponds to -1 value
    @param material_id  ID of the material this poly contains
  */
  MatPoly(int const material_id = -1) : material_id_(material_id) {}

  /*! Destructor */
  ~MatPoly() {}

  /*!
   @brief Assignment operator
   @param source_poly  Material poly to copy from
  */
  MatPoly& operator=(const MatPoly& source_poly) {
    material_id_ = source_poly.material_id_;
    vertex_points_ = source_poly.vertex_points_;
    face_vertices_ = source_poly.face_vertices_;
    dst_tol_ = source_poly.dst_tol_;
    nvertices_ = source_poly.nvertices_;
    nfaces_ = source_poly.nfaces_;
    moments_ = source_poly.moments_;
  }
  
  /*!
    @brief ID of the material this poly contains
    @return Material ID
  */
  int mat_id() const { return material_id_; }

  /*!
    @brief Set the ID of the material contained in this poly
    @param material_id  Material ID to assign to this poly: should be a valid (>=0) value
  */
  void set_mat_id(int const mat_id) {
#ifdef DEBUG
    assert(mat_id >= 0);
#endif
    material_id_ = mat_id;
  }
  
  /*!
   @brief Set the material ID for this poly to the undefined state
  */
  void reset_mat_id() { material_id_ = -1; }

  /*!
   @brief Resets the MatPoly data
  */
  void clear() {
    material_id_ = -1;
    vertex_points_.clear();
    face_vertices_.clear();
    nvertices_ = 0;
    nfaces_ = 0;
    moments_.clear();    
  }

  /*!
   @brief Initialize a 2D polygon from its vertices
   @param poly_points  Vector of coordinates of material polygon's vertices
   listed counter-clockwise
   @param dst_tol  Distance tolerance
  */
  void initialize(const std::vector<Point2>& poly_points,
                  const double dst_tol);
  /*!
   @brief Initialize a 3D polyhedron
   @param vertex_points  Vector of coordinates of material polyhedron's vertices
   @param face_vertices  Faces of the material polyhedron, every face is given by
   IDs of its vertices in counter-clockwise order
   @param dst_tol  Distance tolerance   
  */
  void initialize(const std::vector<Point3>& vertex_points,
                  const std::vector< std::vector<int> >& face_vertices,
                  const double dst_tol);
  
  /*!
    @brief Distance tolerance used by this poly
    @return Distance tolerance
  */
  double dst_tol() const { return dst_tol_; }

  /*!
   @brief Points for all the vertices of the material poly
   @return  Vector of points of matpoly
  */
  const std::vector< Point<D> >& points() const { return vertex_points_; }

  /*!
   @brief Indices of vertices of the material poly's face
   @param face_id  ID of the face of the material poly
   @return  Vector of indices of face's vertices
  */
  const std::vector<int>& face_vertices(int const face_id) const {
#ifdef DEBUG
    assert((face_id >= 0) && (face_id < nfaces_));
#endif
    return face_vertices_[face_id];
  }
  
  /*!
   @brief Coordinates of vertices of the material poly's face
   @param face_id  ID of the face of the material poly
   @return  Vector of coordinates of face's vertices
  */
  std::vector< Point<D> > face_points(int const face_id) const {
#ifdef DEBUG
    assert((face_id >= 0) && (face_id < nfaces_));
#endif
    int nvrts = static_cast<int>(face_vertices_[face_id].size());
    std::vector< Point<D> > fpoints;
    fpoints.reserve(nvrts);
    for (int ivrt = 0; ivrt < nvrts; ivrt++)
      fpoints.push_back(vertex_points_[face_vertices_[face_id][ivrt]]);

    return fpoints;
  }

  /*!
   @brief All indices of vertices of the material poly's faces
   @return  Vector of indices of face's vertices
  */
  const std::vector<std::vector<int>>& face_vertices() const {
    return face_vertices_;
  }
  
  /*!
   @brief Coordinates of the material poly's vertex
   @param vertex_id  ID of the vertex of the material poly
   @return  Coordinates of that vertex
  */
  Point<D> vertex_point(int const vertex_id) const {
#ifdef DEBUG
    assert((vertex_id >= 0) && (vertex_id < nvertices_));
#endif
    return vertex_points_[vertex_id];
  }
  
  /*!
   @brief Coordinates of the centroid of the material poly's face
   @param face_id  ID of the face of the material poly
   @return  Coordinates of the centroid of that face
  */
  Point<D> face_centroid(int const face_id) const;
  
  /*!
   @brief Facetization of the 2D polygon's boundary: simply creates a copy of the polygon
   @param faceted_poly  Pointer to a 2D material polygon object to which
   this polygon is copied to
  */
  void faceted_matpoly(MatPoly<2>* faceted_poly) const;
  /*!
   @brief Facetization of the 3D polyhedron's boundary: uses face centroids to perform
   a triangulation of non-triangular faces
   @param faceted_poly  Pointer to a 3D material polyhedron object to which
   the faceted version of this polyhedron is written to
  */
  void faceted_matpoly(MatPoly<3>* faceted_poly) const;
  
  /*!
   @brief Number of material poly's vertices
   @return  Number of vertices
  */
  int num_vertices() const { return nvertices_; }
  /*!
   @brief Number of material poly's faces
   @return  Number of faces
  */
  int num_faces() const { return nfaces_; }

  /*!
   @brief Moments of material poly, when the method is called for the first time,
   moments will be computed and stored for future use
   @return  Vector of moments; moments[0] is the size, moments[i+1]/moments[0] is i-th
   coordinate of the centroid
  */  
  const std::vector<double>& moments() const {
    if (moments_.empty())
      compute_moments(moments_);
    return moments_;
  }
  
  /*!
   @brief Assigns externally computed moments to this material poly.
   Can also be used to copy moments from the faceted poly to the original one
   (for example, post-decomposition)
   @param moments Externally computed moments, moments[0] is the size, 
   moments[i+1]/moments[0] is i-th coordinate of the centroid
  */  
  void assign_moments(const std::vector<double>& moments) const {
#ifdef DEBUG
    assert(moments.size() == D + 1);
#endif    
    moments_ = moments;
  }

  /*!
    @brief Decomposes this MatPoly into MatPoly's using its centroid.
    If faces of MatPoly are planar, MatPoly's in the decomposition will be convex.
    @param[in] mat_poly MatPoly to decompose
    @param[out] convex_matpolys Vector of MatPoly's: 
    as many MatPoly's as mat_poly has faces will be appended to it.
  */
  void decompose(std::vector< MatPoly<D> >& sub_polys) const;

  /*!
    @brief Facetizes and decomposes this MatPoly into simplex MatPoly's 
    using its centroid.
    @param[in] mat_poly MatPoly to decompose
    @param[out] convex_matpolys Vector of MatPoly's: 
    as many MatPoly's as mat_poly has facets will be appended to it.
  */
  void facetize_decompose(std::vector< MatPoly<D> >& sub_polys) const;

  /*!
    @brief For every face, returns a plane containing that face.
    Important: faces should be planar, so you might need to facetize
    the MatPoly first. In 2D, faces of size below the distance tolerance 
    will be omitted. In 3D, faces should have at least three edges of size
    above the distance tolerance. If the number of valid planes is less 
    than four, empty vector will be returned.
    @param fplanes  Vector of planes for all the faces of this MatPoly
  */
  void face_planes(std::vector< Plane_t<D> >& fplanes) const;

  BoundingBox_t<D> bounding_box() const {
    BoundingBox_t<D> bbox;
    for (int ivrt = 0; ivrt < nvertices_; ivrt++)
      for (int idim = 0; idim < D; idim++) {
        if (vertex_points_[ivrt][idim] < bbox.min[idim])
          bbox.min[idim] = vertex_points_[ivrt][idim];

        if (vertex_points_[ivrt][idim] > bbox.max[idim])
          bbox.max[idim] = vertex_points_[ivrt][idim];        
      }
    return bbox;
  }

 protected:
  /*!
   @brief Computes moments of this material poly
   @param moments Computed moments, moments[0] is the size, 
   moments[i+1]/moments[0] is i-th coordinate of the centroid
  */  
  void compute_moments(std::vector<double>& moments) const;
 private:

  int material_id_;  // material ID of this matpoly
  std::vector< Point<D> > vertex_points_;  // coordinates of vertices
  std::vector< std::vector<int> > face_vertices_;  // vertices of faces
  double dst_tol_; // distance tolerance
  
  int nvertices_ = 0;  // number of vertices
  int nfaces_ = 0;  //number of faces
  mutable std::vector<double> moments_; //moments of this matpoly
};  // class MatPoly

/*!
 @brief Initialize a 2D polygon from its vertices
 @param poly_points  Vector of coordinates of material polygon's vertices
                     listed counter-clockwise
*/
template<>
void MatPoly<2>::initialize(const std::vector<Point2>& poly_points, 
                            double dst_tol) {
  dst_tol_ = dst_tol;
  nvertices_ = static_cast<int>(poly_points.size());
#ifdef DEBUG
  assert(nvertices_ > 2);
#endif
  nfaces_ = nvertices_;
 
  vertex_points_ = poly_points;
  face_vertices_.resize(nfaces_);
  for (int iface = 0; iface < nfaces_; iface++)
    face_vertices_[iface] = { iface, (iface + 1)%nfaces_ };
}

/*!
 @brief Initialize a 3D polyhedron
 @param vertex_points  Vector of coordinates of material polyhedron's vertices
 @param face_vertices  Faces of the material polyhedron, every face is given by
                       IDs of its vertices in counter-clockwise order
 */
template<>
void MatPoly<3>::initialize(const std::vector<Point3>& vertex_points,
                            const std::vector< std::vector<int> >& face_vertices,
                            double dst_tol) {
  dst_tol_ = dst_tol;
  nvertices_ = static_cast<int>(vertex_points.size());
  nfaces_ = static_cast<int>(face_vertices.size());
#ifdef DEBUG
  assert(nvertices_ > 3);
  assert(nfaces_ > 3);
#endif
  
  vertex_points_ = vertex_points;
  face_vertices_ = face_vertices;
}

/*!
  @brief Coordinates of the centroid of the material polygon's face
  @param face_id  ID of the face of the material polygon
  @return  Coordinates of the centroid of that face
*/
template<>  
Point2 MatPoly<2>::face_centroid(int const face_id) const {
#ifdef DEBUG
  assert((face_id >= 0) && (face_id < nfaces_));
#endif
  return 0.5*(vertex_points_[face_id] + vertex_points_[(face_id + 1)%nvertices_]);
}

/*!
  @brief Coordinates of the centroid of the material polyhedron's face
  @param face_id  ID of the face of the material polyhedron
  @return  Coordinates of the centroid of that face
*/
template<>  
Point3 MatPoly<3>::face_centroid(int const face_id) const {
#ifdef DEBUG
  assert((face_id >= 0) && (face_id < nfaces_));
#endif
  std::vector<double> face_moments;
  polygon3d_moments(vertex_points_, face_vertices_[face_id], face_moments, dst_tol_);

  if (face_moments[0] != 0.0)
    return Point3(face_moments[1]/face_moments[0], face_moments[2]/face_moments[0],
                  face_moments[3]/face_moments[0]);
  else {
    int nvrts = static_cast<int>(face_vertices_[face_id].size());
    Point3 gcenter;
    for (int ivrt = 0; ivrt < nvrts; ivrt++)
      gcenter += vertex_points_[face_vertices_[face_id][ivrt]];
    gcenter /= nvrts;

    return gcenter;
  }
}

/*!
 @brief Facetization of the 2D polygon's boundary: simply creates a copy of the polygon
 @param faceted_poly  Pointer to a 2D material polygon object to which
                      this polygon is copied to
*/
template<>
void MatPoly<2>::faceted_matpoly(MatPoly<2>* faceted_poly) const {
  *faceted_poly = *this;
}

/*!
 @brief Facetization of the 3D polyhedron's boundary: uses face centroids to perform
        a triangulation of non-triangular faces
 @param faceted_poly  Pointer to a 3D material polyhedron object to which
                      the faceted version of this polyhedron is written to
*/
template<>
void MatPoly<3>::faceted_matpoly(MatPoly<3>* faceted_poly) const {
  if (material_id_ >= 0)
    faceted_poly->set_mat_id(material_id_);
  else
    faceted_poly->reset_mat_id();

  std::vector<Point3> facetedpoly_vertices = vertex_points_;
  std::vector< std::vector<int> > facetedpoly_faces_;
  facetedpoly_faces_.reserve(nfaces_);
  for (int iface = 0; iface < nfaces_; iface++) {
    int nvrts = static_cast<int>(face_vertices_[iface].size());
    if (nvrts == 3) {
      facetedpoly_faces_.push_back(face_vertices_[iface]);
      continue;
    }
    int icenvrt = static_cast<int>(facetedpoly_vertices.size());
    facetedpoly_vertices.push_back(face_centroid(iface));
    for (int ivrt = 0; ivrt < nvrts; ivrt++)
      facetedpoly_faces_.push_back({ icenvrt, face_vertices_[iface][ivrt],
                                     face_vertices_[iface][(ivrt + 1)%nvrts] });
  }
  
  faceted_poly->initialize(facetedpoly_vertices, facetedpoly_faces_, dst_tol_);
  if (!moments_.empty())
    faceted_poly->assign_moments(moments_);
}

/*!
  @brief Computes moments of this material polygon,
  @param moments Computed moments: moments[0] is area, 
  moments[i+1]/moments[0] is i-th coordinate of the centroid, i=1,2
*/ 
template<>
void MatPoly<2>::compute_moments(std::vector<double>& moments) const {
  moments.assign(3, 0.0);

  for (int ivrt = 0; ivrt < nvertices_; ivrt++) {
    double cur_term = vertex_points_[ivrt][0]*vertex_points_[(ivrt + 1)%nvertices_][1] - 
                      vertex_points_[ivrt][1]*vertex_points_[(ivrt + 1)%nvertices_][0];
    moments[0] += cur_term;
    for (int idim = 0; idim < 2; idim++)
      moments[idim + 1] += cur_term*(
        vertex_points_[ivrt][idim] + vertex_points_[(ivrt + 1)%nvertices_][idim]);
  }
  moments[0] /= 2.0;  
  for (int idim = 0; idim < 2; idim++)
    moments[idim + 1] /= 6.0;
}

/*!
  @brief Computes moments of this material polyhedron
  @param moments Computed moments, moments[0] is volume, 
  moments[i+1]/moments[0] is i-th coordinate of the centroid, i=1,2,3
*/  
template<>
void MatPoly<3>::compute_moments(std::vector<double>& moments) const {
  moments.assign(4, 0.0); 

  for (int iface = 0; iface < nfaces_; iface++) {
    std::vector<double> face_moments;
    polygon3d_moments(vertex_points_, face_vertices_[iface], face_moments, dst_tol_);

    if (face_moments[0] == 0.0)
      continue;

    std::vector<Point3> face_pts = face_points(iface);
    std::vector< std::vector<int> > itri_pts;
    
    int nvrts = face_vertices_[iface].size();
    if (nvrts == 3)
      itri_pts.push_back({0, 1, 2});
    else {
      itri_pts.reserve(nvrts);
      face_pts.push_back(Point3(face_moments[1]/face_moments[0],
        face_moments[2]/face_moments[0], face_moments[3]/face_moments[0]));
      for (int ivrt = 0; ivrt < nvrts; ivrt++)
        itri_pts.push_back({nvrts, ivrt, (ivrt + 1)%nvrts});
    }

    for (int itri = 0; itri < itri_pts.size(); itri++) {
      Vector3 vcp = cross(face_pts[itri_pts[itri][1]] - face_pts[itri_pts[itri][0]], 
                          face_pts[itri_pts[itri][2]] - face_pts[itri_pts[itri][0]]);

      moments[0] += dot(vcp, face_pts[itri_pts[itri][0]].asV());
      for (int idim = 0; idim < 3; idim++)
        for (int ivrt = 0; ivrt < 3; ivrt++)
          moments[idim + 1] += vcp[idim]*pow(face_pts[itri_pts[itri][ivrt]][idim] + 
                                             face_pts[itri_pts[itri][(ivrt + 1)%3]][idim], 2);
    }
  }

  moments[0] /= 6.0;
  for (int idim = 0; idim < 3; idim++)
    moments[idim + 1] /= 48.0;
}

/*!
  @brief Decomposes a 2D MatPoly into triangular MatPoly's using its centroid.
  @param[in] mat_poly MatPoly to decompose
  @param[out] convex_matpolys Vector of MatPoly's: 
  as many MatPoly's as mat_poly has faces will be appended to it.
*/
template <>
void MatPoly<2>::decompose(std::vector< MatPoly<2> >& sub_polys) const {
  std::vector<double> moments;
  if (moments_.empty()) 
    compute_moments(moments_);

  Point2 matpoly_cen;
  for (int ixy = 0; ixy < 2; ixy++)
    matpoly_cen[ixy] = moments_[ixy + 1]/moments_[0];

  int offset = static_cast<int>(sub_polys.size());
  sub_polys.resize(offset + nfaces_);

  for (int iface = 0; iface < nfaces_; iface++) {
    std::vector<Point2> subpoly_points = face_points(iface);
    subpoly_points.push_back(matpoly_cen);
    sub_polys[offset + iface].initialize(subpoly_points, dst_tol_);
  }
}

/*!
  @brief Decomposes a 3D MatPoly into MatPoly's using its centroid.
  If faces of MatPoly are planar, MatPoly's in the decomposition will be convex.
  @param[in] mat_poly MatPoly to decompose
  @param[out] convex_matpolys Vector of MatPoly's: 
  as many MatPoly's as mat_poly has faces will be appended to it.
*/
template <>
void MatPoly<3>::decompose(std::vector< MatPoly<3> >& sub_polys) const {
  if (moments_.empty()) 
    compute_moments(moments_);

  Point3 matpoly_cen;
  for (int ixyz = 0; ixyz < 3; ixyz++)
    matpoly_cen[ixyz] = moments_[ixyz + 1]/moments_[0];

  int offset = static_cast<int>(sub_polys.size());
  sub_polys.resize(offset + nfaces_);

  for (int iface = 0; iface < nfaces_; iface++) {
    int face_nvrts = static_cast<int>(face_vertices_[iface].size());

    std::vector<Point3> subpoly_vrts(face_nvrts + 1);
    std::vector< std::vector<int> > subpoly_faces(face_nvrts + 1);
    for (int ivrt = 0; ivrt < face_nvrts; ivrt++) {
      subpoly_vrts[ivrt] = vertex_points_[face_vertices_[iface][ivrt]];
      subpoly_faces[ivrt] = {face_nvrts, (ivrt + 1)%face_nvrts, ivrt};
    }
    subpoly_vrts[face_nvrts] = matpoly_cen;
    subpoly_faces[face_nvrts].resize(face_nvrts);
    std::iota(subpoly_faces[face_nvrts].begin(), 
              subpoly_faces[face_nvrts].end(), 0);      
    
    sub_polys[offset + iface].initialize(subpoly_vrts, subpoly_faces, dst_tol_);
  }
}

/*!
  @brief Decomposes a 2D MatPoly into triangular MatPoly's using its centroid.
  This method is identical to decompose.
  @param[in] mat_poly MatPoly to decompose
  @param[out] convex_matpolys Vector of MatPoly's: 
  as many MatPoly's as mat_poly has faces will be appended to it.
*/
template <>
void MatPoly<2>::facetize_decompose(std::vector< MatPoly<2> >& sub_polys) const {
  decompose(sub_polys);
}

/*!
  @brief Facetizes and decomposes a 3D MatPoly into tetrahedral MatPoly's 
  using its centroid.
  @param[in] mat_poly MatPoly to decompose
  @param[out] convex_matpolys Vector of MatPoly's: 
  as many MatPoly's as mat_poly has facets will be appended to it.
*/
template <>
void MatPoly<3>::facetize_decompose(std::vector< MatPoly<3> >& sub_polys) const {
  std::vector< std::vector<int> > tet_faces(4);
  for (int ivrt = 0; ivrt < 3; ivrt++)
    tet_faces[ivrt] = {3, (ivrt + 1)%3, ivrt};
  tet_faces[3] = {0, 1, 2};  

  if (moments_.empty()) 
    compute_moments(moments_);

  std::vector<Point3> tet_vrts(4);
  for (int ixyz = 0; ixyz < 3; ixyz++)
    tet_vrts[3][ixyz] = moments_[ixyz + 1]/moments_[0];

  for (int iface = 0; iface < nfaces_; iface++) {
    int face_nvrts = static_cast<int>(face_vertices_[iface].size());
    if (face_nvrts == 3) {
      for (int ivrt = 0; ivrt < 3; ivrt++)
        tet_vrts[ivrt] = vertex_points_[face_vertices_[iface][ivrt]];
      int icur_poly = static_cast<int>(sub_polys.size());
      sub_polys.push_back(MatPoly<3>(material_id_));
      sub_polys[icur_poly].initialize(tet_vrts, tet_faces, dst_tol_);
    }
    else {
      tet_vrts[0] = face_centroid(iface);
      for (int ivrt = 0; ivrt < face_nvrts; ivrt++) {
        tet_vrts[1] = vertex_points_[face_vertices_[iface][ivrt]];
        tet_vrts[2] = vertex_points_[face_vertices_[iface][(ivrt + 1)%face_nvrts]];

        int icur_poly = static_cast<int>(sub_polys.size());
        sub_polys.push_back(MatPoly<3>(material_id_));
        sub_polys[icur_poly].initialize(tet_vrts, tet_faces, dst_tol_);
      }
    }
  }
}

/*!
  @brief Initializes a material polygon that corresponds to (i.e. has
  the same geometry as) a given mesh cell.
  @tparam Mesh_Wrapper A lightweight wrapper to a specific input mesh
  implementation that provides certain functionality

  @param[in] Mesh Mesh wrapper
  @param[in] cellid Index of a cell
  @param[out] mat_poly Pointer to the material polygon to initialize
  @param[in] dst_tol Distance tolerance to be used by that material polygon
*/
template <class Mesh_Wrapper>
void cell_get_matpoly(const Mesh_Wrapper& Mesh,
                      const int cellid,
                      MatPoly<2>* mat_poly,
                      const double dst_tol) {
#ifdef DEBUG                        
  assert(Mesh.space_dimension() == 2);
#endif

  mat_poly->reset_mat_id();
  std::vector<Point2> poly_points;
  Mesh.cell_get_coordinates(cellid, &poly_points);

  mat_poly->initialize(poly_points, dst_tol);
}

/*!
  @brief Initializes a material polyhedron that corresponds to (i.e. has
  the same geometry as) a given mesh cell.
  @tparam Mesh_Wrapper A lightweight wrapper to a specific input mesh
  implementation that provides certain functionality

  @param[in] Mesh Mesh wrapper
  @param[in] cellid Index of a cell
  @param[out] mat_poly Pointer to the material polyhedron to initialize
  @param[in] dst_tol Distance tolerance to be used by that material polyhedron
*/
template <class Mesh_Wrapper>
void cell_get_matpoly(const Mesh_Wrapper& Mesh,
                      const int cellid,
                      MatPoly<3>* mat_poly,
                      const double dst_tol) {
#ifdef DEBUG                        
  assert(Mesh.space_dimension() == 3);
#endif

  mat_poly->reset_mat_id();
  std::vector<Point3> poly_points;
  std::vector< std::vector<int> > poly_faces;

  std::vector<int> cnodes;
  Mesh.cell_get_nodes(cellid, &cnodes);
  int ncnodes = cnodes.size();
  poly_points.resize(ncnodes);
  for (int n = 0; n < ncnodes; ++n)
    Mesh.node_get_coordinates(cnodes[n], &poly_points[n]);

  std::vector<int> cfaces, cfdirs;
  Mesh.cell_get_faces_and_dirs(cellid, &cfaces, &cfdirs);
  int ncfaces = cfaces.size();
  for (int f = 0; f < ncfaces; f++) {
    std::vector<int> fnodes;
    Mesh.face_get_nodes(cfaces[f], &fnodes);
    int nfnodes = fnodes.size();
    
    //Check that the order of nodes is ccw
    if (cfdirs[f] != 1)
      std::reverse(fnodes.begin(), fnodes.end());

    // Get the local indices (in the cell node list) of the face nodes
    std::vector<int> fnodes_local(nfnodes);
    for (int n = 0; n < nfnodes; n++) {
      fnodes_local[n] = std::distance(cnodes.begin(),
        std::find(cnodes.begin(), cnodes.end(), fnodes[n]));

      assert(fnodes_local[n] != ncnodes);
    }
    poly_faces.emplace_back(fnodes_local);
  }

  mat_poly->initialize(poly_points, poly_faces, dst_tol);
}

/*!
  @brief For every face, returns a line containing that face.
  Important: faces with length below the distance tolerance will be
  omitted. If the number of valid lines is less than three, empty
  vector will be returned.
  @param flines  Vector of lines for all the faces of this MatPoly
*/
template <>
void MatPoly<2>::face_planes(std::vector< Plane_t<2> >& flines) const {
  flines.clear();

  for (int iface = 0; iface < nfaces_; iface++) {
    Plane_t<2> face_line;
    face_line.normal[0] = vertex_points_[face_vertices_[iface][1]][1] - 
                          vertex_points_[face_vertices_[iface][0]][1];
    face_line.normal[1] = -(vertex_points_[face_vertices_[iface][1]][0] - 
                            vertex_points_[face_vertices_[iface][0]][0]);

    double normal_len = face_line.normal.norm();
    if (normal_len >= dst_tol_) {
      face_line.normal /= normal_len;
      face_line.dist2origin = 
        -dot(vertex_points_[face_vertices_[iface][0]].asV(), face_line.normal);

      flines.push_back(face_line); 
    }
  }

  if (flines.size() < 3)
    flines.clear();
}

/*!
  @brief For every face, returns a plane containing that face.
  Important: faces should be planar, so you might need to facetize
  the MatPoly first. Faces with less that three edges longer than the distance 
  tolerance will be omitted. If the number of valid planes is less than four,
  empty vector will be returned.
  @param fplanes  Vector of planes for all the faces of this MatPoly
*/
template <>
void MatPoly<3>::face_planes(std::vector< Plane_t<3> >& fplanes) const {
  fplanes.clear();

  for (int iface = 0; iface < nfaces_; iface++) {
    Plane_t<3> face_plane;
    face_plane.normal = polygon3d_normal(vertex_points_, face_vertices_[iface], 
                                         dst_tol_);
    if (face_plane.normal.is_zero(dst_tol_))
      continue;
    
    face_plane.dist2origin = 
      -dot(vertex_points_[face_vertices_[iface][0]].asV(), face_plane.normal);

    fplanes.push_back(face_plane);  
  }

  if (fplanes.size() < 4)
    fplanes.clear();
}

/*!
  @brief Eliminates degenerate faces from a polygon given by a sequence of vertices
  @param[in] poly_points Vertices of the polygon in the counter-clockwise order
  @param[in] dst_tol Distance tolerance
  @param[in] reference_pts Preferred points to snap vertices to
  @return MatPoly with no degeneracies
*/
MatPoly<2> natural_selection(const std::vector< Point<2> >& poly_points,
                             const double dst_tol,
                             const std::vector< Point<2> >* reference_pts = nullptr) {
  MatPoly<2> fit_poly;
  int nvrts = static_cast<int>(poly_points.size());
  std::vector< Point<2> > unique_points;
  unique_points.reserve(nvrts);

  for (int ivrt = 0; ivrt < nvrts; ivrt++) {
    // Check if this point is coincident with the next one,
    // i.e. if the respective edge is degenerate: points 
    // are considered coincident if the distance between them
    // is below the distance tolerance
    if ((poly_points[ivrt] - poly_points[(ivrt + 1)%nvrts]).norm() >= dst_tol) {
      bool reference_match = false;
      if (reference_pts != nullptr) {
        for (int irp = 0; irp < reference_pts->size(); irp++)
          if ((poly_points[ivrt] - (*reference_pts)[irp]).norm() < dst_tol) {
            unique_points.push_back((*reference_pts)[irp]);
            reference_match = true;
            break;  
          }
      }

      if (!reference_match)
        unique_points.push_back(poly_points[ivrt]);
    }
  }
  unique_points.shrink_to_fit();

  // Polygons with less than three vertices are clearly degenerate
  if (unique_points.size() > 2)
    fit_poly.initialize(unique_points, dst_tol);

  return fit_poly;
}

/*!
  @brief Eliminates degeneracies from the boundary representation,
  which includes vertices within distance tolerance from other vertices,
  hanging nodes, and degenerate faces with less than three vertices.
  @param[in] vertex_points Vertices of the polyhedron
  @param[in] face_vertices Faces of the polyhedron, every face is given by
  IDs of its vertices in counter-clockwise order
  @param[in] dst_tol Distance tolerance
  @param[in] reference_pts Preferred points to snap vertices to
  @return MatPoly with no degeneracies
*/
MatPoly<3> natural_selection(const std::vector<Point3>& vertex_points,
                             const std::vector< std::vector<int> >& face_vertices,
                             const double dst_tol,
                             const std::vector< Point<3> >* reference_pts = nullptr) {
  MatPoly<3> fit_poly;
  int nvrts = static_cast<int>(vertex_points.size());
  std::vector<Point3> unique_vrts;
  unique_vrts.reserve(nvrts);
  // We only store unique vertices, so we need a map from the original node indices
  // to unique node indices for when we process faces
  std::vector<int> orig2unique_vrt_ids(nvrts, -1);
  for (int ivrt = 0; ivrt < nvrts; ivrt++) {
    for (int i = 0; i < unique_vrts.size(); i++) {
      // Check if this point is coincident with a stored one: points 
      // are considered coincident if the distance between them is below
      // the distance tolerance
      if ((vertex_points[ivrt] - unique_vrts[i]).norm() < dst_tol) {
        orig2unique_vrt_ids[ivrt] = i;
        break;
      }
    }
    // If the point is unique, we store it
    if (orig2unique_vrt_ids[ivrt] == -1) {
      orig2unique_vrt_ids[ivrt] = static_cast<int>(unique_vrts.size());
      bool reference_match = false;
      if (reference_pts != nullptr) {
        // We check if this point matches any of the reference ones
        for (int irp = 0; irp < reference_pts->size(); irp++) {
          const Point<3>& cur_ref_pt = (*reference_pts)[irp];
          if ((vertex_points[ivrt] - cur_ref_pt).norm() < dst_tol) {            
            unique_vrts.push_back(cur_ref_pt);            
            reference_match = true;
            break; 
          }
        }
      }
      if (!reference_match)
        unique_vrts.push_back(vertex_points[ivrt]);
    }
  }
  unique_vrts.shrink_to_fit();
  nvrts = static_cast<int>(unique_vrts.size());

  // Check if it still has a chance to be at least a tetrahedron
  if (nvrts < 4) return fit_poly;

  //First pass on faces: have them reference only unique vertices
  int nfaces = static_cast<int>(face_vertices.size());
  std::vector< std::vector<int> > face_unique_vrts(nfaces);
  for (int iface = 0; iface < nfaces; iface++) {
    int face_nvrts = static_cast<int>(face_vertices[iface].size());
    for (int ifv = 0; ifv < face_nvrts; ifv++) {
      int cur_vrt_id = orig2unique_vrt_ids[face_vertices[iface][ifv]];
      // There can be only one (c) such index in the list of face indices
      if (std::find(face_unique_vrts[iface].begin(), face_unique_vrts[iface].end(), 
                    cur_vrt_id) == face_unique_vrts[iface].end())
        face_unique_vrts[iface].push_back(cur_vrt_id);
    }
  }

  // Hanging vertices can be a problem, so we first find how many unique
  // connections is in each vertice's network
  std::vector< std::vector<int> > connected_vrt_ids(nvrts);
  std::vector< std::vector< std::pair<int, int> > > vrt_in_faces(nvrts);
  for (int ivrt = 0; ivrt < nvrts; ivrt++) {
    for (int iface = 0; iface < face_unique_vrts.size(); iface++) {
      const std::vector<int>& cur_face_vrts = face_unique_vrts[iface];
      int nface_vrts = static_cast<int>(cur_face_vrts.size());
      int id_in_face = std::distance(cur_face_vrts.begin(),
                                     std::find(cur_face_vrts.begin(),
                                               cur_face_vrts.end(), ivrt));
      
      if (id_in_face != nface_vrts) {
        // Vertex is in the current face
        vrt_in_faces[ivrt].push_back(std::pair<int, int>(iface, id_in_face));
        // We check if it has unknown connections with vertices of the current face
        for (int p0n1 = 0; p0n1 < 2; p0n1++) {
          int adj_vrt_face_id = (nface_vrts + id_in_face + 2*p0n1 - 1)%nface_vrts;
          int adj_vrt_id = face_unique_vrts[iface][adj_vrt_face_id];
          if (std::find(connected_vrt_ids[ivrt].begin(), connected_vrt_ids[ivrt].end(),
                        adj_vrt_id) == connected_vrt_ids[ivrt].end())
            connected_vrt_ids[ivrt].push_back(adj_vrt_id);
        }
      }
    }      
  }

  // Now we dispose of hanger-on vertices that don't have enough connections
  // This also affects the connections of vertices which they networked with,
  // so we need to check if those get dragged down the cliff as well
  bool dropped_connections;
  std::vector<int> hanging_vrts_ids;
  do {
    dropped_connections = false;
    for (int ivrt = 0; ivrt < nvrts; ivrt++) {
      int ncv = static_cast<int>(connected_vrt_ids[ivrt].size());
      if (ncv < 3) {
        hanging_vrts_ids.push_back(ivrt);
        // Others break connections first
        for (int icv = 0; icv < ncv; icv++) {
          int icvrt = connected_vrt_ids[ivrt][icv];
          connected_vrt_ids[icvrt].erase(
            std::find(connected_vrt_ids[icvrt].begin(), 
                      connected_vrt_ids[icvrt].end(), ivrt));
        }
        if (ncv == 1) dropped_connections = true;

        bool third_wheel = false;
        for (int icf = 0; icf < vrt_in_faces[ivrt].size(); icf++) {
          int iface = vrt_in_faces[ivrt][icf].first;
          int hvrt_fid = vrt_in_faces[ivrt][icf].second;
          // Check if the two others can connect by themselves
          if (ncv == 2) {
            int nface_vrts = static_cast<int>(face_unique_vrts[iface].size());
            for (int icv = 0; icv < 2; icv++)
              if ( (face_unique_vrts[iface][(hvrt_fid + 1)%nface_vrts] ==
                    connected_vrt_ids[ivrt][(icv + 1)%2]) &&
                   (face_unique_vrts[iface][(nface_vrts + hvrt_fid - 1)%nface_vrts] ==
                    connected_vrt_ids[ivrt][icv]) ) {
                third_wheel = true;
                for (int i = 0; i < 2; i++)
                  connected_vrt_ids[connected_vrt_ids[ivrt][i]].push_back(
                    connected_vrt_ids[ivrt][(i + 1)%2]);
                // And the hanger-on now has no one to cling to
                connected_vrt_ids[ivrt].clear();
                ncv = 0;
              }
          }
          // Gets evicted from the hosting face
          face_unique_vrts[iface].erase(face_unique_vrts[iface].begin() + 
                                        vrt_in_faces[ivrt][icf].second);
          // Tell the following face vertices about their change of address
          for (int ifv = vrt_in_faces[ivrt][icf].second; 
                   ifv < face_unique_vrts[iface].size(); ifv++) {
            int icur_vrt = face_unique_vrts[iface][ifv];
            for (int i = 0; i < vrt_in_faces[icur_vrt].size(); i++)
              if (vrt_in_faces[icur_vrt][i].first == iface) {
                vrt_in_faces[icur_vrt][i].second--;
                break;
              }
          }

        }      
        if (!third_wheel) dropped_connections = true;
      }
    }
  } while (dropped_connections);

  // Filter out degenerate faces
  int face_id = 0;
  while (face_id < face_unique_vrts.size()) {
    // Check if this face is still at least a triangle
    if (face_unique_vrts[face_id].size() > 2) face_id++;
    else {
      //Dispose of the degenerate face
      face_unique_vrts.erase(face_unique_vrts.begin() + face_id);  
    }
  }
  
  // Check if it managed to remain at least a tetrahedron
  nfaces = static_cast<int>(face_unique_vrts.size());
  if (nfaces < 4) return fit_poly;

  // Only the surviving vertices get remembered
  if (!hanging_vrts_ids.empty()) {
    std::vector<int> unique2fit_vrt_ids(nvrts, -1);
    std::vector<Point3> fit_vrts;
    fit_vrts.reserve(nvrts - hanging_vrts_ids.size());
    for (int ivrt = 0; ivrt < nvrts; ivrt++) {
      // Check if this point is still around
      if (std::find(hanging_vrts_ids.begin(), hanging_vrts_ids.end(), ivrt) ==
          hanging_vrts_ids.end()) {
          unique2fit_vrt_ids[ivrt] = fit_vrts.size();
          fit_vrts.push_back(unique_vrts[ivrt]);
      }
    }
    unique_vrts = fit_vrts;  

    std::vector< std::vector<int> > pruned_faces(nfaces);
    for (int iface = 0; iface < nfaces; iface++) {
      int face_nverts = static_cast<int>(face_unique_vrts[iface].size());
      pruned_faces[iface].resize(face_nverts);
      for (int ifv = 0; ifv < face_nverts; ifv++) 
        pruned_faces[iface][ifv] = unique2fit_vrt_ids[face_unique_vrts[iface][ifv]];
    }
    face_unique_vrts = pruned_faces;
  }

  // And they lived happily ever after
  fit_poly.initialize(unique_vrts, face_unique_vrts, dst_tol);
  return fit_poly;
}

/*!
@brief Checks if a given point is interior wrt to a given polyhedron. 
Note: boundary points are no considered to be interior.
@param[in] mat_poly A given polyhedron
@param[in] pt A given point
@param[in] convex_poly Flag to indicate if the given MatPoly is convex:
if set to false, it will be decomposed into tetrahedrons
@return True if the point is in the interior of MatPoly
*/
bool point_inside_matpoly(const MatPoly<3> mat_poly,
                          const Point<3>& pt, 
                          double dst_tol,
                          bool convex_poly) {
  std::vector< MatPoly<3> > convex_polys;
  if (convex_poly)
    convex_polys.push_back(mat_poly);
  else 
    mat_poly.facetize_decompose(convex_polys);

  for (int icp = 0; icp < convex_polys.size(); icp++) {
    const std::vector< Point<3> >& poly_pts = convex_polys[icp].points();

    bool pt_inside_cur_poly = true;
    for (int iface = 0; iface < convex_polys[icp].num_faces(); iface++) {
      const std::vector<int>& iface_vrts = convex_polys[icp].face_vertices(iface);
      Vector3 pt2vrt_vec;
      double pt2vrt_dst = 0.0;
      //Find a face vertex that is the farthest from the given point
      for (int ifvrt = 0; ifvrt < iface_vrts.size(); ifvrt++) {
        Vector3 cur_pt2vrt_vec = 
          convex_polys[icp].vertex_point(iface_vrts[ifvrt]) - pt;
        double cur_vec_norm = cur_pt2vrt_vec.norm();
        if (cur_vec_norm > pt2vrt_dst) {
          pt2vrt_vec = cur_pt2vrt_vec;
          pt2vrt_dst = cur_vec_norm;
        }
      }
      pt2vrt_vec /= pt2vrt_dst;

      Vector3 face_normal = polygon3d_normal(poly_pts, iface_vrts, dst_tol);
      if (face_normal.is_zero(dst_tol))
        continue;

      double dist2plane = dot(pt2vrt_vec, face_normal);
      //Check the signed distance to the plane
      if (dist2plane <= dst_tol) {
        pt_inside_cur_poly = false;
        break;
      }
    }  

    if (pt_inside_cur_poly)
      return true;
  }
  return false;
}

}  // namespace Tangram

#endif  // TANGRAM_MATPOLY_H_
