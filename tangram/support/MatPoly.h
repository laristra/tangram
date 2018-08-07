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
#include "tangram/support/Point.h"

namespace Tangram {

double polygon3d_area(const std::vector<Point3>& points,
                      const std::vector<int>& poly_vertices) {
  double poly_area = 0.0;

  int nvrts = static_cast<int>(poly_vertices.size());
  if (nvrts == 3) {
    poly_area = 0.5*cross(points[poly_vertices[1]] - points[poly_vertices[0]], 
      points[poly_vertices[2]] - points[poly_vertices[0]]).norm();

    if (std::fabs(poly_area) <= std::numeric_limits<double>::epsilon())
      return 0.0;
  }

  Point3 gcenter;
  for (int ivrt = 0; ivrt < nvrts; ivrt++)
    gcenter += points[poly_vertices[ivrt]];
  gcenter /= nvrts;
  
  for (int ivrt = 0; ivrt < nvrts; ivrt++) {
    int ifv = poly_vertices[ivrt];
    int isv = poly_vertices[(ivrt + 1)%nvrts];
    double tri_size = 
      0.5*cross(points[isv] - points[ifv], gcenter - points[ifv]).norm();

    if (std::fabs(tri_size) > std::numeric_limits<double>::epsilon())
      poly_area += tri_size;
  }

  if (std::fabs(poly_area) <= std::numeric_limits<double>::epsilon())
    return 0.0;

  return poly_area;
}

void polygon3d_moments(const std::vector<Point3>& points,
                       const std::vector<int>& poly_vertices,
                       std::vector<double>& poly_moments) {
  poly_moments.assign(4, 0.0);

  int nvrts = static_cast<int>(poly_vertices.size());
  if (nvrts == 3) {
    poly_moments[0] = 0.5*cross(points[poly_vertices[1]] - points[poly_vertices[0]], 
      points[poly_vertices[2]] - points[poly_vertices[0]]).norm();

    if (std::fabs(poly_moments[0]) <= std::numeric_limits<double>::epsilon())
      poly_moments[0] = 0.0;
    else {
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
    double tri_size = 
      0.5*cross(points[isv] - points[ifv], gcenter - points[ifv]).norm();

    if (std::fabs(tri_size) > std::numeric_limits<double>::epsilon()) {
      poly_moments[0] += tri_size;
      for (int idim = 0; idim < 3; idim++)
        poly_moments[idim + 1] += 
          tri_size*(points[ifv][idim] + points[isv][idim] + gcenter[idim])/3.0;
    }
  }

  if (std::fabs(poly_moments[0]) <= std::numeric_limits<double>::epsilon())
    poly_moments.assign(4, 0.0);
}

Vector<3> polygon3d_normal(const std::vector<Point3>& points,
                           const std::vector<int>& poly_vertices) {
  std::vector<double> poly_moments;
  polygon3d_moments(points, poly_vertices, poly_moments);

  Vector3 poly_normal(0.0, 0.0, 0.0);
  if (poly_moments[0] <= std::numeric_limits<double>::epsilon())
    return poly_normal;

  Point3 poly_centroid;
  for (int idim = 0; idim < 3; idim++)
    poly_centroid[idim] = poly_moments[idim + 1]/poly_moments[0];

  int nvrts = static_cast<int>(poly_vertices.size());
  for (int ivrt = 0; ivrt < nvrts; ivrt++) {
    Vector3 tri_normal = cross(points[poly_vertices[ivrt]] - poly_centroid,
                               points[poly_vertices[(ivrt + 1)%nvrts]] - poly_centroid);
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
  */
  void initialize(const std::vector<Point2>& poly_points);
  /*!
   @brief Initialize a 3D polyhedron
   @param vertex_points  Vector of coordinates of material polyhedron's vertices
   @param face_vertices  Faces of the material polyhedron, every face is given by
   IDs of its vertices in counter-clockwise order
  */
  void initialize(const std::vector<Point3>& vertex_points,
                  const std::vector< std::vector<int> >& face_vertices);
  
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
    the MatPoly first. Faces with area below machine epsilon will be
    omitted. If the number of valid planes is less than four, empty
    vector will be returned.
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
void MatPoly<2>::initialize(const std::vector<Point2>& poly_points) {
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
                            const std::vector< std::vector<int> >& face_vertices) {
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
  polygon3d_moments(vertex_points_, face_vertices_[face_id], face_moments);

  if (std::fabs(face_moments[0]) > std::numeric_limits<double>::epsilon())
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
  
  faceted_poly->initialize(facetedpoly_vertices, facetedpoly_faces_);
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
    polygon3d_moments(vertex_points_, face_vertices_[iface], face_moments);

    if (std::fabs(face_moments[0]) <= std::numeric_limits<double>::epsilon())
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
      if (vcp.norm() <= std::numeric_limits<double>::epsilon())
        continue;

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
    sub_polys[offset + iface].initialize(subpoly_points);
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
    
    sub_polys[offset + iface].initialize(subpoly_vrts, subpoly_faces);
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
      sub_polys[icur_poly].initialize(tet_vrts, tet_faces);
    }
    else {
      tet_vrts[0] = face_centroid(iface);
      for (int ivrt = 0; ivrt < face_nvrts; ivrt++) {
        tet_vrts[1] = vertex_points_[face_vertices_[iface][ivrt]];
        tet_vrts[2] = vertex_points_[face_vertices_[iface][(ivrt + 1)%face_nvrts]];

        int icur_poly = static_cast<int>(sub_polys.size());
        sub_polys.push_back(MatPoly<3>(material_id_));
        sub_polys[icur_poly].initialize(tet_vrts, tet_faces);
      }
    }
  }
}

template <class Mesh_Wrapper>
void cell_get_matpoly(const Mesh_Wrapper& Mesh,
                      int const cellid,
                      MatPoly<2>* mat_poly) {
#ifdef DEBUG                        
  assert(Mesh.space_dimension() == 2);
#endif

  mat_poly->reset_mat_id();
  std::vector<Point2> poly_points;
  Mesh.cell_get_coordinates(cellid, &poly_points);

  mat_poly->initialize(poly_points);
}

template <class Mesh_Wrapper>
void cell_get_matpoly(const Mesh_Wrapper& Mesh,
                      int const cellid,
                      MatPoly<3>* mat_poly) {
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

  mat_poly->initialize(poly_points, poly_faces);
}

/*!
  @brief For every face, returns a line containing that face.
  Important: faces with length below machine epsilon will be
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
    if (normal_len > std::numeric_limits<double>::epsilon()) {
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
  the MatPoly first. Faces with area below machine epsilon will be
  omitted. If the number of valid planes is less than four, empty
  vector will be returned.
  @param fplanes  Vector of planes for all the faces of this MatPoly
*/
template <>
void MatPoly<3>::face_planes(std::vector< Plane_t<3> >& fplanes) const {
  fplanes.clear();

  for (int iface = 0; iface < nfaces_; iface++) {
    Plane_t<3> face_plane;
<<<<<<< HEAD
    face_plane.normal = polygon3d_normal(vertex_points_, 
                                         face_vertices_[iface]);
    if (face_plane.normal.is_zero())
=======
    int nfvrts = static_cast<int>(face_vertices_[iface].size());

    double normal_len = 0.0;
    std::vector<int> itri_pts(3);
    for (int ivrt = 0; ivrt < nfvrts; ivrt++) {
      itri_pts = { face_vertices_[iface][(ivrt + 1)%nfvrts],
        face_vertices_[iface][(ivrt + 2)%nfvrts], face_vertices_[iface][ivrt] };

      Vector3 cur_normal = cross(vertex_points_[itri_pts[1]] - vertex_points_[itri_pts[0]], 
                                 vertex_points_[itri_pts[2]] - vertex_points_[itri_pts[0]]);
      double cur_normal_len = cur_normal.norm();

      if (cur_normal_len > normal_len) {
        face_plane.normal = cur_normal;
        normal_len = cur_normal_len;
      }                      
    }

    // 0.5*normal_len is the area of the triangle formed 
    // by two vectors in the cross product
    if (normal_len <= 2*std::numeric_limits<double>::epsilon())
>>>>>>> f67d12b5d67716e8a691b51fc28b0a9b367d7191
      continue;
    
    face_plane.dist2origin = 
      -dot(vertex_points_[face_vertices_[iface][0]].asV(), face_plane.normal);

    fplanes.push_back(face_plane);  
  }

  if (fplanes.size() < 4)
    fplanes.clear();
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
                          const Point<3>& pt, bool convex_poly=false) {
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

      Vector3 face_normal = polygon3d_normal(poly_pts, iface_vrts);
      if (face_normal.is_zero())
        continue;

      double fnormal_prj = dot(pt2vrt_vec, face_normal);
      //Check the sign of the projection onto the normal
      if (fnormal_prj <= std::numeric_limits<double>::epsilon()) {
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
