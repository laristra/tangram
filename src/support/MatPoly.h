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

#include "tangram/support/tangram.h"
#include "tangram/support/Point.h"

namespace Tangram {

template <int D>
class MatPoly {
 public:
  /*!
    @brief Constructor, undefined material ID corresponds to -1 value
    @param material_id  ID of the material this poly contains
  */
  MatPoly(int const material_id = -1) : material_id_(material_id) { }

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
   @brief Moments of material poly
   @return  Vector of moments; moments[0] is the size, moments[i+1]/moments[0] is i-th
   coordinate of the centroid
  */  
  std::vector<double> moments();
 private:

  int material_id_;  // material ID of this matpoly
  std::vector< Point<D> > vertex_points_;  // coordinates of vertices
  std::vector< std::vector<int> > face_vertices_;  // vertices of faces
  
  int nvertices_;  // number of vertices
  int nfaces_;  //number of faces
};  // class MatPoly

/*!
 @brief Initialize a 2D polygon from its vertices
 @param poly_points  Vector of coordinates of material polygon's vertices
                     listed counter-clockwise
*/
template<>
void MatPoly<2>::initialize(const std::vector<Point2>& poly_points) {
  nvertices_ = (int) poly_points.size();
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
  nvertices_ = (int) vertex_points.size();
  nfaces_ = (int) face_vertices.size();
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
  Point3 centroid;

  int nvrts = (int) face_vertices_[face_id].size();
  if (nvrts == 3) {
    for (int ivrt = 0; ivrt < 3; ivrt++)
      centroid += vertex_points_[face_vertices_[face_id][ivrt]];
    centroid /= 3.0;
    
    return centroid;
  }

  Point3 gcenter;
  for (int ivrt = 0; ivrt < nvrts; ivrt++)
    gcenter += vertex_points_[face_vertices_[face_id][ivrt]];
  gcenter /= nvrts;
  
  double size = 0.0;
  for (int ivrt = 0; ivrt < nvrts; ivrt++) {
    int ifv = face_vertices_[face_id][ivrt];
    int isv = face_vertices_[face_id][(ivrt + 1)%nvrts];
    Vector3 vec0 = vertex_points_[isv] - vertex_points_[ifv];
    Vector3 vec1 = gcenter - vertex_points_[ifv];
    double tri_size = 0.5*cross(vec0, vec1).norm();
    size += tri_size;
    Point3 tri_centroid = (vertex_points_[ifv] + vertex_points_[isv] + gcenter)/3.0;
    centroid += tri_size*tri_centroid;
  }
  centroid /= size;

  return centroid;
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
  faceted_poly->set_mat_id(material_id_);
  
  std::vector<Point3> facetedpoly_vertices = vertex_points_;
  std::vector< std::vector<int> > facetedpoly_faces_;
  facetedpoly_faces_.reserve(nfaces_);
  for (int iface = 0; iface < nfaces_; iface++) {
    int nvrts = (int) face_vertices_[iface].size();
    if (nvrts == 3) {
      facetedpoly_faces_.push_back(face_vertices_[iface]);
      continue;
    }
    int icenvrt = (int) facetedpoly_vertices.size();
    facetedpoly_vertices.push_back(face_centroid(iface));
    for (int ivrt = 0; ivrt < nvrts; ivrt++)
      facetedpoly_faces_.push_back({ icenvrt, face_vertices_[iface][ivrt],
                                     face_vertices_[iface][(ivrt + 1)%nvrts] });
  }
  
  faceted_poly->initialize(facetedpoly_vertices, facetedpoly_faces_);
}

/*!
  @brief Moments of material polygon
  @return  Vector of moments; moments[0] is area, moments[i+1]/moments[0] is i-th
  coordinate of the centroid, i=1,2
*/   
template<>
std::vector<double> MatPoly<2>::moments() {
  std::vector<double> mpoly_moments(3, 0.0);
  for (int ivrt = 0; ivrt < nvertices_; ivrt++) {
    double cur_term = vertex_points_[ivrt][0]*vertex_points_[(ivrt + 1)%nvertices_][1] - 
                      vertex_points_[ivrt][1]*vertex_points_[(ivrt + 1)%nvertices_][0];
    mpoly_moments[0] += cur_term;
    for (int idim = 0; idim < 2; idim++)
      mpoly_moments[idim + 1] += cur_term*(
        vertex_points_[ivrt][idim] + vertex_points_[(ivrt + 1)%nvertices_][idim]);
  }
  mpoly_moments[0] /= 2.0;  
  for (int idim = 0; idim < 2; idim++)
    mpoly_moments[idim + 1] /= 6.0;

  return mpoly_moments;
}

/*!
  @brief Moments of material polyhedron
  @return  Vector of moments; moments[0] is volume, moments[i+1]/moments[0] is i-th
  coordinate of the centroid, i=1,2,3
*/   
template<>
std::vector<double> MatPoly<3>::moments() {
  std::vector<double> mpoly_moments(4, 0.0);
  MatPoly<3> faceted_poly;
  faceted_matpoly(&faceted_poly);
  int nfaces = faceted_poly.num_faces();
  for (int iface = 0; iface < nfaces; iface++) {
    std::vector<Point3> tri_pts;
    tri_pts.reserve(3);
    for (int ivrt = 0; ivrt < 3; ivrt++)
      tri_pts.push_back(faceted_poly.points()[faceted_poly.face_vertices(iface)[ivrt]]);
    Vector3 vcp = cross(tri_pts[1] - tri_pts[0], tri_pts[2] - tri_pts[0]);
    mpoly_moments[0] += dot(vcp, tri_pts[0].asV());
    for (int idim = 0; idim < 3; idim++)
      for (int ivrt = 0; ivrt < 3; ivrt++)
        mpoly_moments[idim + 1] += vcp[idim]*pow(tri_pts[ivrt][idim] + tri_pts[(ivrt + 1)%3][idim], 2);
  }
  mpoly_moments[0] /= 6.0;
  for (int idim = 0; idim < 3; idim++)
    mpoly_moments[idim + 1] /= 48.0;

  return mpoly_moments;  
}
}  // namespace Tangram

#endif  // TANGRAM_MATPOLY_H_
