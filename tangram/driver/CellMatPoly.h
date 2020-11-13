/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#ifndef TANGRAM_CELLMATPOLY_H_
#define TANGRAM_CELLMATPOLY_H_

#include <vector>
#include <array>
#include <algorithm>
#include <climits>
#include <sstream>

#include "tangram/support/tangram.h"
#include "tangram/support/MatPoly.h"

namespace Tangram {

// When this is in relatively good shape, we can write 1, 2 and 3
// dimensional specializations

/*
  Note that CellMatPoly expects the coordinates of coincident vertices of
  different material poly's to match exactly. If the coordinates are not 
  exactly equal, the respective vertices will be treated as separate nodes
  of CellMatPoly, which will also affect the topology: faces will not be 
  recognized as coincident and cell will not be recognized as adjacent.
*/

template <int D>
class CellMatPoly {
 public:
  /*!
    @brief Constructor with known cell ID
    @param cellid  ID of the cell
  */
  explicit CellMatPoly(int const cellid = -1) : cellid_(cellid) {}

  /*! Destructor */
  ~CellMatPoly() = default;

  /*!
    @brief Cell ID for which this object describes material polygons
    @return Cell ID
  */
  int cell() const {return cellid_;}

  /*!
    @brief Set the cell ID for which this object describes material polygons
    @param ID of the cell
  */
  void set_cell(int const cellid) {cellid_ = cellid;}

  /*!
    @brief Number of material polygons in cell
    @return Number of material polygons in cell which IN PRINCIPLE could
    be greater than the number of materials in the cell (for example,
    a layered reconstruction with Mat 1, Mat 2 and Mat 1)
   */
  int num_matpolys() const {return materialids_.size();}

  /*!
    @brief Number of distinct materials in cell
    @return Number of unique materials
  */
  int num_materials() const { return static_cast<int>(cell_materialids_.size()); }

  /*!
    @brief IDs of distinct materials in cell
    @return IDs of unique materials
  */
  const std::vector<int>& cell_matids() const { return cell_materialids_; }

  /*!
    @brief Check if cell contains a material with the given ID
    @param mat_id ID of the material
    @return True if such material is in cell
  */
  bool is_cell_material(const int mat_id) const { 
    return !(std::find(cell_materialids_.begin(), 
                       cell_materialids_.end(), mat_id) == 
             cell_materialids_.end()); 
  }

  /*!
    @brief Which material do material polygons in cell contain?
    @return Vector of IDs of material for all polygons in cell
  */
  const std::vector<int>& matpoly_matids() const { return materialids_; }
  
  /*!
    @brief Which material does a material polygon in cell contain?
    @param matpoly_id    Local polygon ID in the cell
    @return ID of material in the polygon
  */
  int matpoly_matid(int const matpoly_id) const {
    assert(num_matpolys_ > 0);
    assert(matpoly_id < num_matpolys_);
    return materialids_[matpoly_id];
  }

  /*!
    @brief Volume of material polygon in cell
    @param matpoly_id    Local polygon ID in the cell
    @return Volume of material polygon
  */
  double matpoly_volume(int const matpoly_id) const {
    assert(num_matpolys_ > 0);
    assert(matpoly_id < num_matpolys_);
    return matpoly_volumes_[matpoly_id];
  }

  /*!
    @brief Centroid of material polygon in cell
    @param matpoly_id    Local polygon ID in the cell
    @return Centroid of material polygon
  */
  Point<D> matpoly_centroid(int const matpoly_id) const {
    assert(num_matpolys_ > 0);
    assert(matpoly_id < num_matpolys_);
    return matpoly_centroids_[matpoly_id];
  }

  /*!
    @brief Vertices of the material polygon in cell
    @param matpoly_id    Local polygon ID in the cell
    @return              Vector of vertices of matpoly
  */
  std::vector<int> matpoly_vertices(int const matpoly_id) const {
    assert(num_matpolys_ > 0);
    assert(matpoly_id < num_matpolys_);
    std::vector<int> matverts;
    std::vector<int> const& matfaces = matpoly_faces_[matpoly_id];
    int nf = matfaces.size();
    for (int jf = 0; jf < nf; jf++) {
      int f = matfaces[jf];
      std::vector<int> const& fverts = matface_vertices_[f];
      if (D != 2) {  // order does not matter
        for (auto const& fv : fverts)
          if (std::find(matverts.begin(), matverts.end(), fv) == matverts.end())
            matverts.push_back(fv);
      } else {  // order is important
        if (matpoly_facedirs_[matpoly_id][jf] == 1) {  // natural order
          for (auto const& fv : fverts)
            if (std::find(matverts.begin(), matverts.end(), fv) == matverts.end())
              matverts.push_back(fv);
        } else {  // reverse the order
          for (auto it = fverts.rbegin(); it != fverts.rend(); it++) {
            auto const &fv = *it;
            if (std::find(matverts.begin(), matverts.end(), fv) == matverts.end())
              matverts.push_back(fv);
          }
        }
      }
    }
    return matverts;
  }


  /*!
    @brief Points of the material polygon in cell
    @param matpoly_id    Local polygon ID in the cell
    @return              Vector of points of matpoly
  */
  std::vector<Point<D>> matpoly_points(int const matpoly_id)  const {
    assert(num_matpolys_ > 0);
    assert(matpoly_id < num_matpolys_);
    std::vector<int> const& matpolyverts = matpoly_vertices(matpoly_id);
    int nv = matpolyverts.size();
    std::vector<Point<D>> matpolypnts(nv);
    for (int i = 0; i < nv; ++i)
      matpolypnts[i] = matvertex_points_[matpolyverts[i]];
    return matpolypnts;
  }

  /*!
    @brief Distance tolerance of CellMatPoly
    @return Distance tolerance for all contained MatPoly's
  */
  double dst_tol() {return dst_tol_;}  

  /*!
    @brief Faces of the material polygon in cell
    @param matpoly_id    Local polygon ID in the cell
    @return              Local IDs of matpoly faces (unique within cell)
  */
  std::vector<int> const& matpoly_faces(int matpoly_id) const {
    assert(num_matpolys_ > 0);    
    assert(matpoly_id < num_matpolys_);
    return matpoly_faces_[matpoly_id];
  }

  /*!
    @brief Number of unique faces among all material polygons
    @return Number of faces
  */
  int num_matfaces() const {
    return num_matfaces_;
  }
  
  /*!
    @brief Points of the material polygon in cell
    @param matface_id    Local ID of material polygon face in the cell
    @param nv                 Number of points
    @param fpoints            Vector of points
  */
  std::vector<Point<D>> matface_points(int const matface_id) const {
    assert(num_matpolys_ > 0);
    assert(matface_id < num_matfaces_);
    std::vector<int> const& matfaceverts = matface_vertices(matface_id);
    int nv = matfaceverts.size();
    std::vector<Point<D>> matfacepnts(nv);
    for (int i = 0; i < nv; ++i)
      matfacepnts[i] = matvertex_points_[matfaceverts[i]];
    return matfacepnts;
  }

  /*!
    @brief "Vertices" of the material polygon in cell
    @param matface_id    Local ID of material polygon face in the cell
    @param nv                 Number of vertices
    @param fvertices          Vector of local vertex IDs (shared across all
                              matpolys IN THE CELL)
  */
  std::vector<int> const& matface_vertices(int const matface_id) const {
    assert(num_matpolys_ > 0);
    assert(matface_id < num_matfaces_);
    return matface_vertices_[matface_id];
  }

  /*!
    @brief Is material polygon face shared between two material polygons IN THE CELL
    @param matface_id    Local ID of material polygon face in the cell
    @return True or False
  */
  bool  matface_is_interface(int const matface_id) const {
    assert(num_matpolys_ > 0);
    assert(matface_id < num_matfaces_);
    return (matface_matpolys_[matface_id][0] >= 0 &&
            matface_matpolys_[matface_id][1] >= 0);
  }

  /*!
    @brief Material polygons in cell on either side of material polygon face
    @param matface_id    Local ID of material polygon face in the cell
    @param matpoly_id0   Local ID of material polygon behind the face i.e. away from the normal (can be -1)
    @param matpoly_id1   Local ID of material polygon in front of the face i.e. in the direction of the normal (can be -1)

    In 1D/2D, matpoly_id0 is to the left of the matface and
    matpoly_id1 is to the right

    CAVEAT: We cannot retrieve matpolys across cell boundaries
  */
  void  matface_matpolys(int const matface_id,
                         int *matpoly_id0, int *matpoly_id1) const {
    assert(num_matpolys_ > 0);
    assert(matface_id < num_matfaces_);
    *matpoly_id0 = matface_matpolys_[matface_id][0];
    *matpoly_id1 = matface_matpolys_[matface_id][1];
  }

  /*!
    @brief Kind of mesh entity that material face lies on
    @param matface_id    Local ID of material polygon face in the cell
    @return Kind of mesh entity (can be FACE if its on the boundary of a cell or CELL if its in the interior of the cell)
  */
  Entity_kind matface_parent_kind(int const matface_id) const {
    assert(num_matpolys_ > 0);
    assert(matface_id < num_matfaces_);
    return matface_parentkind_[matface_id];
  }

  /*!
    @brief ID of the mesh entity that material face lies on
    @param matface_id    Local ID of material polygon face in the cell
    @return ID of mesh entity (global to the mesh? Local to the cell?)
  */
  int matface_parent_id(int const matface_id) const {
    assert(num_matpolys_ > 0);
    assert(matface_id < num_matfaces_);
    return matface_parentid_[matface_id];
  }

  /*!
    @brief Number of unique vertices among all material polygons
    @return Number of vertices
  */
  int num_matvertices() const {
    return num_matverts_;
  }
  
  Point<D> const& matvertex_point(int const matvert_id) const {
    assert(num_matpolys_ > 0);
    assert(matvert_id < num_matverts_);
    return matvertex_points_[matvert_id];
  }

  /*!
    @brief Kind of mesh entity that material vertex lies on
    @param matvert_localid    Local ID of material polygon vertex in the cell
    @return Kind of mesh entity (can be VERTEX/EDGE/FACE/CELL)
  */
  Entity_kind matvertex_parent_kind(int const matvert_id) const {
    assert(num_matpolys_ > 0);
    assert(matvert_id < num_matverts_);
    return matvertex_parentkind_[matvert_id];
  }

  /*!
    @brief ID of the mesh entity that material face lies on or in
    @param matvert_localid    Local ID of material polygon vertex in the cell
    @return ID of mesh entity (global to the mesh? Local to the cell?)
  */
  int matvertex_parent_id(int const matvert_id) const {
    assert(num_matpolys_ > 0);
    assert(matvert_id < num_matverts_);
    return matvertex_parentid_[matvert_id];
  }

  /*!
    @brief Add a material "polygon" to 1D cell
    @param matid              ID of material in polygon
    @param face_group_id      ID of the face group that material polygon is part of
    @param matpoly_point0     First point of material polygon
    @param matpoly_point1     Second point of material polygon
    @param point0_parentkind  Entity_kind of parent entity on which point 0 lies (can be UNKNOWN)
    @param point1_parentkind  Entity_kind of parent entity on which point 1 lies (can be UNKNOWN)
    @param point0_parentid    Local ID in cell of parent entity on which point 0 lies (can be -1)
    @param point1_parentid    Local ID in cell of parent entity on which point 1 lies (can be -1)
  */
  void add_matpoly(int matid,
                   int face_group_id,
                   double matpoly_point0, double matpoly_point1,
                   double dst_tol,
                   Entity_kind point0_parentkind,
                   Entity_kind point1_parentkind,
                   int point0_parentid, int point1_parentid);


  /*!
    @brief Add a material polygon to 2D cell
    @param matid              ID of material in polygon
    @param face_group_id      ID of the face group that material polygon is part of    
    @param numverts           Number of vertices of material polygon
    @param matpoly_points     List of points of material polygon listed ccw (numverts long)
    @param points_parentkind  Entity_kind of parent entity on which point lies (can be empty or numverts long)
    @param points_parentid    Local ID in cell of parent entity on which point lies (can be empty or numverts long)
  */
  void add_matpoly(int matid,
                   int face_group_id,
                   int numverts, Point<2> const * matpoly_points,
                   double dst_tol,
                   Entity_kind const * points_parentkind,
                   int const * points_parentid,
                   Entity_kind const * faces_parentkind,
                   int const * faces_parentid);

  /*!
    @brief Add a material polyhedron to 3D cell
    @param matid              ID of material in polygon
    @param face_group_id      ID of the face group that material polygon is part of    
    @param numvertices        Number of vertices of material polygon
    @param matpoly_points     List of points of material polygon listed ccw (numverts long)
    @param points_parentkind  Entity_kind of parent entity on which each point lies (can be empty or numverts long)
    @param points_parentid    Local ID in cell of parent entity on which each point lies (can be empty or numverts long)
    @param numfaces           Number of faces of material polygon
    @param numfaceverts       Number of vertices for each face of polyhedron (numfaces long)
    @param faces_point_ids    Index of face points in matpoly_points list (linear or flattened array)
    @param faces_parentkind   Entity_kind of parent entity on which each face lies (can be empty or numfaces long)
    @param faces_parentid     Local ID in cell of parent entity on which point lies (can be empty or numfaces long)
  */
  void add_matpoly(int matid,
                   int face_group_id,
                   int numverts, Point<3> const * matpoly_points,
                   double dst_tol,
                   Entity_kind const * points_parentkind,
                   int const * points_parentid,
                   int numfaces, int const * numfaceverts,
                   int const * faces_point_ids,
                   Entity_kind const * faces_parentkind,
                   int const * faces_parentid);
  
  /*!
   @brief Add a MatPoly to a cell
   @param mat_poly  MatPoly to add
  */
  void add_matpoly(const MatPoly<D>& mat_poly);

  /*!
   @brief Extracts material poly as a MatPoly object
   @param matpoly_id  ID of the material poly
   @return  Corresponding MatPoly object
  */
  MatPoly<D> get_ith_matpoly(int matpoly_id) const {
    std::stringstream ss;
    ss << "\nTangram does NOT support interface reconstruction for dimension " << D << std::endl;
    throw std::runtime_error(ss.str());
  }
  
  /*!
   @brief Extracts all polys containing a particular material
          as a vector of MatPoly objects
   @param mat_id  ID of the material for material polys to contain
   @return  Vector of MatPoly objects containing that material,
            vector is empty if no polys with that material are present
  */
  std::vector<MatPoly<D>> get_matpolys(int mat_id) const {
    std::vector<MatPoly<D>> mat_polys;
    for (int ipoly = 0; ipoly < num_matpolys_; ipoly++) {
      if (materialids_[ipoly] == mat_id)
        mat_polys.push_back(get_ith_matpoly(ipoly));
    }
    return mat_polys;
  }

  /*!
   @brief Extracts all polys belonging to a particular face group
          and containing a given material as a vector of MatPoly objects
   @param face_group_id  ID of the face group
   @param mat_id  ID of the material for material polys to contain
   @return  Vector of MatPoly objects belonging to a face group,
            vector is empty if no polys of that group are present
  */
  std::vector<MatPoly<D>> get_face_group_matpolys(int const face_group_id,
                                                  int const mat_id) const {
    assert(!face_group_ids_.empty());

    std::vector<MatPoly<D>> mat_polys;
    for (int ipoly = 0; ipoly < num_matpolys_; ipoly++) {
      if ((face_group_ids_[ipoly] == face_group_id) &&
          (materialids_[ipoly] == mat_id))
        mat_polys.push_back(get_ith_matpoly(ipoly));
    }
    return mat_polys;
  }

  /*!
   @brief Assigns externally computed aggregated moments for a distinct material
   @param mat_id  ID of the material
   @param moments Externally computed moments, moments[0] is the size, 
   moments[i+1]/moments[0] is i-th coordinate of the centroid
  */  
  void assign_material_moments(const int mat_id,
                               const std::vector<double>& moments) const {
    assert(moments.size() == D + 1);
    int cmp_mat_id = std::distance(cell_materialids_.begin(),
                                   std::find(cell_materialids_.begin(), 
                                             cell_materialids_.end(), mat_id));
    assert(cmp_mat_id < num_materials());

    if (static_cast<int>(material_moments_.size()) <= cmp_mat_id)
      material_moments_.resize(cmp_mat_id + 1);

    material_moments_[cmp_mat_id] = moments;
  }

  /*!
   @brief Aggregated moments of a distinct material, when the method is called 
   for the first time, moments will be computed and stored for future use
   @param mat_id  ID of the material
   @return  Vector of moments for the respective material
  */  
  const std::vector<double>& material_moments(const int mat_id) const {
    int cmp_mat_id = std::distance(cell_materialids_.begin(),
                                   std::find(cell_materialids_.begin(), 
                                             cell_materialids_.end(), mat_id));
    assert(cmp_mat_id < num_materials());

    if (static_cast<int>(material_moments_.size()) <= cmp_mat_id) {
      material_moments_.resize(cmp_mat_id + 1);
      compute_material_moments(cmp_mat_id);
    }
    else if (material_moments_[cmp_mat_id].empty())
      compute_material_moments(cmp_mat_id);

    return material_moments_[cmp_mat_id];
  }

 protected:
  /*!
   @brief Computes aggregated moments for a cell's material
   @param cmp_mat_id Local ID of a distinct cell's material
  */  
  void compute_material_moments(int cmp_mat_id) const {
    assert(unsigned(cmp_mat_id) < material_moments_.size());
    material_moments_[cmp_mat_id].assign(D + 1, 0.0);

    int mat_id = cell_materialids_[cmp_mat_id];

    for (int ipoly = 0; ipoly < num_matpolys_; ipoly++)
      if (materialids_[ipoly] == mat_id) {
        material_moments_[cmp_mat_id][0] += matpoly_volumes_[ipoly];
        for (int idim = 0; idim < D; idim++)
          material_moments_[cmp_mat_id][idim + 1] += 
            matpoly_volumes_[ipoly]*matpoly_centroids_[ipoly][idim];
      }
  }

 private:
  int cellid_ = -1;  // ID of the cell which we are describing

  // Topology
  int num_matpolys_ = 0;  // Number of material polygons
  //                      // Redundant but convenient!
  std::vector<int> materialids_;  // Material IDs of matpolys (can be repeated)
  std::vector<int> cell_materialids_;  // IDs of distinct materials in this CellMatPoly
  std::vector<int> face_group_ids_;  // IDs of faces associated with matpolys (can be repeated)
  mutable std::vector< std::vector<double> > material_moments_;  // Aggregated moments for distinct material
  std::vector<std::vector<int>> matpoly_faces_;  // IDs of faces of matpolys
  std::vector<std::vector<int>> matpoly_facedirs_;  // Dirs of faces of matpolys
  //                          // 1 means face normal points out of matpoly
  //                          // 0 means face normal points in to matpoly

  int num_matfaces_ = 0;  // Number of unique material polygon faces
  //                      // Redundant but convenient!
  std::vector<std::vector<int>> matface_vertices_;  // vertices of faces
  std::vector<std::array<int, 2>> matface_matpolys_;  // indices of matpolys that a face is connected to

  std::vector<Entity_kind> matface_parentkind_;  // Kind of parent mesh entity for polyhedron faces - FACE/CELL
  std::vector<int> matface_parentid_; // ID of parent mesh entity for polyhedral faces

  int num_matverts_ = 0;  // Number of unique material polygon vertices
  //                      // Redundant but convenient
  std::vector<Entity_kind> matvertex_parentkind_;  // Kind of parent mesh entity for polyhedron vertices - VERTEX/EDGE/FACE/CELL
  std::vector<int> matvertex_parentid_;  // ID of parent mesh entity for vertices

  // Geometry
  std::vector<Point<D>> matvertex_points_;  // coordinates of vertices
  std::vector<double> matpoly_volumes_;  // precomputed volumes of matpolys
  std::vector<Point<D>> matpoly_centroids_;  // precomputed centroids of matpolys
  double dst_tol_ = 0.0;
};  // class CellMatPoly


/*!
  @brief Add a material "polygon" to 1D cell
  @param matid              ID of material in polygon
  @param face_group_id      ID of the face group that material polygon is part of  
  @param matpoly_point0     First point of material polygon
  @param matpoly_point1     Second point of material polygon
  @param point0_parentkind  Entity_kind of parent entity on which point 0 lies (can be UNKNOWN)
  @param point1_parentkind  Entity_kind of parent entity on which point 1 lies (can be UNKNOWN)
  @param point0_parentid    Local ID in cell of parent entity on which point 0 lies (can be -1)
  @param point1_parentid    Local ID in cell of parent entity on which point 1 lies (can be -1)
*/
template<int D>
void CellMatPoly<D>::add_matpoly(int const matid,
                                 int const face_group_id,
                                 double const matpoly_point0,
                                 double const matpoly_point1,
                                 double dst_tol,
                                 Entity_kind const point0_parentkind,
                                 Entity_kind const point1_parentkind,
                                 int const point0_parentid,
                                 int const point1_parentid) {

  assert(D == 1);
  if (num_matpolys_ == 0)
    dst_tol_ = dst_tol;
  else {
    assert(dst_tol_ == dst_tol);
  }

  if (face_group_id == -1) {
    assert(face_group_ids_.empty());
  }
  else {
    assert(face_group_ids_.size() == unsigned(num_matpolys_));
    face_group_ids_.push_back(face_group_id);
  }

  int new_matpoly_id = num_matpolys_;
  if (std::find(cell_materialids_.begin(), cell_materialids_.end(), matid) == 
      cell_materialids_.end())
    cell_materialids_.push_back(matid);

  materialids_.push_back(matid);
  matpoly_faces_.resize(num_matpolys_+1);
  matpoly_facedirs_.resize(num_matpolys_+1);

  // Add matpoly points if they are not already there
  // In 1D, faces are the same as vertices/points

  bool found = false;
  int nv_ini = num_matverts_;
  for (int i = 0; i < nv_ini; i++) {
    Point<1>& p = matvertex_points_[i];
    if (std::fabs(p[0] - matpoly_point0) < dst_tol_) {
      found = true;  // This point is already in CellMatPoly
      matpoly_faces_[new_matpoly_id].push_back(i);
      matpoly_facedirs_[new_matpoly_id].push_back(0);

      assert(matface_matpolys_[i][1] == -1);  // verify slot is empty
      matface_matpolys_[i][1] = new_matpoly_id;  // new matpoly is to right of "face"
      break;
    }
  }
  if (!found) {  // Introduce this point into CellMatPoly
    matvertex_points_.push_back(Point<1>(matpoly_point0));
    matvertex_parentkind_.push_back(point0_parentkind);
    matvertex_parentid_.push_back(point0_parentid);

    int vert0 = num_matverts_, face0 = num_matverts_;
    matpoly_faces_[new_matpoly_id].push_back(face0);
    matpoly_facedirs_[new_matpoly_id].push_back(0);

    int new_matface_id = num_matfaces_;
    matface_matpolys_.emplace_back(std::array<int, 2>({-1, -1}));
    matface_matpolys_[new_matface_id][1] = new_matpoly_id;  // new matpoly is to right of "face"

    matface_vertices_.resize(num_matfaces_+1);
    matface_vertices_[new_matface_id].push_back(vert0);
    matface_parentkind_.push_back(point0_parentkind);
    matface_parentid_.push_back(point0_parentid);

    num_matverts_++;
    num_matfaces_++;
  }

  found = false;
  for (int i = 0; i < nv_ini; i++) {
    Point<1>& p = matvertex_points_[i];
    if (std::fabs(p[0] - matpoly_point1) < dst_tol_) {
      found = true;  // This point is already in CellMatPoly
      matpoly_faces_[new_matpoly_id].push_back(i);
      matpoly_facedirs_[new_matpoly_id].push_back(1);

      assert(matface_matpolys_[i][0] == -1);  // verify slot is empty
      matface_matpolys_[i][0] = new_matpoly_id;  // new matpoly is to left of "face"
      break;
    }
  }
  if (!found) {  // Introduce this point into CellMatPoly
    matvertex_points_.push_back(Point<1>(matpoly_point1));
    matvertex_parentkind_.push_back(point1_parentkind);
    matvertex_parentid_.push_back(point1_parentid);

    int vert1 = num_matverts_, face1 = num_matverts_;
    matpoly_faces_[new_matpoly_id].push_back(face1);
    matpoly_facedirs_[new_matpoly_id].push_back(1);

    int new_matface_id = num_matfaces_;

    matface_matpolys_.emplace_back(std::array<int, 2>({-1, -1}));
    matface_matpolys_[new_matface_id][0] = new_matpoly_id;  // new matpoly is to left of "face"

    matface_vertices_.resize(num_matfaces_+1);
    matface_vertices_[new_matface_id].push_back(vert1);
    matface_parentkind_.push_back(point1_parentkind);
    matface_parentid_.push_back(point1_parentid);

    num_matverts_++;
    num_matfaces_++;
  }

  Point<D> vec = matpoly_point1 - matpoly_point0;
  matpoly_volumes_.push_back(vec[0]);

  Point<D> cen = (matpoly_point0 + matpoly_point1)/2.0;
  matpoly_centroids_.push_back(cen);

  num_matpolys_++;
}  // add_matpoly for 1D



/*!
  @brief Add a material polygon to 2D cell
  @param matid              ID of material in polygon
  @param face_group_id      ID of the face group that material polygon is part of  
  @param numverts           Number of vertices
  @param matpoly_points     List of points of material polygon listed ccw (numverts long)
  @param points_parentkind  Entity_kind of parent entity on which point lies (can be empty or numverts long)
  @param points_parentid    Local ID in cell of parent entity on which point lies (can be empty or numverts long)
*/
template<int D>
void CellMatPoly<D>::add_matpoly(int const matid,
                                 int const face_group_id,
                                 int const numverts,
                                 Point<2> const * const matpoly_points,
                                 double dst_tol,
                                 Entity_kind const * const points_parentkind,
                                 int const * const points_parentid,
                                 Entity_kind const * const faces_parentkind,
                                 int const * const faces_parentid) {
  assert(D == 2);
  if (num_matpolys_ == 0) {
    dst_tol_ = dst_tol;
  }
  // avoid empty conditional branch in release mode
#ifndef NDEBUG
  else {
    assert(dst_tol_ == dst_tol);
  }  

  if (face_group_id == -1) {
    assert(face_group_ids_.empty());
  }
  else {
    assert(face_group_ids_.size() == unsigned(num_matpolys_));
#else
  if (face_group_id != -1) {
#endif
    face_group_ids_.push_back(face_group_id);
  }

  int new_matpoly_id = num_matpolys_;
  if (std::find(cell_materialids_.begin(), cell_materialids_.end(), matid) == 
      cell_materialids_.end())
    cell_materialids_.push_back(matid);

  materialids_.push_back(matid);

  std::vector<int> matverts(numverts);

  int nv_ini = num_matverts_;

  // Add matpoly points if they are not already there

  for (int j = 0; j < numverts; j++) {
    Point<D> const& pmat = matpoly_points[j];

    bool found = false;
    for (int i = 0; i < nv_ini; i++) {  // Note, nv remains number of original pnts
      Point<D>& p = matvertex_points_[i];
      // We use the provided distance tolerance to identify coincident nodes.
      if ((p - pmat).norm() < dst_tol_) {
        found = true;  // This point is already in CellMatPoly
        matverts[j] = i;
        break;
      }
    }

    if (!found) {  // Introduce this point into CellMatPoly
      matvertex_points_.push_back(pmat);
      Entity_kind parentkind = points_parentkind ? points_parentkind[j] :
          Entity_kind::UNKNOWN_KIND;
      matvertex_parentkind_.push_back(parentkind);
      int parentid = points_parentid ? points_parentid[j] : -1;
      matvertex_parentid_.push_back(parentid);
      matverts[j] = matvertex_points_.size()-1;
      num_matverts_++;
    }
  }

  // Add faces
  int numfaces = numverts;  // In 2D, number of "faces", which are really edges,
  //                        // is the same as number of vertices

  std::vector<int> matfaces(numfaces);
  std::vector<int> matfacedirs(numfaces);

  int nmf_ini = num_matfaces_;
  for (int j = 0; j < numfaces; j++) {
    std::vector<int> mfverts(2);
    mfverts[0] = matverts[j];
    mfverts[1] = matverts[(j+1)%numverts];

    bool found = false;
    for (int i = 0; i < nmf_ini; i++) {
      std::vector<int> const& mfverts2 = matface_vertices_[i];

      if (mfverts[0] == mfverts2[1] && mfverts[1] == mfverts2[0]) {
        found = true;
        matfaces[j] = i;
        matfacedirs[j] = 0;
        assert(matface_matpolys_[i][1] == -1);  // verify slot is empty
        matface_matpolys_[i][1] = new_matpoly_id;
        break;
      }
    }

    if (!found) {  // create a face
      matface_vertices_.push_back(mfverts);
      Entity_kind parentkind = faces_parentkind ? faces_parentkind[j] :
          Entity_kind::UNKNOWN_KIND;
      matface_parentkind_.push_back(parentkind);
      int parentid = faces_parentid ? faces_parentid[j] : -1;
      matface_parentid_.push_back(parentid);
      matfaces[j] = num_matfaces_;
      matfacedirs[j] = 1;

      matface_matpolys_.emplace_back(std::array<int, 2>({-1, -1}));
      matface_matpolys_[num_matfaces_][0] = new_matpoly_id;
      num_matfaces_++;
    }
  }

  matpoly_faces_.push_back(matfaces);
  matpoly_facedirs_.push_back(matfacedirs);


  // Geometric center of polygon; Assume that polygon is NOT so
  // non-convex that the center is outside the polygon

  Point<D> cen;
  for (int j = 0; j < numverts; j++)
    cen += matpoly_points[j];
  cen /= numverts;

  double volume = 0;
  Point<D> centroid;
  for (int j = 0; j < numverts; j++) {
    Vector<D> vec0 = matpoly_points[(j+1)%numverts] - matpoly_points[j];
    Vector<D> vec1 = cen - matpoly_points[j];
    double trivolume = 0.5*cross(vec0, vec1);
    volume += trivolume;
    Point<D> tricen = (matpoly_points[j] + matpoly_points[(j+1)%numverts] +
                       cen)/3.0;
    centroid += trivolume*tricen;
  }
  matpoly_volumes_.push_back(volume);

  centroid /= volume;
  matpoly_centroids_.push_back(centroid);

  num_matpolys_++;
}  // add_matpoly for 2D


/*!
  @brief Add a material polyhedron to 3D cell
  @param matid              ID of material in polygon
  @param face_group_id      ID of the face group that material polygon is part of  
  @param matpoly_points     List of points of material polygon listed ccw (numverts long)
  @param points_parentkind  Entity_kind of parent entity on which each point lies (can be empty or numverts long)
  @param points_parentid    Local ID in cell of parent entity on which each point lies (can be empty or numverts long)
  @param numfacepoints      Number of vertices for each face of polyhedron (numfaces long)
  @param faces_point_ids    Index of face points in matpoly_points list (linear or flattened array)
  @param faces_parentkind   Entity_kind of parent entity on which each face lies (can be empty or numfaces long)
  @param faces_parentid     Local ID in cell of parent entity on which point lies (can be empty or numfaces long)

  It is assumed that the vertices of the faces are specified such that
  the face normal points out of the material polyhedron being created
*/

template<int D>
void CellMatPoly<D>::add_matpoly(int const matid,
                                 int const face_group_id,
                                 int const numverts,
                                 Point<3> const * const matpoly_points,
                                 double dst_tol,
                                 Entity_kind const * const points_parentkind,
                                 int const * const points_parentid,
                                 int numfaces, int const * const numfaceverts,
                                 int const * const faces_point_ids,
                                 Entity_kind const * const faces_parentkind,
                                 int const * const faces_parentid) {

  assert(D == 3);
  if (num_matpolys_ == 0) {
    dst_tol_ = dst_tol;
  }
  // avoid empty conditional branch on release mode
#ifndef NDEBUG
  else {
    assert(dst_tol_ == dst_tol);
  }

  if (face_group_id == -1) {
    assert(face_group_ids_.empty());
  }
  else {
    assert(face_group_ids_.size() == unsigned(num_matpolys_));
#else
  if (face_group_id != -1) {
#endif
    face_group_ids_.push_back(face_group_id);
  }

  int new_matpoly_id = num_matpolys_;
  if (std::find(cell_materialids_.begin(), cell_materialids_.end(), matid) == 
      cell_materialids_.end())
    cell_materialids_.push_back(matid);


  materialids_.push_back(matid);

  std::vector<int> matverts(numverts);

  int nv_ini = matvertex_points_.size();

  // Add matpoly points if they are not already there

  for (int j = 0; j < numverts; j++) {
    Point<D> const & pmat = matpoly_points[j];

    bool found = false;
    
    for (int i = 0; i < nv_ini; i++) {
      Point<D>& p = matvertex_points_[i];
      // We use the provided distance tolerance to identify coincident nodes.
      if ((p - pmat).norm() < dst_tol_) {
        found = true;  // This point is already in CellMatPoly
        matverts[j] = i;
        break;
      }
    }

    if (!found) {  // Introduce this point into CellMatPoly
      matvertex_points_.push_back(pmat);
      Entity_kind parentkind = points_parentkind ? points_parentkind[j] :
          Entity_kind::UNKNOWN_KIND;
      matvertex_parentkind_.push_back(parentkind);
      int parentid = points_parentid ? points_parentid[j] : -1;
      matvertex_parentid_.push_back(parentid);
      matverts[j] = num_matverts_;
      num_matverts_++;
    }
  }

  // Add faces

  std::vector<int> matfaces(numfaces);
  std::vector<int> matfacedirs(numfaces);

  int nmf_ini = num_matfaces_;

  int offset = 0;
  for (int j = 0; j < numfaces; j++) {
    int nmfv = numfaceverts[j];

    std::vector<int> mfverts(nmfv);
    for (int i = 0; i < nmfv; i++)
      mfverts[i] = matverts[faces_point_ids[offset+i]];
    
    offset += nmfv;

    // First check if an existing face matches the points of this face.
    // Since we stipulated that the faces of the input material
    // polyhedron to this routine must be specified such that their
    // normal points out of the polyhedron, a face, if it exists, must
    // point out of a previously created polyhedron and into the one
    // that we are trying to create. Thus we check for a matching face
    // with vertices reversed

    bool found = false;
    for (int i = 0; i < nmf_ini; i++) {
      
      // Get vertices of existing face

      std::vector<int> const& mfverts2 = matface_vertices_[i];
      int nmfv2 = mfverts2.size();
      if (nmfv != nmfv2) continue;  // No match if number of verts don't match

      // Check if all vertices of the two faces match

      for (int k = 0; k < nmfv; k++) {
        if (mfverts[0] == mfverts2[k]) {

	  // Found one matching vertex
          // Check if the rest of the face vertices match (in reverse
          // since the two material polyhedra will use the face in
          // opposite directions)

	  bool allmatch = true;
          for (int m = 1; m < nmfv; m++) {
            if (mfverts[m] != mfverts2[(k-m+nmfv)%nmfv]) {
	      allmatch = false;
              break;
	    }
          }

	  // Found a matching face
          if (allmatch) {
	    found = true;
	    matfaces[j] = i;
	    matfacedirs[j] = 0;
	    assert(matface_matpolys_[i][1] == -1);
	    matface_matpolys_[i][1] = new_matpoly_id;
	    break;
	  }
        }
      }

      if (found)
	break;
    }  // for (int i = 0; i < ncmf_ini; i++)

    if (!found) {  // create a new matface
      matface_vertices_.push_back(mfverts);
      Entity_kind parentkind = faces_parentkind ? faces_parentkind[j] :
          Entity_kind::UNKNOWN_KIND;
      matface_parentkind_.push_back(parentkind);
      int parentid = faces_parentid ? faces_parentid[j] : -1;
      matface_parentid_.push_back(parentid);
      matfaces[j] = num_matfaces_;
      matfacedirs[j] = 1;

      matface_matpolys_.emplace_back(std::array<int, 2>({-1, -1}));
      matface_matpolys_[num_matfaces_][0] = new_matpoly_id;
      num_matfaces_++;
    }
  }

  matpoly_faces_.push_back(matfaces);
  matpoly_facedirs_.push_back(matfacedirs);


  // Compute geometric center of points. Assume that the cell is NOT so
  // distorted that this point is outside the cell

  Point<3> cen;
  for (int j = 0; j < numverts; j++)
    cen += matpoly_points[j];
  cen /= numverts;

  double volume = 0.0;
  Point<D> centroid;
  for (int j = 0; j < numfaces; j++) {
    Point<D> fcen;
    int f = matfaces[j];
    int dir = matfacedirs[j];
    int nmfv = matface_vertices_[f].size();
    for (int k = 0; k < nmfv; k++) {
      int k2 = matface_vertices_[f][k];
      fcen += matvertex_points_[k2];
    }
    fcen /= nmfv;

    for (int k = 0; k < nmfv; k++) {
      int v1 = matface_vertices_[f][k];
      int v2 = (dir == 1) ? matface_vertices_[f][(k+1)%nmfv] :
	matface_vertices_[f][(k-1+nmfv)%nmfv];

      Point<3> p1 = matvertex_points_[v1];
      Point<3> p2 = matvertex_points_[v2];
      Vector<3> vec0 =  p2 - p1;
      Vector<3> vec1 = fcen - p1;
      Vector<3> vec2 = cen - p1;
      Vector<3> inward_normal = -1.0*cross(vec0, vec1);
      double tetvolume = dot(inward_normal, vec2)/6.0;
      volume += tetvolume;
      Point<3> tetcentroid = (p1 + p2 + fcen + cen)/4.0;
      centroid += tetcentroid*tetvolume;
    }
  }

  matpoly_volumes_.push_back(volume);

  centroid /= volume;
  matpoly_centroids_.push_back(centroid);

  num_matpolys_++;
}  // add_matpoly for 3D

/*!
  @brief Add a MatPoly to a 2D cell
  @param mat_poly  2D MatPoly to add
*/
template <>
inline
void CellMatPoly<2>::add_matpoly(const MatPoly<2>& mat_poly) {
  add_matpoly(mat_poly.mat_id(), mat_poly.face_group_id(), mat_poly.num_vertices(), 
              &mat_poly.points()[0], mat_poly.dst_tol(), nullptr, nullptr,
              nullptr, nullptr);
}

/*!
  @brief Add a MatPoly to a 3D cell
  @param mat_poly  3D MatPoly to add
*/
template <>
inline
void CellMatPoly<3>::add_matpoly(const MatPoly<3>& mat_poly) {
  // Flatten the face vertices
  int nfaces = mat_poly.num_faces();
  std::vector<int> nface_vrts(nfaces);
  std::vector<int> faces_vrts;
  for (int iface = 0; iface < nfaces; iface++) {
    const std::vector<int>& face_ivrts = mat_poly.face_vertices(iface);
    nface_vrts[iface] = face_ivrts.size();
    faces_vrts.insert(faces_vrts.end(), face_ivrts.begin(), face_ivrts.end());
  }

  add_matpoly(mat_poly.mat_id(), mat_poly.face_group_id(), mat_poly.num_vertices(), 
              &mat_poly.points()[0], mat_poly.dst_tol(), nullptr, nullptr,
              nfaces, &nface_vrts[0], &faces_vrts[0], nullptr, nullptr);
}

/*!
  @brief Extracts 2D material polygon as a MatPoly object
  @param matpoly_id  ID of the material poly
  @return  Corresponding MatPoly object
*/
template<>
inline
MatPoly<2> CellMatPoly<2>::get_ith_matpoly(int matpoly_id) const {
#ifndef NDEBUG
  assert((matpoly_id >= 0) && (matpoly_id < num_matpolys_));
#endif
  std::vector<Point<2>> mp_pts = matpoly_points(matpoly_id);
  MatPoly<2> matpoly(matpoly_matid(matpoly_id));
  matpoly.initialize(mp_pts, dst_tol_);
  if (!face_group_ids_.empty())
    matpoly.set_face_group_id(face_group_ids_[matpoly_id]);
  
  return matpoly;
}

/*!
  @brief Extracts 3D material polyhedron as a MatPoly object
  @param matpoly_id  ID of the material poly
  @return  Corresponding MatPoly object
*/
template<>
inline
MatPoly<3> CellMatPoly<3>::get_ith_matpoly(int matpoly_id) const {
#ifndef NDEBUG
  assert((matpoly_id >= 0) && (matpoly_id < num_matpolys_));
#endif
  const std::vector<int>& mp_vrt_ids = matpoly_vertices(matpoly_id);
  int nvrts = static_cast<int>(mp_vrt_ids.size());
  std::vector<Point<3>> mp_pts;
  mp_pts.reserve(nvrts);
  for (int ivrt = 0; ivrt < nvrts; ivrt++)
    mp_pts.push_back(matvertex_points_[mp_vrt_ids[ivrt]]);
  
  const std::vector<int>& mp_faces = matpoly_faces(matpoly_id);
  int nfaces = static_cast<int>(mp_faces.size());
  std::vector<std::vector<int>> mf_vrts(nfaces);

  for (int iface = 0; iface < nfaces; iface++) {
    mf_vrts[iface] = matface_vertices(mp_faces[iface]);

    if (matpoly_facedirs_[matpoly_id][iface] == 0)
      std::reverse(mf_vrts[iface].begin(), mf_vrts[iface].end());

    for (auto&& vertex : mf_vrts[iface]) {
      int local_vrt_id =
        std::distance(mp_vrt_ids.begin(),
          std::find(mp_vrt_ids.begin(), mp_vrt_ids.end(), vertex));

      vertex = local_vrt_id;
    }
  }

  MatPoly<3> matpoly(matpoly_matid(matpoly_id));
  matpoly.initialize(mp_pts, mf_vrts, dst_tol_);
  if (!face_group_ids_.empty())
    matpoly.set_face_group_id(face_group_ids_[matpoly_id]);
  
  return matpoly;
}
  
}  // namespace Tangram

#endif  // TANGRAM_CELLMATPOLY_H_
