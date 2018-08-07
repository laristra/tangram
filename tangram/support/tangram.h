/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#ifndef TANGRAM_H_
#define TANGRAM_H_

#include <cfloat>
#include <vector>
#include <algorithm>
#include <memory>
#include <limits>

#ifdef THRUST

#include "thrust/device_vector.h"
#include "thrust/iterator/counting_iterator.h"
#include "thrust/transform.h"

#else  // no thrust

#include <boost/iterator/counting_iterator.hpp>

#endif

<<<<<<< HEAD
=======
#include <vector>
#include <algorithm>
#include <memory>
#include <limits>
>>>>>>> f67d12b5d67716e8a691b51fc28b0a9b367d7191
#include "tangram/support/Vector.h"

/*
  @file tangram.h
  @brief Several utility types and functions within the Tangram namespace.
 */

/*
  @namespace Tangram
  The Tangram namespace houses all of the code within Tangram.

  Cells (aka zones/elements) are the highest dimension entities in a mesh
  Nodes (aka vertices) are lowest dimension entities in a mesh
  Faces in a 3D mesh are 2D entities, in a 2D mesh are 1D entities
  BOUNDARY_FACE is a special type of entity that is need so that process
  kernels can define composite vectors (see src/data_structures) on
  exterior boundary faces of the mesh only

  Wedges are special subcell entities that are a simplicial
  decomposition of cell. In 3D, a wedge is tetrahedron formed by one
  point of the edge, the midpoint of the edge, the "center" of the
  face and the "center" of the cell volume. In 2D, a wedge is a
  triangle formed by an end-point of the edge, the mid-point of the
  edge and the center of the cell. In 1D, wedges are lines, that are
  formed by the endpoint of the cell and the midpoint of the
  cell. There are two wedges associated with an edge of cell face in
  3D. The two wedges associated with an edge and face pair form a larger
  simplex called a side

  Corners are also subcell entities that are associated uniquely with
  a node of a cell. Each corner is the union of all the wedges incident
  upon that node in the cell

  Facets are the boundary entity between two wedges in adjacent
  cells. In 3D, a facet is a triangular subface of the cell face
  shared by two wedges in adjacent cells. In 2D, a facet is half of
  an edge that is shared by two wedges in adjacent cells
 */

namespace Tangram {
// TODO:  Right now we're relying on the fact that this enum is
//        identical to Jali::Entity_kind.  Need to fix this.
/// The type of mesh entity.
enum Entity_kind {
  ALL_KIND = -3,     /*!< All possible types */
  ANY_KIND = -2,     /*!< Any of the possible types */
  UNKNOWN_KIND = -1, /*!< Usually indicates an error */
  NODE = 0,
  EDGE,
  FACE,
  CELL,
  SIDE,
  WEDGE,
  CORNER,
  FACET,
  BOUNDARY_FACE
};

const int NUM_ENTITY_KINDS = 8;

// Parallel status of entity
/// The parallel type of a given entity.
enum Entity_type {
  TYPE_UNKNOWN = -1,
  DELETED = 0,
  PARALLEL_OWNED = 1,   /*!< Owned by this processor */
  PARALLEL_GHOST = 2,   /*!< Owned by another processor */
  BOUNDARY_GHOST = 3,   /*!< Ghost/Virtual entity on boundary */
  ALL  = 4              /*!< PARALLEL_OWNED + PARALLEL_GHOST + BOUNDARY_GHOST */
};

/// Element (cell topology) type
enum Element_type {
  UNKNOWN_TOPOLOGY = 0,
  TRI,
  QUAD,
  POLYGON,
  TET,
  PRISM,
  PYRAMID,
  HEX,
  POLYHEDRON
};


#ifdef THRUST

template<typename T>
    using vector = thrust::device_vector<T>;

template<typename T>
    using pointer = thrust::device_ptr<T>;

typedef thrust::counting_iterator<int> counting_iterator;
inline counting_iterator make_counting_iterator(int const i) {
  return thrust::make_counting_iterator(i);
}

template<typename InputIterator, typename OutputIterator,
         typename UnaryFunction>
inline OutputIterator transform(InputIterator first, InputIterator last,
                                OutputIterator result, UnaryFunction op) {
  return thrust::transform(first, last, result, op);
}

template<typename InputIterator1, typename InputIterator2,
         typename OutputIterator, typename BinaryFunction>
inline OutputIterator transform(InputIterator1 first1, InputIterator1 last1,
                                InputIterator2 first2, OutputIterator result,
                                BinaryFunction op) {
  return thrust::transform(first1, last1, first2, result, op);
}

template<typename InputIterator, typename UnaryFunction>
inline void for_each(InputIterator first, InputIterator last,
                              UnaryFunction f) {
  thrust::for_each(first, last, f);
}

#else  // no thrust

template<typename T>
    using vector = std::vector<T>;

template<typename T>
    using pointer = std::shared_ptr<T>;

typedef boost::counting_iterator<int> counting_iterator;
inline counting_iterator make_counting_iterator(int const i) {
  return boost::make_counting_iterator<int>(i);
}

template<typename InputIterator, typename OutputIterator,
    typename UnaryFunction>
inline OutputIterator transform(InputIterator first, InputIterator last,
                                OutputIterator result, UnaryFunction op) {
  return std::transform(first, last, result, op);
}

template<typename InputIterator1, typename InputIterator2,
         typename OutputIterator, typename BinaryFunction>
inline OutputIterator transform(InputIterator1 first1, InputIterator1 last1,
                                InputIterator2 first2, OutputIterator result,
                                BinaryFunction op) {
  return std::transform(first1, last1, first2, result, op);
}

template<typename InputIterator, typename UnaryFunction>
inline void for_each(InputIterator first, InputIterator last,
                     UnaryFunction f) {
  std::for_each(first, last, f);
}


#endif

struct Weights_t {
  int entityID;
  std::vector<double> weights;
};

template <int D>
struct Plane_t {
  Vector<D> normal;
  double  dist2origin; // Distance from the plane (line) to the origin.
                       // If P is a point of a plane (line) and PO is the vector
                       // from P to the origin, then dist2origin = dot(PO, normal)
};

// Reverse plane normal
template<int D>
inline
const Plane_t<D> operator-(const Plane_t<D>& plane) {
  Plane_t<D> opp_dir_plane = {.normal = -plane.normal, 
                              .dist2origin = -plane.dist2origin};
  return opp_dir_plane;
}

template <int D>
class MatPoly;

template <int D>
struct MatPolySet_t {
  std::vector< MatPoly<D> > matpolys;  // MatPoly's in the set
  std::vector<double> moments;         // Aggregated moments of all MatPoly's
                                       // in the set.
};

template <int D>
struct HalfSpaceSets_t {
  MatPolySet_t<D> lower_halfspace_set;  // Set of MatPoly's below the line/plane
  MatPolySet_t<D> upper_halfspace_set;  // Set of MatPoly's above the line/plane
};

/* In the context of interface reconstruction methods based on the nested 
   dissections algorithm, the iterative methods are used to find the 
   position of the cutting plane. 
   VOF uses an iterative method to find the cutting distance for a given normal, 
   in which case the objective function is volume and fun_eps is volume tolerance,
   while arg_eps is unused.
   MOF, in addition to that, iteratively solves for the normal. In this case,
   the tolerances are used not for the objective function, but the argument,
   and arg_eps will be imposed on the change in polar angles between consecutive
   steps. Hence, MOF can use one instance of IterativeMethodTolerances_t with
   fun_eps used for finding the cutting distance, and arg_eps used for finding
   the direction of the normal.
*/   
struct IterativeMethodTolerances_t {
  int max_num_iter;   // Max number of iterations
  double arg_eps;     // Tolerance on the arguments of the function
  double fun_eps;     // Tolerance on the value of the function
};

template <int D>
struct BoundingBox_t {
  double min[D];
  double max[D];

  BoundingBox_t() {
    for (int idim = 0; idim < D; idim++) {
      min[idim] = DBL_MAX; max[idim] = -DBL_MAX;
    }
  }
};

template <int D>
inline bool overlapping_boxes(const BoundingBox_t<D>& bb1, 
                              const BoundingBox_t<D>& bb2,
                              double eps = std::numeric_limits<double>::epsilon()) {
  for (int idim = 0; idim < D; idim++)
    if ((bb1.min[idim] > bb2.max[idim] - eps) || 
        (bb1.max[idim] < bb2.min[idim] + eps))
      return false;

  return true;
}

// Check if two floating point values are equal up to the machine precision
inline bool is_equal(const double fp1, const double fp2) {
  return ( std::fabs(fp1 - fp2) < std::numeric_limits<double>::epsilon() );
}

}  // namespace Tangram

#endif  // TANGRAM_H_
