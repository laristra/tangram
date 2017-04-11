/*
Copyright (c) 2017, Los Alamos National Security, LLC
All rights reserved.

Copyright 2017. Los Alamos National Security, LLC. This software was produced
under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National
Laboratory (LANL), which is operated by Los Alamos National Security, LLC for
the U.S. Department of Energy. The U.S. Government has rights to use,
reproduce, and distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS
NATIONAL SECURITY, LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY
LIABILITY FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
derivative works, such modified software should be clearly marked, so as not to
confuse it with the version available from LANL.

Additionally, redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. Neither the name of Los Alamos National Security, LLC, Los Alamos
   National Laboratory, LANL, the U.S. Government, nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND
CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL
SECURITY, LLC OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
*/



#ifndef SRC_DRIVER_DRIVER_H_
#define SRC_DRIVER_DRIVER_H_

#include <sys/time.h>

#include <algorithm>
#include <vector>
#include <iterator>
#include <string>
#include <utility>
#include <iostream>
#include <type_traits>

#ifdef ENABLE_MPI
#include "mpi.h"
#endif

#include "tangram/support/Point.h"
#include "tangram/support/tangram.h"
#include "tangram/driver/CellMatPoly.h"

/*!
  @file driver.h
  @brief Example driver for computing interface reconstruction

  This should serve as a good example for how to write your own driver routine
  and data structures.
*/

namespace Tangram {


/*!
  @class Driver "driver.h"
  @brief Driver provides the API for reconstructing interfaces in
  multimaterial cells and answering queries about said interfaces

  @tparam  CellReconstruct   A class that can reconstruct the pure material
  polygons in a cell given the volume fractions in the cell (and if needed
  its neighbors)
  @tparam  Mesh_Wrapper  A lightweight wrapper to a specific input mesh
  implementation that provides certain functionality.
*/

template <template <class, int> class CellInterfaceReconstructor,
    int Dim,
    class Mesh_Wrapper>
class Driver {
 public:
  /*!
    @brief Constructor for running the interface reconstruction driver.
    @param[in] Mesh @c Wrapper to the source mesh.
  */
  explicit Driver(Mesh_Wrapper const& Mesh)
      : mesh_(Mesh)
  {}

  /// Copy constructor (disabled)
  Driver(const Driver &) = delete;

  /// Assignment operator (disabled)
  Driver & operator = (const Driver &) = delete;

  /// Destructor
  ~Driver() {}

  /*!
    @brief Get the dimensionality of the meshes.
    @return The dimensionality of the meshes.
  */
  unsigned int dim() const {
    return mesh_.space_dimension();
  }

  /*!
    @brief Send volume fraction (and optional centroid) data to the Driver
    @param cell_num_mats        Number of materials in each cell
    @param cell_mat_ids         Material IDs in each cell; length of
    sum(cell_num_mats)
    @param cell_mat_volfracs    Volume fractions of materials in each cell
    @param cell_mat_centroids   Centroids of materials in each cell
    @TODO actually use the centroid data if supplied
  */
  void set_volume_fractions(std::vector<int> const& cell_num_mats,
                            std::vector<int> const& cell_mat_ids,
                            std::vector<double> const& cell_mat_volfracs,
                            std::vector<Point<Dim>>
                            const& cell_mat_centroids = {}) {
    cell_num_mats_ = cell_num_mats;
    cell_mat_ids_ = cell_mat_ids;
    cell_mat_volfracs_ = cell_mat_volfracs;

    int nc = mesh_.num_entities(Tangram::Entity_kind::CELL);

    cell_mat_offsets_.resize(nc);
    cell_mat_offsets_[0] = 0;
    for (int c = 1; c < nc; c++)
      cell_mat_offsets_[c] = cell_mat_offsets_[c-1] + cell_num_mats_[c-1];
  }


  /*!
    @brief Perform the reconstruction
  */
  void reconstruct() {
    int comm_rank = 0;

#ifdef ENABLE_MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &comm_rank);
#endif

    if (comm_rank == 0) std::printf("in Driver::run()...\n");
    {
      float tot_seconds = 0.0, tot_seconds_srch = 0.0,
          tot_seconds_xsect = 0.0, tot_seconds_interp = 0.0;
      struct timeval begin_timeval, end_timeval, diff_timeval;

      gettimeofday(&begin_timeval, 0);

      // Instantiate the interface reconstructor class that will
      // compute the interfaces and compute the pure material submesh
      // in each cell
      CellInterfaceReconstructor<Mesh_Wrapper, Dim>
          reconstructor(mesh_);

      // Tell the reconstructor what materials are in each cell and
      // what their volume fractions are
      reconstructor.set_volume_fractions(cell_num_mats_, cell_mat_ids_,
                                         cell_mat_volfracs_);

      Tangram::transform(mesh_.begin(CELL, PARALLEL_OWNED),
                         mesh_.end(CELL, PARALLEL_OWNED),
                         cellmatpolys_.begin(), reconstructor);

      gettimeofday(&end_timeval, 0);
      timersub(&end_timeval, &begin_timeval, &diff_timeval);
      tot_seconds = diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;

      std::cout << "Transform Time Rank " << comm_rank << " (s): " <<
          tot_seconds << std::endl;
    }
  }

  // /*!
  //   @brief Get a handle to an object that describes the pure
  //   material polygon topology in the cell

  //   @param  cellid  The cell for which we want the pure material polygon info

  //   @return An object that can be queried for details of the material
  //   polygons in the cell
  // */

  // CellMatPoly const& cell_matpoly_data(int const cellid) {
  //   return cellmatpolys_[cellid];
  // }

  // /*!
  //   @brief Number of material polygons in cell
  //   @param cellid  ID of the cell we are querying
  //   @return Number of material polygons in cell which IN PRINCIPLE could
  //   be greater than the number of materials in the cell (for example,
  //   a layered reconstruction with Mat 1, Mat 2 and Mat 1)
  //  */

  // int num_matpolys(int const cellid) {
  // }

  // /*!
  //   @brief Which material does a material polygon in cell contain?
  //   @param cellid             ID of the cell we are querying
  //   @param matpoly_localid    Local polygon ID in the cell
  //   @return ID of material in the polygon
  // */

  // int matpoly_matid(int const cellid, int const matpoly_localid) const;

  // /*!
  //   @brief Volume of material polygon in cell
  //   @param cellid             ID of the cell we are querying
  //   @param matpoly_localid    Local polygon ID in the cell
  //   @return Volume of material polygon
  // */
  // double  matpoly_volume(int const cellid, int const matpoly_localid) const;

  // /*!
  //   @brief Centroid of material polygon in cell
  //   @param cellid             ID of the cell we are querying
  //   @param matpoly_localid    Local polygon ID in the cell
  //   @return Centroid of material polygon
  // */
  // Point<Dim> matpoly_centroid(int const cellid, int const matpoly_localid) const;

  // /*!
  //   @brief Points of the material polygon in cell
  //   @param cellid             ID of the cell we are querying
  //   @param matpoly_localid    Local polygon ID in the cell
  //   @param nv                 Number of points
  //   @param mpoints            Vector of points
  // */
  // void matpoly_points(int const cellid, int const matpoly_localid, int *nv,
  //                     std::vector<Point<Dim>> *mpoints) const;

  // /*!
  //   @brief Faces of the material polygon in cell
  //   @param cellid             ID of the cell we are querying
  //   @param matpoly_localid    Local polygon ID in the cell
  //   @param nf                 Number of faces
  //   @param mfaceids           Local IDs of material polygon faces (unique within
  //                             cell)
  // */
  // void matpoly_faces(int const cellid, int matpoly_localid,
  //                    int *nf, std::vector<int> *mfaceids) const;

  // /*!
  //   @brief Points of the material polygon in cell
  //   @param cellid             ID of the cell we are querying
  //   @param matface_localid    Local ID of material polygon face in the cell
  //   @param nv                 Number of points
  //   @param fpoints            Vector of points
  // */
  // void matface_points(int const cellid, int const matface_localid, int *nv,
  //                     std::vector<Point<D>> *fpoints) const;

  // /*!
  //   @brief "Vertices" of the material polygon in cell
  //   @param cellid             ID of the cell we are querying
  //   @param matface_localid    Local ID of material polygon face in the cell
  //   @param nv                 Number of vertices
  //   @param fvertices          Vector of local vertex IDs (shared across all
  //                             matpolys IN THE CELL)
  // */
  // void  matface_vertices(int const cellid, int const matface_localid,
  //                          int *nv, std::vector<int> *mvertex_localids) const;

  // /*!
  //   @brief Is material polygon face shared between two material polygons IN THE CELL
  //   @param cellid             ID of the cell we are querying
  //   @param matface_localid    Local ID of material polygon face in the cell
  //   @return True or False
  // */
  // bool  matface_is_interface(int const cellid, int const matface_localid) const;

  // /*!
  //   @brief Material polygons in cell on either side of material polygon face
  //   @param cellid             ID of the cell we are querying
  //   @param matface_localid    Local ID of material polygon face in the cell
  //   @param matpoly_localid0   Local ID of material polygon behind the face i.e. away from the normal (can be -1)
  //   @param matpoly_localid1   Local ID of material polygon in front of the face i.e. in the direction of the normal (can be -1)

  //   In 1D/2D, matpoly_localid0 is to the left of the matface and
  //   matpoly_localid1 is to the right

  //   CAVEAT: We cannot retrieve matpolys across cell boundaries
  // */
  // void  matface_matpolys(int const cellid, int const matface_localid,
  //                          int *matpoly_localid0, int *matpoly_localid1) const;

  // /*!
  //   @brief Kind of mesh entity that material face lies on
  //   @param cellid             ID of the cell we are querying
  //   @param matface_localid    Local ID of material polygon face in the cell
  //   @return Kind of mesh entity (can be FACE if its on the boundary of a cell or CELL if its in the interior of the cell)
  // */
  // Entity_kind matface_parent_kind(int const cellid, int const matface_localid) const;

  // /*!
  //   @brief ID of the mesh entity that material face lies on
  //   @param cellid             ID of the cell we are querying
  //   @param matface_localid    Local ID of material polygon face in the cell
  //   @return ID of mesh entity (global to the mesh? Local to the cell?)
  // */
  // int matface_parent_id(int const cellid, int const matface_localid) const;

  // /*!
  //   @brief Kind of mesh entity that material vertex lies on
  //   @param cellid             ID of the cell we are querying
  //   @param matvert_localid    Local ID of material polygon vertex in the cell
  //   @return Kind of mesh entity (can be FACE if its on the boundary of a cell or CELL if its in the interior of the cell)
  // */
  // Entity_kind  matvert_parent_kind(int const cellid, int const matvertid) const; // can be VERTEX/EDGE/FACE/CELL

  // /*!
  //   @brief ID of the mesh entity that material face lies on or in
  //   @param cellid             ID of the cell we are querying
  //   @param matvert_localid    Local ID of material polygon vertex in the cell
  //   @return ID of mesh entity (global to the mesh? Local to the cell?)
  // */
  // int  matvert_parent_id(int const cellid, int const matvertid) const;


 private:
  Mesh_Wrapper const& mesh_;
  std::vector<int> cell_num_mats_;
  std::vector<int> cell_mat_ids_;
  std::vector<double> cell_mat_volfracs_;
  std::vector<std::vector<double>> cell_mat_centroids_;
  std::vector<int> cell_mat_offsets_;

  // For now we are going to store a CellMatPoly data structure for all cells
  // even pure cells. We could store it only for material
  std::vector<CellMatPoly<Dim>> cellmatpolys_;

  bool reconstruction_done_ = false;
};  // class Driver



}  // namespace Tangram

#endif  // SRC_DRIVER_DRIVER_H_
