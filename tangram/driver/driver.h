/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#ifndef TANGRAM_DRIVER_H_
#define TANGRAM_DRIVER_H_

#include <sys/time.h>

#include <algorithm>
#include <vector>
#include <iterator>
#include <string>
#include <utility>
#include <iostream>
#include <type_traits>

#ifdef WONTON_ENABLE_MPI
#include "mpi.h"
#endif

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

  @tparam  CellInterfaceReconstructor   A class that can reconstruct the pure material
  polygons in a cell given the volume fractions in the cell (and if needed
  its neighbors)
  @tparam  Mesh_Wrapper  A lightweight wrapper to a specific input mesh
  implementation that provides certain functionality.
*/

template <template <class, int, class, class> class CellInterfaceReconstructor,
    int Dim,
    class Mesh_Wrapper,
    class MatPoly_Splitter=void,
    class MatPoly_Clipper=void>
class Driver {
 public:
  /*!
    @brief Constructor for running the interface reconstruction driver.
    @param[in] Mesh @c Wrapper to the source mesh.
  */
  explicit Driver(Mesh_Wrapper const& Mesh,
                  const std::vector<IterativeMethodTolerances_t>& ims_tols,
                  const bool all_convex)
      : mesh_(Mesh), ims_tols_(ims_tols), all_convex_(all_convex) { }

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
    @brief Used iterative methods tolerances
    @return  Tolerances for iterative methods,
    here ims_tols_[0] correspond to methods for volumes
    and ims_tols_[1] correspond to methods for centroids.
    In particular, ims_tols_[0].arg_eps is a negligible
    change in cutting distance, ims_tols_[0].fun_eps is a
    negligible discrepancy in volume, ims_tols_[1].arg_eps
    is a negligible change in the cutting plane orientation,
    and ims_tols_[1].fun_eps is a negligible distance between
    actual and reference centroids. The change in cutting plane
    orientation is defined as the norm of change of the cutting
    plane's normal, which is expressed in polar coordinates (angles).
  */
  const std::vector<IterativeMethodTolerances_t>&
  iterative_methods_tolerances() const {
    return ims_tols_;
  }

  
  bool is_reconstruction_done() const { return reconstruction_done_; }

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
    int nc = mesh_.num_entities(Tangram::Entity_kind::CELL);
    assert(cell_num_mats.size() == unsigned(nc));
    
    cell_num_mats_ = cell_num_mats;
    cell_mat_ids_.clear();
    cell_mat_volfracs_.clear();
    cell_mat_centroids_.clear();

    // ims_tols_[0] are tolerances for 0-th moments, 
    // i.e. ims_tols_[0].fun_eps is the tolerance on absolute value of volume
    double volume_tol = ims_tols_[0].fun_eps;

    bool insufficient_vol_tol = false;
    // We only add materials with volumes above the tolerance
    int offset = 0;
    for (int icell = 0; icell < nc; icell++) {
      double cell_volume = mesh_.cell_volume(icell);  
      int ncmats = cell_num_mats_[icell];

      for (int icmat = 0; icmat < ncmats; icmat++) {
        double mat_volume = cell_volume*cell_mat_volfracs[offset + icmat];
        // quick fix for Marc Charest's segfault on very tiny material volumes,
        // do not update until a permanent solution is set (shift and scale).
        if (not(mat_volume < volume_tol)) {
          cell_mat_ids_.push_back(cell_mat_ids[offset + icmat]);
          cell_mat_volfracs_.push_back(cell_mat_volfracs[offset + icmat]);
          if (!cell_mat_centroids.empty())
            cell_mat_centroids_.push_back(cell_mat_centroids[offset + icmat]);

          // If specified volume tolerance is not too small, ensure that the volume
          // of the reconstructed material poly will not drop below the tolerance
          if (mat_volume <= 2*volume_tol)
            insufficient_vol_tol = true;
        }
        else
           cell_num_mats_[icell]--;
      }
      offset += ncmats;
    }

    if (insufficient_vol_tol) {
      auto& ims_tols_r = const_cast<std::vector<IterativeMethodTolerances_t>&>(ims_tols_);
      ims_tols_r[0].fun_eps /= 2;
    }
  }

  /*!
    @brief Perform the reconstruction
  */
  void reconstruct(Wonton::Executor_type const *executor = nullptr) {
#if !defined(NDEBUG) && defined(VERBOSE_OUTPUT)    
    int comm_rank = 0;
    int world_size = 1;

#ifdef WONTON_ENABLE_MPI

    MPI_Comm mycomm = MPI_COMM_NULL;
    auto mpiexecutor = dynamic_cast<Wonton::MPIExecutor_type const *>(executor);
    if (mpiexecutor && mpiexecutor->mpicomm != MPI_COMM_NULL) {
      mycomm = mpiexecutor->mpicomm;
      MPI_Comm_rank(mycomm, &comm_rank);
      MPI_Comm_size(mycomm, &world_size);
    }
#endif
#endif

    assert(!cell_num_mats_.empty());
    int ncells = mesh_.num_entities(Tangram::Entity_kind::CELL);
    
#if !defined(NDEBUG) && defined(VERBOSE_OUTPUT)
    if (comm_rank == 0) std::printf("in Driver::run()...\n");
#endif
    {
      // Instantiate the interface reconstructor class that will
      // compute the interfaces and compute the pure material submesh
      // in each cell
      CellInterfaceReconstructor<Mesh_Wrapper, Dim, MatPoly_Splitter, MatPoly_Clipper>
        reconstructor(mesh_, ims_tols_, all_convex_);

      // Tell the reconstructor what materials are in each cell and
      // what their volume fractions are
      reconstructor.set_volume_fractions(cell_num_mats_, cell_mat_ids_,
                                         cell_mat_volfracs_, cell_mat_centroids_);

#if !defined(NDEBUG) && defined(VERBOSE_OUTPUT)
      struct timeval begin_timeval, end_timeval, diff_timeval;
      float tot_seconds = 0.0;

      float xmat_cells_seconds = 0.0;
      int count = 2;

      struct timeval xmat_begin_timeval, xmat_end_timeval, xmat_diff_timeval;

      gettimeofday(&begin_timeval, 0);
#endif

      //Normally, we only need CellMatPoly's for multi-material cells,
      //so we first find their indices and group MMC's based on the number
      //of contained materials.
      //Because we only store indices of MMC's, iMMCs[0] vector corresponds
      //to two-material cells, and iMMCs[i] vector corresponds to MMC's with
      //(i+2) materials. If the partition contains MMC's with up to n_max
      //materials, the size of iMMCs vector is therefore (n_max-1).
      std::vector<std::vector<int>> iMMCs;
      for (int icell = 0; icell < ncells; icell++) {
        int nb_mats = cell_num_mats_[icell];
        int nb_mmcs = static_cast<int>(iMMCs.size());
        if (nb_mats < 2)
          continue;
        else if (nb_mats - 1 > nb_mmcs) {
          iMMCs.resize(nb_mats - 1);
        }
        iMMCs[nb_mats - 2].push_back(icell);
      }
      cellmatpolys_.resize(ncells);

      //Reconstructor is set to operate on multi-material cells only.
      //To improve load balancing, we operate on the cells with the same
      //number of materials at a time
      for (auto&& mm_cells : iMMCs) {
        int const nMMCs = mm_cells.size();
        if (nMMCs == 0)
          continue;

#if !defined(NDEBUG) && defined(VERBOSE_OUTPUT)        
        if (world_size == 1)
          gettimeofday(&xmat_begin_timeval, 0);
#endif

        reconstructor.set_cell_indices_to_operate_on(mm_cells);

        //If reconstruction is performed for a single cell, we do not use transform:
        //this allows to use transform inside the reconstructor (e.g. to run
        //multiple threads for material order permutations in MOF)
        if (nMMCs == 1)
          cellmatpolys_[mm_cells[0]] = reconstructor(0);
        else {
          Wonton::vector<std::shared_ptr<CellMatPoly<Dim>>> MMCs_cellmatpolys(nMMCs);

          Wonton::transform(Wonton::make_counting_iterator(0),
                            Wonton::make_counting_iterator(nMMCs),
                             MMCs_cellmatpolys.begin(), reconstructor);
          for (int immc = 0; immc < nMMCs; immc++)
            cellmatpolys_[mm_cells[immc]] = MMCs_cellmatpolys[immc];
        }

#if !defined(NDEBUG) && defined(VERBOSE_OUTPUT)
        if (world_size == 1) {
          gettimeofday(&xmat_end_timeval, 0);
          timersub(&xmat_end_timeval, &xmat_begin_timeval, &xmat_diff_timeval);
          xmat_cells_seconds = xmat_diff_timeval.tv_sec + 1.0E-6*xmat_diff_timeval.tv_usec;

          std::cout << "Transform Time for " << nMMCs << " "
                    << count++ << "-material cells (s): "
                    << xmat_cells_seconds << ", mean time per cell (s): "
                    << xmat_cells_seconds/nMMCs << std::endl;
        }
#endif
      }

#if !defined(NDEBUG) && defined(VERBOSE_OUTPUT)
      gettimeofday(&end_timeval, 0);
      timersub(&end_timeval, &begin_timeval, &diff_timeval);
      tot_seconds = diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;

      float max_transform_time = tot_seconds;
#ifdef WONTON_ENABLE_MPI
      if (world_size > 1) {
        MPI_Allreduce(&tot_seconds, &max_transform_time, 1, MPI_FLOAT, MPI_MAX,
          mycomm);
      }
#endif

      if (comm_rank == 0)
        std::cout << "Max Transform Time over " << world_size << " Ranks (s): " <<
          max_transform_time << std::endl;
#endif
    }
    
    reconstruction_done_ = true;
  }

  // /*!
  //   @brief Get a handle to an object that describes the pure
  //   material polygon topology in the cell

  //   @param  cellid  The cell for which we want the pure material polygon info

  //   @return An object that can be queried for details of the material
  //   polygons in the cell
  // */

  CellMatPoly<Dim> const& cell_matpoly_data(int const cellid) const {
    return *(cellmatpolys_[cellid].get());
  }
  
  const std::vector<std::shared_ptr<CellMatPoly<Dim>>>
    cell_matpoly_ptrs() const { return cellmatpolys_; }

 private:
  Mesh_Wrapper const& mesh_;
  const std::vector<IterativeMethodTolerances_t> ims_tols_; //Tolerances for iterative methods
                                                            //ims_tols_[0] for methods dealing with volumes
                                                            //ims_tols_[1] for methods dealing with centroids
  const bool all_convex_;
  std::vector<int> cell_num_mats_;
  std::vector<int> cell_mat_ids_;
  std::vector<double> cell_mat_volfracs_;
  std::vector<Point<Dim>> cell_mat_centroids_;

  std::vector<std::shared_ptr<CellMatPoly<Dim>>> cellmatpolys_;

  bool reconstruction_done_ = false;
};  // class Driver



}  // namespace Tangram

#endif  // TANGRAM_DRIVER_H_
