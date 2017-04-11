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

#ifndef SRC_RECONSTRUCT_SLIC_H_
#define SRC_RECONSTRUCT_SLIC_H_

#include <vector>

#include "tangram/support/Point.h"
#include "tangram/driver/CellMatPoly.h"

/*!
 @file SLIC.h
 @brief Simple implemenation of the crude Piecewise Linear Interface
 Reconstruction algorithm.

 Here, we use only vertical, "y"-aligned interfaces
*/

namespace Tangram {

/*!
  @class SLIC "SLIC.h"
  @brief Calculates the interface and constructs CellMatPoly for the SLIC
  algorithm.

  @tparam Mesh_Wrapper A lightweight wrapper to a specific input mesh
  implementation that provides certain functionality
  @tparam Dim The spatial dimension of the problem.
*/

template <class Mesh_Wrapper, int Dim>
class SLIC {
 public:
  /*!
    @brief Constructor performing a SLIC algorithm for interface reconstruction.
  */
  explicit SLIC(const Mesh_Wrapper & Mesh) : mesh_(Mesh) {
    // For now
    assert(Dim == 3);
  }

  /*!
    @brief Pass in the volume fraction data for use in the reconstruction.
    @param[in] cell_num_mats A vector of length (num_cells) specifying the
    number of materials in each cell.
    @param[in] cell_mat_ids A vector of length (sum(cell_num_mats)) specifying
    the ID of each material in each cell
    @param[in] cell_mat_volfracs A vector of length(sum(cell_num_mats))
    specifying the volume fraction of each material in each cell.
  */
  void set_volume_fractions(std::vector<int> const& cell_num_mats,
                            std::vector<int> const& cell_mat_ids,
                            std::vector<double> const& cell_mat_volfracs) {
    cell_num_mats_ = cell_num_mats;
    cell_mat_ids_ = cell_mat_ids_;
    cell_mat_volfracs_ = cell_mat_volfracs_;

    auto nc = mesh_.num_entities(Entity_kind::CELL);
    cell_mat_offsets_.resize(nc);
    cell_mat_offsets_[0] = 0;
    for (int c(1); c < nc; ++c)
      cell_mat_offsets_[c] = cell_mat_offsets_[c-1] + cell_num_mats_[c-1];
  }

  /*!
    @brief Given a cell index, calculate the CellMatPoly for this reconstruction
  */
  CellMatPoly<Dim> operator()(const int cellID) const {
    CellMatPoly<Dim> cellpoly(cellID);

    // We are assuming here - for now - that cells are rectangular prisms for
    // simplicity
    auto numMats = cell_num_mats_[cellID];
    auto iStart = cell_mat_offsets_[cellID];

    std::cout << "Working on cell " << cellID;
    std::cout << "  " << numMats << " " << iStart << std::endl;

    return cellpoly;
  }
private:
  const Mesh_Wrapper & mesh_;
  std::vector<int> cell_num_mats_;
  std::vector<int> cell_mat_ids_;
  std::vector<double> cell_mat_volfracs_;
  std::vector<int> cell_mat_offsets_;
};  // class SLIC
}  // namespace Tangram


#endif  // SRC_RECONSTRUCT_SLIC_H_
