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

#include "tangram/driver/write_to_gmv.h"
#include "tangram/driver/CellMatPoly.h"
#include "tangram/wrappers/mesh/jali/jali_mesh_wrapper.h"
#include "tangram/support/tangram.h"

#include "gtest/gtest.h"
#include "mpi.h"

#include "Mesh.hh"             // Jali
#include "MeshFactory.hh"      // Jali

/// Test GMV Output of CellMatPoly structure for two-dimensional cells

TEST(WriteCellMatPoly, Mesh2D) {

  Jali::MeshFactory mf(MPI_COMM_WORLD);

  if (Jali::framework_available(Jali::MSTK))
    mf.framework(Jali::MSTK);
  mf.included_entities({Jali::Entity_kind::EDGE, Jali::Entity_kind::FACE});

  // Make a 4-cell two-dimensional mesh
  std::shared_ptr<Jali::Mesh> mesh = mf(0.0, 0.0, 2.0, 2.0, 2, 2);
  Tangram::Jali_Mesh_Wrapper mesh_wrapper(*mesh);

  // Make a 2-material CellMatPoly object for cell 0

  Tangram::CellMatPoly<2> cellmatpoly;

  // Cell info
  int cellid = 0;
  cellmatpoly.set_cell(cellid);

  std::vector<int> cnodes(4);
  mesh->cell_get_nodes(cellid, &cnodes);

  // Two material polygons 
  //
  //   o----------o
  //   |         /|
  //   |        / |
  //   |       /  |
  //   |      /   |
  //   o-----o----o


  // Specify the first polygon 

  std::vector<Tangram::Point<2>> points0(4);
  std::vector<Tangram::Entity_kind> vparentkind0(4);
  std::vector<int> vparentid0(4);

  points0[0] = Tangram::Point<2>(0.0, 0.0);
  vparentkind0[0] = Tangram::Entity_kind::NODE;
  vparentid0[0] = 0;                            // parent is node 0
  points0[1] = Tangram::Point<2>(0.5, 0.0);
  vparentkind0[1] = Tangram::Entity_kind::FACE;
  vparentid0[1] = 0;                            // parent is face 0
  points0[2] = Tangram::Point<2>(1.0, 1.0);
  vparentkind0[2] = Tangram::Entity_kind::NODE;
  vparentid0[2] = 4;                            // parent is node 4
  points0[3] = Tangram::Point<2>(0.0, 1.0);
  vparentkind0[3] = Tangram::Entity_kind::NODE;
  vparentid0[3] = 3;                            // parent is node 3

  cellmatpoly.add_matpoly(0, 4, &(points0[0]), &(vparentkind0[0]),
                          &(vparentid0[0]), nullptr, nullptr);


  // Specify the second polygon with parent info only for points

  std::vector<Tangram::Point<2>> points1(3);
  std::vector<Tangram::Entity_kind> vparentkind1(3);
  std::vector<int> vparentid1(3);

  points1[0] = Tangram::Point<2>(0.5, 0.0);
  vparentkind1[0] = Tangram::Entity_kind::FACE;
  vparentid1[0] = 0;                            // parent is face 0
  points1[1] = Tangram::Point<2>(1.0, 0.0);
  vparentkind1[1] = Tangram::Entity_kind::NODE;
  vparentid1[1] = 1;                            // parent is node 1
  points1[2] = Tangram::Point<2>(1.0, 1.0);
  vparentkind1[2] = Tangram::Entity_kind::NODE;
  vparentid1[2] = 4;                            // parent is node 2

  cellmatpoly.add_matpoly(1, 3, &(points1[0]), &(vparentkind1[0]),
                          &(vparentid1[0]), nullptr, nullptr);

  // Most cells have 1 material
  std::vector<int> cell_num_mats(4, 1);
  cell_num_mats[0] = 2;  // cell 0 has two materials

  std::vector<int> cell_mat_ids;
  cell_mat_ids.push_back(0);
  cell_mat_ids.push_back(1);
  cell_mat_ids.push_back(1);
  cell_mat_ids.push_back(0);
  cell_mat_ids.push_back(1);

  // Most cells don't have a CellMatPoly structure
  std::vector<std::shared_ptr<Tangram::CellMatPoly<2>>> cellmatpoly_list(4);
  cellmatpoly_list[0] = std::shared_ptr<CellMatPoly<2>>(&cellmatpoly);

  Tangram::write_to_gmv<Tangram::Jali_Mesh_Wrapper, 2>(mesh_wrapper, 2,
                                                       cell_num_mats,
                                                       cell_mat_ids,
                                                       cellmatpoly_list,
                                                       "interface_2D.gmv");
}
