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

#include "tangram/driver/CellMatPoly.h"

#include "gtest/gtest.h"
#include "mpi.h"

#include "Mesh.hh"
#include "MeshFactory.hh"

#include "tangram/support/tangram.h"

/// Test the CellMatPoly structure for one-dimensional cells

TEST(CellMatPoly, Mesh1D) {

  Jali::MeshFactory mf(MPI_COMM_WORLD);

  if (Jali::framework_available(Jali::Simple))
    mf.framework(Jali::Simple);
  mf.included_entities({Jali::Entity_kind::EDGE, Jali::Entity_kind::FACE});

  // Make a 1-cell one-dimensional mesh
  std::shared_ptr<Jali::Mesh> mesh = mf(0.0, 1.0, 2);

  // Make a 2-material CellMatPoly object for the cell

  Tangram::CellMatPoly<1> cellmatpoly;

  // Cell info
  int cellid = 0;
  cellmatpoly.set_cell(cellid);

  std::vector<int> cnodes(2);
  mesh->cell_get_nodes(cellid, &cnodes);

  std::vector<Tangram::Point<1>> cpoints(3);
  double x0, x1;
  mesh->node_get_coordinates(cnodes[0], &x0);
  cpoints[0] = Tangram::Point<1>(x0);
  mesh->node_get_coordinates(cnodes[1], &x1);
  cpoints[1] = Tangram::Point<1>(x1);
  cpoints[2] = cpoints[0] + 0.4*(cpoints[1]-cpoints[0]);

  // Add a two material polygons with material ID 0 and 1 for the cell

  cellmatpoly.add_matpoly(0, cpoints[0][0], cpoints[2][0],
                          Tangram::Entity_kind::NODE,
                          Tangram::Entity_kind::CELL,
                          cnodes[0], cellid);

  cellmatpoly.add_matpoly(1, cpoints[2][0], cpoints[1][0],
                          Tangram::Entity_kind::CELL,
                          Tangram::Entity_kind::NODE,
                          cellid, cnodes[1]);


  // Verify the matpoly info in cell 0

  ASSERT_EQ(cellid, cellmatpoly.cell());

  // There should be two material polygons, with material IDs 0 and 1

  ASSERT_EQ(2, cellmatpoly.num_matpolys());
  ASSERT_EQ(0, cellmatpoly.matpoly_matid(0));
  ASSERT_EQ(1, cellmatpoly.matpoly_matid(1));

  // Verify coordinates of vertices of matpoly 0 and 1
  std::vector<Tangram::Point<1>> mpoints;
  mpoints = cellmatpoly.matpoly_points(0);
  ASSERT_NEAR(cpoints[0][0], mpoints[0][0], 1e-16);
  ASSERT_NEAR(cpoints[2][0], mpoints[1][0], 1.0e-16);

  mpoints = cellmatpoly.matpoly_points(1);
  ASSERT_NEAR(cpoints[2][0], mpoints[0][0], 1.0e-16);
  ASSERT_NEAR(cpoints[1][0], mpoints[1][0], 1.0e-16);

  // Verify parent info for matpoly vertices

  ASSERT_EQ(Tangram::Entity_kind::NODE, cellmatpoly.matvertex_parent_kind(0));
  ASSERT_EQ(cnodes[0], cellmatpoly.matvertex_parent_id(0));
  ASSERT_EQ(Tangram::Entity_kind::CELL, cellmatpoly.matvertex_parent_kind(1));
  ASSERT_EQ(cellid, cellmatpoly.matvertex_parent_id(0));
  ASSERT_EQ(Tangram::Entity_kind::NODE, cellmatpoly.matvertex_parent_kind(2));
  ASSERT_EQ(cnodes[1], cellmatpoly.matvertex_parent_id(2));


  // Verify info for faces of matpoly 0 and 1

  std::vector<int> const& mfaces = cellmatpoly.matpoly_faces(0);
  ASSERT_EQ(2, mfaces.size());

  std::vector<Tangram::Point<1>> mfpoints = cellmatpoly.matface_points(0);
  ASSERT_EQ(1, mfpoints.size());
  ASSERT_NEAR(cpoints[0][0], mfpoints[0][0], 1.0e-16);

  std::vector<int> const& mfverts0 = cellmatpoly.matface_vertices(0);
  ASSERT_EQ(1, mfverts0.size());
  ASSERT_EQ(0, mfverts0[0]);

  mfpoints = cellmatpoly.matface_points(1);
  ASSERT_EQ(1, mfpoints.size());
  ASSERT_NEAR(cpoints[2][0], mfpoints[0][0], 1.0e-16);

  std::vector<int> const& mfverts1 = cellmatpoly.matface_vertices(1);
  ASSERT_EQ(1, mfverts1.size());
  ASSERT_EQ(1, mfverts1[0]);

  mfpoints = cellmatpoly.matface_points(2);
  ASSERT_EQ(1, mfpoints.size());
  ASSERT_NEAR(cpoints[1][0], mfpoints[0][0], 1.0e-16);

  std::vector<int> const& mfverts2 = cellmatpoly.matface_vertices(2);
  ASSERT_EQ(1, mfverts2.size());
  ASSERT_EQ(2, mfverts2[0]);

  // Are faces being correctly identified as being on the interface or not?

  ASSERT_TRUE(!cellmatpoly.matface_is_interface(0));
  ASSERT_TRUE(cellmatpoly.matface_is_interface(1));
  ASSERT_TRUE(!cellmatpoly.matface_is_interface(2));

  // Are faces connected to the correct matpolys?
  int matpoly0, matpoly1;
  cellmatpoly.matface_matpolys(0, &matpoly0, &matpoly1);
  ASSERT_EQ(-1, matpoly0);
  ASSERT_EQ(0, matpoly1);
  cellmatpoly.matface_matpolys(1, &matpoly0, &matpoly1);
  ASSERT_EQ(0, matpoly0);
  ASSERT_EQ(1, matpoly1);
  cellmatpoly.matface_matpolys(2, &matpoly0, &matpoly1);
  ASSERT_EQ(1, matpoly0);
  ASSERT_EQ(-1, matpoly1);

  // Parent info of faces
  ASSERT_EQ(Tangram::Entity_kind::NODE, cellmatpoly.matface_parent_kind(0));
  ASSERT_EQ(0, cellmatpoly.matface_parent_id(0));
  ASSERT_EQ(Tangram::Entity_kind::CELL, cellmatpoly.matface_parent_kind(1));
  ASSERT_EQ(0, cellmatpoly.matface_parent_id(1));
  ASSERT_EQ(Tangram::Entity_kind::NODE, cellmatpoly.matface_parent_kind(0));
  ASSERT_EQ(1, cellmatpoly.matface_parent_id(2));

  // Verify volumes

  ASSERT_NEAR(0.4*mesh->cell_volume(0), cellmatpoly.matpoly_volume(0),
              1.0e-08);
  ASSERT_NEAR(0.6*mesh->cell_volume(0), cellmatpoly.matpoly_volume(1),
              1.0e-08);

  // Verify centroids

  Tangram::Point<1> expcen0((x0 + (x0 + 0.4*(x1-x0)))/2.0);
  Tangram::Point<1> expcen1(((x0 + 0.4*(x1-x0)) + x1)/2.0);

  Tangram::Point<1> mcen0 = cellmatpoly.matpoly_centroid(0);
  ASSERT_NEAR(expcen0[0], mcen0[0], 1.0e-08);

  Tangram::Point<1> mcen1 = cellmatpoly.matpoly_centroid(1);
  ASSERT_NEAR(expcen1[0], mcen1[0], 1.0e-08);
}
