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

#include <stdlib.h>

#include <memory>
#include <vector>
#include <fstream>
#include <ostream>
#include <iterator>
#include <algorithm>
#include <string>

#include "tangram/support/tangram.h"
#include "tangram/simple_mesh/simple_mesh.h"
#include "tangram/wrappers/mesh/simple_mesh/simple_mesh_wrapper.h"
#include "tangram/driver/driver.h"
#include "tangram/reconstruct/SLIC.h"

using Tangram::Simple_Mesh;
using Tangram::Simple_Mesh_Wrapper;
using Tangram::Driver;
using Tangram::SLIC;


std::vector<std::vector<double>> inputData(const int probNum) {
  // Problem 0 - diamond
  // Problem 1 - eye
  // Output is a vector of vectors of volume fractions for each material
  std::vector<std::vector<double>> ret;

  auto fname = (probNum == 0) ?
    "diamond_20x20_vfracs.txt" : "eye_20x20_vfracs.txt";

  std::ifstream infile(fname);

  // Stupid parsing - there is certainly a better way
  if (infile.is_open()) {
    std::string line;
    getline(infile, line);
    int numMats = std::stoi(line);

    ret.resize(numMats);

    while (getline(infile, line)) {
      std::istringstream iss(line);

      std::vector<double> vfracsCell;
      std::copy(std::istream_iterator<double>(iss),
                std::istream_iterator<double>(),
                std::back_inserter(vfracsCell));

      assert(vfracsCell.size() == numMats);
      for (int iMat(0); iMat < numMats; ++iMat)
        ret[iMat].push_back(vfracsCell[iMat]);
    }
  } else {
    std::cout << "Failed to open file " << fname << std::endl;
    abort();
  }
  infile.close();

  return ret;
}

int main(int argc, char** argv) {
  // Read the input data
  // TODO - error checking on argv
  auto vfracs = inputData(atoi(argv[1]));

  auto numMats = vfracs.size();
  auto ncells = vfracs[0].size();
  // Both demo problems have grids of this size.
  int nx(20), ny(20);
  assert(ncells == nx*ny);

  // Simple_Mesh is only 3d, so we fake it
  auto mymesh = std::make_shared<Simple_Mesh>(0.0, 0.0, 0.0,
                                              1.0, 1.0, 1.0,
                                              nx, nx, 1);

  Simple_Mesh_Wrapper mymeshWrapper(*mymesh);

  // Build the driver
  Driver<SLIC, 3, Simple_Mesh_Wrapper> d(mymeshWrapper);

  // Load the volume fractions
  // I'm going to be dumb here - all cells will have all materials, even
  // if their volume fractions are zeros
  std::vector<int> cell_num_mats(ncells, numMats);
  std::vector<int> cell_mat_ids(ncells*numMats);
  std::vector<double> cell_mat_volfracs(ncells*numMats);
  for (int icell(0); icell < ncells; ++icell)
    for (int iMat(0); iMat < numMats; ++iMat) {
      cell_mat_ids[icell*numMats + iMat] = iMat;
      cell_mat_volfracs[icell*numMats + iMat] = vfracs[iMat][icell];
    }
  d.set_volume_fractions(cell_num_mats, cell_mat_ids, cell_mat_volfracs);
  d.reconstruct();

  // Do some dumb dumping of data for python plotting
  for (int c(0); c < ncells; ++c) {
    auto matpoly = d.cell_matpoly_data(c);
    auto numPolys = matpoly.num_matpolys();
    for (int ipoly(0); ipoly < numPolys; ++ipoly) {
      auto matID = matpoly.matpoly_matid(ipoly);
      auto nodeCoords = matpoly.matpoly_points(ipoly);
      std::cout << matID << " ";
      for (auto p : nodeCoords) std::cout << p << " ";
      std::cout << std::endl;
    }
  }

  return 0;
}
