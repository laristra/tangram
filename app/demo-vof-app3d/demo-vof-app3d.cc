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
#include <sstream>

#include "tangram/support/tangram.h"
#include "tangram/simple_mesh/simple_mesh.h"
#include "tangram/wrappers/mesh/simple_mesh/simple_mesh_wrapper.h"
#include "tangram/intersect/split_r3d.h"
#include "tangram/driver/driver.h"
#include "tangram/reconstruct/VOF.h"
#include "tangram/driver/write_to_gmv.h"

using Tangram::Simple_Mesh;
using Tangram::Simple_Mesh_Wrapper;
using Tangram::Driver;
using Tangram::VOF;

#include <set>
std::vector<std::vector<double>> inputData(std::string fname) {
  // Output is a vector of vectors of volume fractions for each material
  std::vector<std::vector<double>> ret;
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
#ifdef ENABLE_MPI
  MPI_Init(&argc, &argv);
#endif
  // Read the input data
  // TODO - error checking on argv
  std::string fname = std::string("3d_diamond_6x6x6_vfracs.txt");
  std::string out_fname = std::string("3d_diamond_6x6x6.gmv");
  int nx(6), ny(6), nz(6);

  if (argc > 1) {
    if (argc < 6) {
      std::ostringstream os;
      os << std::endl <<
      "Correct usage: demo-vof-app <nx> <ny> <nz> " << 
      "<vol_fractions_filename> <out_gmv_filename>" << std::endl;
      throw std::runtime_error(os.str());
    }
    nx = atoi(argv[1]); ny = atoi(argv[2]); nz = atoi(argv[3]);
    
    fname = argv[4];
    if (argc > 5) out_fname = argv[5];
    else out_fname = fname + ".gmv";
  }

  auto vfracs = inputData(fname);

  auto numMats = vfracs.size();
  auto ncells = vfracs[0].size();

  assert(ncells == nx*ny*nz);

  // Simple_Mesh is only 3d, so we fake it
  auto mymesh = std::make_shared<Simple_Mesh>(0.0, 0.0, 0.0,
                                              1.0, 1.0, 1.0,
                                              nx, ny, nz);

  Simple_Mesh_Wrapper mymeshWrapper(*mymesh);

  // Volume fraction tolerance
  Tangram::IterativeMethodTolerances_t im_tols = {
    .max_num_iter = 1000, .arg_eps = 1.0e-13, .fun_eps = 1.0e-13};
  // Build the driver
  Driver<VOF, 3, Simple_Mesh_Wrapper, Tangram::SplitR3D, Tangram::ClipR3D> 
    d(mymeshWrapper, im_tols, true);

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

  std::vector<std::shared_ptr<Tangram::CellMatPoly<3>>> cellmatpoly_list = d.cell_matpoly_ptrs();

  // Confirm there are no degenerate faces
  for (int icell = 0; icell < ncells; icell++)
    if (cellmatpoly_list[icell] != nullptr) {
      const Tangram::CellMatPoly<3>& cellmatpoly = *cellmatpoly_list[icell];
      int npolys = cellmatpoly.num_matpolys();
      if (npolys == 1)
        std::cout << "Cell #" << icell << " has only one MatPoly!" << std::endl;
      for (int ipoly = 0; ipoly < npolys; ipoly++) {
        const std::vector<int>& mp_faces = cellmatpoly.matpoly_faces(ipoly);
        for (int iface = 0; iface < mp_faces.size(); iface++) {
          const std::vector<int>& face_vrts = cellmatpoly.matface_vertices(mp_faces[iface]);
          std::set<int> unique_face_vrts(face_vrts.begin(), face_vrts.end());
          if (unique_face_vrts.size() != face_vrts.size()) {
            std::cout << "Cell #" << icell << ", MatPoly #" << ipoly << ", side #" << iface <<
              ": repeated node indices!" << std::endl;
          }
        }
      }
    }

  // Filter out materials with volume fractions below tolerance
  std::vector<int> nzvf_cell_num_mats = cell_num_mats, nzvf_cell_mat_ids;
  int offset = 0, nzvf_offset = 0;
  for (int icell = 0; icell < ncells; icell++) {
    int ncmats = cell_num_mats[icell];
    for (int icmat = 0; icmat < ncmats; icmat++)
      if (cell_mat_volfracs[offset + icmat] > im_tols.fun_eps)
        nzvf_cell_mat_ids.push_back(cell_mat_ids[offset + icmat]);
      else
        nzvf_cell_num_mats[icell]--;

    //Create MatPoly's for single-material cells
    if (nzvf_cell_num_mats[icell] == 1) {
      assert(cellmatpoly_list[icell] == nullptr);
      std::shared_ptr< Tangram::CellMatPoly<3> > 
        cmp_ptr(new Tangram::CellMatPoly<3>(icell));
      Tangram::MatPoly<3> cell_matpoly;
      cell_get_matpoly(mymeshWrapper, icell, &cell_matpoly);
      cell_matpoly.set_mat_id(nzvf_cell_mat_ids[nzvf_offset]);
      cmp_ptr->add_matpoly(cell_matpoly);
      cellmatpoly_list[icell] = cmp_ptr;
    }
    
    offset += ncmats;
    nzvf_offset += nzvf_cell_num_mats[icell];
  }

  write_to_gmv(cellmatpoly_list, out_fname);

#ifdef ENABLE_MPI
  MPI_Finalize();
#endif
  return 0;
}
