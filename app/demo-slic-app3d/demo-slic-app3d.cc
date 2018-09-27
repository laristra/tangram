/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
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
#include "tangram/reconstruct/SLIC.h"
#include "tangram/driver/write_to_gmv.h"

using Tangram::Simple_Mesh;
using Tangram::Simple_Mesh_Wrapper;
using Tangram::Driver;
using Tangram::SLIC;


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
  if (argc > 1) {
    fname = argv[1];
    if (argc > 2) out_fname = argv[2];
    else out_fname = fname + ".gmv";
  }

  auto vfracs = inputData(fname);

  auto numMats = vfracs.size();
  auto ncells = vfracs[0].size();
  // Both demo problems have grids of this size.
  int nx(6), ny(6), nz(6);
  assert(ncells == nx*ny*nz);

  // Simple_Mesh is only 3d, so we fake it
  auto mymesh = std::make_shared<Simple_Mesh>(0.0, 0.0, 0.0,
                                              1.0, 1.0, 1.0,
                                              nx, ny, nz);

  Simple_Mesh_Wrapper mymeshWrapper(*mymesh);

  // Volume fraction tolerance
  std::vector<Tangram::IterativeMethodTolerances_t> ims_tols(1);
  ims_tols[0] = {.max_num_iter = 1000, .arg_eps = 1.0e-15, .fun_eps = 1.0e-15};

  // Build the driver
  Driver<SLIC, 3, Simple_Mesh_Wrapper, Tangram::SplitR3D> d(mymeshWrapper, ims_tols, true);

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
  write_to_gmv(cellmatpoly_list, out_fname);

#ifdef ENABLE_MPI
  MPI_Finalize();
#endif
  return 0;
}
