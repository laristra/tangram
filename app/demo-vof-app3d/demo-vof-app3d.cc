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

/* Demo app for a rectangular grid
   and a given material data.
   Uses 3D SimpleMesh.
   Generates an nx x ny x nz mesh, loads material data from file,
   performs interface reconstruction, and outputs material 
   polygons to a gmv file */

#include <set>

void read_mat_data(const std::string& mesh_data_fname,
                   int nx, int ny, int nz,
                   std::vector<int>& cell_num_mats,
                   std::vector<int>& cell_mat_ids,
                   std::vector<double>& cell_mat_volfracs) {
  std::ifstream os(mesh_data_fname.c_str(), std::ifstream::binary);
  if (!os.good()) {
    std::ostringstream os;
    os << std::endl << "Cannot open " << mesh_data_fname <<
      " for binary input" << std::endl;
    throw std::runtime_error(os.str());
  }

  int data_dim;
  os.read(reinterpret_cast<char *>(&data_dim), sizeof(int));
  assert(data_dim == 3);
  int ncells;
  os.read(reinterpret_cast<char *>(&ncells), sizeof(int));
  cell_num_mats.resize(ncells);
  
  if (nx*ny*nz != ncells) {
    std::ostringstream os;
    os << std::endl << "Material data is provided for a mesh with " << ncells <<
      " instead of " << nx*ny*nz << " cells!" << std::endl;
    throw std::runtime_error(os.str());
  }
  
  std::vector<std::vector<int>> icell_mats(ncells);
  int nmatpoly = 0;
  for (int iz = 0; iz < nz; iz++)
    for (int iy = 0; iy < ny; iy++)
      for (int ix = 0; ix < nx; ix++) {
        int icell = iz*nx*ny + iy*nx + ix;
        os.read(reinterpret_cast<char *>(&cell_num_mats[icell]), sizeof(int));

        nmatpoly += cell_num_mats[icell];
        icell_mats[icell].resize(cell_num_mats[icell]);
        for (int im = 0; im < cell_num_mats[icell]; im++)
          os.read(reinterpret_cast<char *>(&icell_mats[icell][im]), sizeof(int));
      }

  std::vector<int> offset(ncells, 0);
  for (int icell = 0; icell < ncells - 1; icell++)
    offset[icell + 1] = offset[icell] + cell_num_mats[icell];
  cell_mat_ids.resize(nmatpoly);
  for (int icell = 0; icell < ncells; icell++)
    std::copy(icell_mats[icell].begin(), icell_mats[icell].end(),
              cell_mat_ids.begin() + offset[icell]);

  cell_mat_volfracs.resize(nmatpoly);
  for (int iz = 0; iz < nz; iz++)
    for (int iy = 0; iy < ny; iy++)
      for (int ix = 0; ix < nx; ix++) {
        int icell = iz*nx*ny + iy*nx + ix;
        if (cell_num_mats[icell] == 1) {
          cell_mat_volfracs[offset[icell]] = 1.0;
          continue;
        }
        for (int im = 0; im < cell_num_mats[icell]; im++)
          os.read(reinterpret_cast<char *>(&cell_mat_volfracs[offset[icell] + im]), sizeof(double));
      }
  
  os.close();
}

int main(int argc, char** argv) {
#ifdef ENABLE_MPI
  MPI_Init(&argc, &argv);
#endif

  if ((argc < 5) || (argc > 6)) {
      std::ostringstream os;
      os << std::endl <<
      "Correct usage: demo-vof-app <nx> <ny> <nz> " << 
      "<mat_data_filename> <out_gmv_filename>" << std::endl;
      throw std::runtime_error(os.str());
  }

  int nx = atoi(argv[1]);
  int ny = atoi(argv[2]);
  int nz = atoi(argv[3]);

  std::string in_data_fname = argv[4];
  std::string out_gmv_fname;
  if (argc > 5) out_gmv_fname = argv[5];
  else out_gmv_fname = in_data_fname + ".gmv";

  std::vector<double> xbnds = {0.0, 1.0};
  std::vector<double> ybnds = {0.0, 1.0};
  std::vector<double> zbnds = {0.0, 1.0};

  Tangram::Simple_Mesh mesh(xbnds[0], ybnds[0], zbnds[0],
                            xbnds[1], ybnds[1], zbnds[1],
                            nx, ny, nz);
  Tangram::Simple_Mesh_Wrapper mesh_wrapper(mesh);

  int ncells = mesh_wrapper.num_owned_cells();

  std::vector<int> cell_num_mats;
  std::vector<int> cell_mat_ids;
  std::vector<double> cell_mat_volfracs;
  read_mat_data(in_data_fname, nx, ny, nz, cell_num_mats, cell_mat_ids,
                cell_mat_volfracs);

  // Volume fraction tolerance
  Tangram::IterativeMethodTolerances_t im_tols = {
    .max_num_iter = 1000, .arg_eps = 1.0e-13, .fun_eps = 1.0e-13};

  // Build the driver
  Tangram::Driver<Tangram::VOF, 3, Tangram::Simple_Mesh_Wrapper, 
                  Tangram::SplitR3D, Tangram::ClipR3D> 
    vof_driver(mesh_wrapper, im_tols, true);

  vof_driver.set_volume_fractions(cell_num_mats, cell_mat_ids, cell_mat_volfracs);
  vof_driver.reconstruct();    

  std::vector<std::shared_ptr<Tangram::CellMatPoly<3>>> cellmatpoly_list = 
    vof_driver.cell_matpoly_ptrs();

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
      cell_get_matpoly(mesh_wrapper, icell, &cell_matpoly);
      cell_matpoly.set_mat_id(nzvf_cell_mat_ids[nzvf_offset]);
      cmp_ptr->add_matpoly(cell_matpoly);
      cellmatpoly_list[icell] = cmp_ptr;
    }
    
    offset += ncmats;
    nzvf_offset += nzvf_cell_num_mats[icell];
  }

  write_to_gmv(cellmatpoly_list, out_gmv_fname);

#ifdef ENABLE_MPI
  MPI_Finalize();
#endif
  return 0;
}
