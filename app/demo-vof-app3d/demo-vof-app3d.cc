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

#include "mpi.h"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "tangram/support/tangram.h"
#include "tangram/wrappers/mesh/jali/jali_mesh_wrapper.h"

#include "tangram/intersect/split_r3d.h"
#include "tangram/driver/driver.h"
#include "tangram/reconstruct/VOF.h"
#include "tangram/driver/write_to_gmv.h"
#include "app/include/read_material_data.h"

/* Demo app for an unstructured 3D mesh
   and a given material data.
   Uses Jali.
   Reads mesh and material data from file,
   performs interface reconstruction, and outputs material 
   polygons to a gmv file */

#include <set>

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  int world_size = 1;
  MPI_Comm_size(comm, &world_size);
  if (world_size > 1)
    throw std::runtime_error("This app is designed to run in serial!");

  if ((argc < 3) || (argc > 4)) {
      std::ostringstream os;
      os << std::endl <<
      "Correct usage: demo-vof-app <mat_data_filename> " << 
      "<base_mesh_file> <out_gmv_filename>" << std::endl;
      throw std::runtime_error(os.str());
  }

  std::string in_data_fname = argv[1];
  std::string out_gmv_fname;
  if (argc > 3) out_gmv_fname = argv[3];
  else out_gmv_fname = in_data_fname + ".gmv";

  Jali::MeshFactory mesh_factory(comm);
  mesh_factory.framework(Jali::MSTK);
  mesh_factory.included_entities({Jali::Entity_kind::EDGE, Jali::Entity_kind::FACE});
  std::shared_ptr<Jali::Mesh> mesh = mesh_factory(argv[2]);
  
  assert(mesh != nullptr);
  Tangram::Jali_Mesh_Wrapper mesh_wrapper(*mesh, true, false, false);

  int ncells = mesh_wrapper.num_owned_cells();

  std::vector<int> cell_num_mats;
  std::vector<int> cell_mat_ids;
  std::vector<double> cell_mat_volfracs;
  std::vector<Tangram::Point3> cell_mat_centroids;
  read_material_data<Tangram::Jali_Mesh_Wrapper, 3>(mesh_wrapper, in_data_fname, 
    cell_num_mats, cell_mat_ids, cell_mat_volfracs, cell_mat_centroids);

  // Volume fraction tolerance
  Tangram::IterativeMethodTolerances_t im_tols = {
    .max_num_iter = 1000, .arg_eps = 1.0e-13, .fun_eps = 1.0e-13};

  // Build the driver
  Tangram::Driver<Tangram::VOF, 3, Tangram::Jali_Mesh_Wrapper, 
                  Tangram::SplitR3D, Tangram::ClipR3D> 
    vof_driver(mesh_wrapper, im_tols, false);

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
