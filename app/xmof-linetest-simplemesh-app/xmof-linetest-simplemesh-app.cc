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

#include <fstream>
#include "tangram/support/tangram.h"
#include "tangram/wrappers/mesh/simple_mesh/simple_mesh_wrapper.h"
#include "tangram/simple_mesh/simple_mesh.h"
#include "tangram/driver/driver.h"
#include "tangram/reconstruct/xmof2D_wrapper.h"
#include "simple_vector.h"
#include "tangram/driver/write_to_gmv.h"

/* A simple test for a rectangular grid
   with a single linear material interface.
   Uses 2D SimpleMesh.
   Generates an nx x ny mesh, loads materials data from file,
   performs interface reconstruction and finds the max
   Hausdorff distance between recovered and reference
   material interface segments */

/* Refence material interface is given in the form
   y = mat_int_slope*x + mat_int_shift */
double mat_int_slope = 1.5;
double mat_int_shift = -0.25;
double deps = 1.0e-14;

void read_mat_data(const std::string& mesh_data_fname,
                   int nx, int ny,
                   std::vector<int>& cell_num_mats,
                   std::vector<int>& cell_mat_ids,
                   std::vector<double>& cell_mat_volfracs,
                   std::vector<Tangram::Point2>& cell_mat_centroids) {
  std::ifstream os(mesh_data_fname.c_str(), std::ifstream::binary);
  if (!os.good()) {
    std::ostringstream os;
    os << std::endl << "Cannot open " << mesh_data_fname <<
      " for binary input" << std::endl;
    throw XMOF2D::Exception(os.str());
  }

  int data_dim;
  os.read(reinterpret_cast<char *>(&data_dim), sizeof(int));
  assert(data_dim == 2);
  int ncells;
  os.read(reinterpret_cast<char *>(&ncells), sizeof(int));
  cell_num_mats.resize(ncells);
  
  if (nx*ny != ncells) {
    std::ostringstream os;
    os << std::endl << "Material data is provided for a mesh with " << ncells <<
      " instead of " << nx*ny << " cells!" << std::endl;
    throw XMOF2D::Exception(os.str());
  }
  
  std::vector<std::vector<int>> icell_mats(ncells);
  int nmatpoly = 0;
  for (int iy = 0; iy < ny; iy++)
    for (int ix = 0; ix < nx; ix++) {
      int icell = iy*nx + ix;
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
  for (int iy = 0; iy < ny; iy++)
    for (int ix = 0; ix < nx; ix++) {
      int icell = iy*nx + ix;
      if (cell_num_mats[icell] == 1) {
        cell_mat_volfracs[offset[icell]] = 1.0;
        continue;
      }
      for (int im = 0; im < cell_num_mats[icell]; im++)
        os.read(reinterpret_cast<char *>(&cell_mat_volfracs[offset[icell] + im]), sizeof(double));
    }
  cell_mat_centroids.resize(nmatpoly);
  for (int iy = 0; iy < ny; iy++)
    for (int ix = 0; ix < nx; ix++) {
      int icell = iy*nx + ix;
      if (cell_num_mats[icell] == 1) {
        cell_mat_centroids[offset[icell]] = Tangram::Point2(DBL_MAX, DBL_MAX);
        continue;
      }
      for (int im = 0; im < cell_num_mats[icell]; im++) {
        double cen_x, cen_y;
        os.read(reinterpret_cast<char *>(&cen_x), sizeof(double));
        os.read(reinterpret_cast<char *>(&cen_y), sizeof(double));
        cell_mat_centroids[offset[icell] + im] = Tangram::Point2(cen_x, cen_y);
      }
    }
  
  os.close();
}

int PosWRTLine(const XMOF2D::Point2D& p, double line_slope, double line_shift, double eps) {
  std::vector<double> lort_vec = {line_slope, -1.0};
  std::vector<double> p2l_vec = {p.x, p.y - line_shift};
  double prj = XMOF2D::ddot(lort_vec, p2l_vec);
  if (std::fabs(prj) < eps)
    return 0;
  return std::signbit(prj) ? -1 : 1;
}

int main(int argc, char** argv) {
#ifdef ENABLE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
#endif
  if (argc != 4) {
    std::ostringstream os;
    os << std::endl <<
    "Correct usage: xmof-linetest-simplemesh-app <mat_data_file> <nx> <ny>" << std::endl;
    throw XMOF2D::Exception(os.str());
  }

  std::string in_data_fname = argv[1];
  int nx = atoi(argv[2]);
  int ny = atoi(argv[3]);
  std::vector<double> xbnds = {0.0, 1.0};
  std::vector<double> ybnds = {0.0, 1.0};
 
  Tangram::Simple_Mesh mesh(xbnds[0], ybnds[0],
                            xbnds[1], ybnds[1],
                            nx, ny);
  Tangram::Simple_Mesh_Wrapper mesh_wrapper(mesh);

  int ncells = mesh_wrapper.num_owned_cells();

  std::vector<int> cell_num_mats;
  std::vector<int> cell_mat_ids;
  std::vector<double> cell_mat_volfracs;
  std::vector<Tangram::Point2> cell_mat_centroids;
  read_mat_data(in_data_fname, nx, ny, cell_num_mats,
                cell_mat_ids, cell_mat_volfracs, cell_mat_centroids);
  
  // Distance(angle) and volume fraction tolerances
  Tangram::IterativeMethodTolerances_t im_tols = {
    .max_num_iter = 1000, .arg_eps = 1.0e-14, .fun_eps = 1.0e-15};

  Tangram::Driver<Tangram::XMOF2D_Wrapper, 2,
    Tangram::Simple_Mesh_Wrapper> xmof_driver(mesh_wrapper, im_tols, true);
  
  xmof_driver.set_volume_fractions(cell_num_mats, cell_mat_ids, cell_mat_volfracs, cell_mat_centroids);
  xmof_driver.reconstruct();
  
  const std::vector<std::shared_ptr<Tangram::CellMatPoly<2>>>&
    cellmatpoly_list = xmof_driver.cell_matpoly_ptrs();
  
  Tangram::write_to_gmv(mesh_wrapper, 2, cell_num_mats, cell_mat_ids,
                        cellmatpoly_list, "out.gmv");
  
  std::vector<XMOF2D::Point2D> ref_line = {
    XMOF2D::Point2D(0.0, mat_int_shift),
    XMOF2D::Point2D(1.0, mat_int_slope + mat_int_shift)};
  
  double max_hausdorff = 0.0;
  int ncmp_int = 0;
  for (int icell = 0; icell < ncells; icell++) {
    std::vector<int> ifaces, fdirs;
    mesh_wrapper.cell_get_faces_and_dirs(icell, &ifaces, &fdirs);
    int nsides = (int) ifaces.size();
    bool should_be_mmc = false;
    std::vector<XMOF2D::Point2D> ref_int_pts;
    for (int iside = 0; iside < nsides; iside++) {
      std::vector<int> inodes;
      mesh_wrapper.face_get_nodes(ifaces[iside], &inodes);
      std::vector<XMOF2D::Point2D> side_vrts(2);
      for (int in = 0; in < 2; in++) {
        Tangram::Point2 node_coord;
        mesh_wrapper.node_get_coordinates(inodes[in], &node_coord);
        side_vrts[in] = XMOF2D::Point2D(node_coord[0], node_coord[1]);
      }
      std::vector<int> vrts_pos(2);
      for (int ivrt = 0; ivrt < 2; ivrt++)
        vrts_pos[ivrt] = PosWRTLine(side_vrts[ivrt], mat_int_slope, mat_int_shift, deps);
      
      if (!vrts_pos[0] || !vrts_pos[1]) {
        XMOF2D::Point2D ref_int = vrts_pos[0] ? side_vrts[1] : side_vrts[0];
        bool new_int_pt = true;
        for (int iip = 0; iip < ref_int_pts.size(); iip++)
          if (XMOF2D::distance(ref_int, ref_int_pts[iip]) < deps) {
            new_int_pt = false;
            break;
          }
        if (new_int_pt) {
          ref_int_pts.push_back(ref_int);
          continue;
        }
      }
      else if (vrts_pos[0] != vrts_pos[1]) {
        should_be_mmc = true;
        
        XMOF2D::Point2D ref_int = LineLineIntersect(ref_line, side_vrts);
        ref_int_pts.push_back(ref_int);
      }
    }
   
    const Tangram::CellMatPoly<2>* cell_mat_poly_ptr = cellmatpoly_list[icell].get(); 
    if (should_be_mmc) {
      if (ref_int_pts.size() != 2) {
        std::ostringstream os;
        os << std::endl <<
        "Interface line should have two intersections with a multi-material cell!" <<
        std::endl;
        throw XMOF2D::Exception(os.str());
      }
      
      int nfaces = cell_mat_poly_ptr->num_matfaces();
      std::vector<int> iintfaces;
      for (int iface = 0; iface < nfaces; iface++)
        if (cell_mat_poly_ptr->matface_is_interface(iface))
          iintfaces.push_back(iface);
      if (iintfaces.size() != 1)
        std::cout << "Cell #" << icell << " has " << iintfaces.size() <<
        " material interfaces instead of one!";
      else {
        ncmp_int++;
        std::vector<Tangram::Point2> int_face_vrts =
          cell_mat_poly_ptr->matface_points(iintfaces[0]);
      
        std::vector<XMOF2D::Segment> cmp_segs = {
          XMOF2D::Segment(ref_int_pts),
          XMOF2D::Segment(XMOF2D::Point2D(int_face_vrts[0][0], int_face_vrts[0][1]),
                          XMOF2D::Point2D(int_face_vrts[1][0], int_face_vrts[1][1])) };
        
        for (int iseg = 0; iseg < 2; iseg++)
          for (int ivrt = 0; ivrt < 2; ivrt++) {
            double cur_dist = cmp_segs[iseg].dist(cmp_segs[(iseg + 1)%2][ivrt]);
            if (cur_dist > max_hausdorff)
              max_hausdorff = cur_dist;
          }
      }
    }
    else {
      if (cell_mat_poly_ptr && (cell_mat_poly_ptr->num_matpolys() != 1))
        std::cout << "Cell #" << icell << " should be single-material, but contains " <<
        cell_mat_poly_ptr->num_matpolys() << " materials!" << std::endl;
    }
  }
  std::cout << "Checked Hausdorff distance for " << ncmp_int <<
    " in-cell material interfaces" << std::endl;
  std::cout << "Max Hausdorff distance between actual and reference " <<
    "material interface segments -> " << max_hausdorff << std::endl;
#ifdef ENABLE_MPI
  MPI_Finalize();
#endif
}
