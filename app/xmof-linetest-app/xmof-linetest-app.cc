//
//  xmof-linetest-app.cc
//  Tangram + XMOF2D
//
//  Created by Evgeny Kikinzon.
//  Copyright © 2017 Los Alamos National Laboratory. All rights reserved.
//

#include <fstream>
#include "mpi.h"
#include "Mesh.hh"
#include "MeshFactory.hh"
#include "tangram/support/tangram.h"
#include "tangram/wrappers/mesh/jali/jali_mesh_wrapper.h"
#include "tangram/driver/driver.h"
#include "tangram/reconstruct/xmof2D_wrapper.h"
#include "simple_vector.h"

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
      int icell = ix*ny + iy;
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
      int icell = ix*ny + iy;
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
      int icell = ix*ny + iy;
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
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  
  if (argc != 4) {
    std::ostringstream os;
    os << std::endl <<
    "Correct usage: xmof-linetest-app <mat_data_file> <nx> <ny>" << std::endl;
    throw XMOF2D::Exception(os.str());
  }

  std::string in_data_fname = argv[1];
  int nx = atoi(argv[2]);
  int ny = atoi(argv[3]);
  std::vector<double> xbnds = {0.0, 1.0};
  std::vector<double> ybnds = {0.0, 1.0};
 
  Jali::MeshFactory mesh_factory(comm);
  mesh_factory.framework(Jali::MSTK);
  mesh_factory.included_entities({Jali::Entity_kind::EDGE, Jali::Entity_kind::FACE});
  std::shared_ptr<Jali::Mesh> mesh =
    mesh_factory(xbnds[0], ybnds[0], xbnds[1], ybnds[1], nx, ny);
  
  assert(mesh != nullptr);
  Tangram::Jali_Mesh_Wrapper mesh_wrapper(*mesh, true, false, false);
  int ncells = mesh_wrapper.num_owned_cells();

  std::vector<int> cell_num_mats;
  std::vector<int> cell_mat_ids;
  std::vector<double> cell_mat_volfracs;
  std::vector<Tangram::Point2> cell_mat_centroids;
  read_mat_data(in_data_fname, nx, ny, cell_num_mats,
                cell_mat_ids, cell_mat_volfracs, cell_mat_centroids);
  
  Tangram::Driver<Tangram::XMOF2D_Wrapper, 2,
    Tangram::Jali_Mesh_Wrapper> xmof_driver(mesh_wrapper);
  
  xmof_driver.set_volume_fractions(cell_num_mats, cell_mat_ids, cell_mat_volfracs, cell_mat_centroids);
  xmof_driver.reconstruct();
  
  const std::vector<std::shared_ptr<Tangram::CellMatPoly<2>>>&
    cellmatpoly_list = xmof_driver.cell_matpoly_ptrs();
  
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
        "Interfce line should have two intersections with a multi-material cell!" <<
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

  MPI_Finalize();
}
