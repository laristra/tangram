/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/


#include <fstream>
#include <string>
#include <vector>
#include <sstream>

#include "mpi.h"

// Jali includes
#include "Mesh.hh"
#include "MeshFactory.hh"

// tangram includes
#include "tangram/support/tangram.h"
#include "tangram/driver/driver.h"
#include "tangram/driver/write_to_gmv.h"
#include "tangram/reconstruct/xmof2D_wrapper.h"

// wonton includes
#include "wonton/mesh/jali/jali_mesh_wrapper.h"

// XMOF includes
#include "simple_vector.h"

template <class Mesh_Wrapper>
void read_mat_data(const Mesh_Wrapper& Mesh,
                   const std::string& mesh_data_fname,
                   std::vector<int>& cell_num_mats,
                   std::vector<int>& cell_mat_ids,
                   std::vector<double>& cell_mat_volfracs,
                   std::vector<Tangram::Point2>& cell_mat_centroids) {
  cell_mat_ids.clear();
  cell_mat_volfracs.clear();
  cell_mat_centroids.clear();
  std::ifstream os(mesh_data_fname.c_str(), std::ifstream::binary);
  if (!os.good()) {
    std::ostringstream os;
    os << std::endl << "Cannot open " << mesh_data_fname <<
      " for binary input" << std::endl;
    throw XMOF2D::Exception(os.str());
  }
  int nrank_cells = Mesh.num_owned_cells() + Mesh.num_ghost_cells();
  cell_num_mats.resize(nrank_cells);
  std::vector<std::vector<int>> on_rank_mat_ids(nrank_cells);
  std::vector<std::vector<double>> on_rank_mat_volfracs(nrank_cells);
  std::vector<std::vector<Tangram::Point2>> on_rank_mat_centroids(nrank_cells);
  int data_dim;
  os.read(reinterpret_cast<char *>(&data_dim), sizeof(int));
  assert(data_dim == 2);
  int ncells;
  os.read(reinterpret_cast<char *>(&ncells), sizeof(int));
  std::vector<int> rank_cells_gid(nrank_cells);
  for (int irc = 0; irc < nrank_cells; irc++)
    rank_cells_gid[irc] = Mesh.get_global_id(irc, Tangram::Entity_kind::CELL);

  std::vector<int> on_rank_ids(ncells);
  std::vector<int> mesh_cell_num_mats(ncells);
  for (int icell = 0; icell < ncells; icell++) {
    on_rank_ids[icell] = std::distance(rank_cells_gid.begin(),
      std::find(rank_cells_gid.begin(), rank_cells_gid.end(), icell));

    int on_rank_id = on_rank_ids[icell];
    bool on_rank = (on_rank_id < nrank_cells);
    os.read(reinterpret_cast<char *>(&mesh_cell_num_mats[icell]), sizeof(int));
    if (on_rank) {
      cell_num_mats[on_rank_id] = mesh_cell_num_mats[icell];
      on_rank_mat_ids[on_rank_id].resize(cell_num_mats[on_rank_id]);
      on_rank_mat_volfracs[on_rank_id].resize(cell_num_mats[on_rank_id]);
      on_rank_mat_centroids.reserve(cell_num_mats[on_rank_id]);
      if (cell_num_mats[on_rank_id] == 1) {
        on_rank_mat_volfracs[on_rank_id][0] = 1.0;
        Tangram::Point2 cur_cell_cen;
        Mesh.cell_centroid(on_rank_id, &cur_cell_cen);
        on_rank_mat_centroids[on_rank_id].push_back(cur_cell_cen);
      }
    }
    for (int im = 0; im < mesh_cell_num_mats[icell]; im++) {
      int imat;
      os.read(reinterpret_cast<char *>(&imat), sizeof(int));
      if (on_rank)
        on_rank_mat_ids[on_rank_id][im] = imat;
    }
  }
  for (int icell = 0; icell < ncells; icell++) {
    int on_rank_id = on_rank_ids[icell];
    bool on_rank = (on_rank_id < nrank_cells);
    if (mesh_cell_num_mats[icell] > 1)
      for (int im = 0; im < mesh_cell_num_mats[icell]; im++) {
        double vfrac;
        os.read(reinterpret_cast<char *>(&vfrac), sizeof(double));
        if (on_rank)
          on_rank_mat_volfracs[on_rank_id][im] = vfrac;
      }
  }
  for (int icell = 0; icell < ncells; icell++) {
    int on_rank_id = on_rank_ids[icell];
    bool on_rank = (on_rank_id < nrank_cells);
    if (mesh_cell_num_mats[icell] > 1)
      for (int im = 0; im < mesh_cell_num_mats[icell]; im++) {
        double cen_x, cen_y;
        os.read(reinterpret_cast<char *>(&cen_x), sizeof(double));
        os.read(reinterpret_cast<char *>(&cen_y), sizeof(double));
        if (on_rank)
          on_rank_mat_centroids[on_rank_id].push_back(Tangram::Point2(cen_x, cen_y));
      }
  }
  os.close();
  for (int irc = 0; irc < nrank_cells; irc++) {
    cell_mat_ids.insert(cell_mat_ids.end(), on_rank_mat_ids[irc].begin(),
                        on_rank_mat_ids[irc].end());
    cell_mat_volfracs.insert(cell_mat_volfracs.end(), on_rank_mat_volfracs[irc].begin(),
                             on_rank_mat_volfracs[irc].end());
    cell_mat_centroids.insert(cell_mat_centroids.end(), on_rank_mat_centroids[irc].begin(),
                              on_rank_mat_centroids[irc].end());
  }
}

int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;
  if (argc != 3) {
    std::ostringstream os;
    os << std::endl <<
    "Correct usage: demo-xmof-app <mat_data_file> <base_mesh_file>" << std::endl;
    throw XMOF2D::Exception(os.str());
  }

  int comm_rank = 0;
  int world_size = 1;
  MPI_Comm_rank(comm, &comm_rank);
  MPI_Comm_size(comm, &world_size);

  if ((world_size > 1) && (comm_rank == 0)) {
    std::string err_msg = "Distributed runs are currently DISABLED ";
    err_msg += "as they require pre-partitioned material data files!";
    throw XMOF2D::Exception(err_msg);
  }

  std::string in_data_fname = argv[1];
  std::string out_gmv_fname = in_data_fname;
  out_gmv_fname.resize(out_gmv_fname.size() - 4);
  std::size_t path_end = out_gmv_fname.rfind('/');
  if (path_end != std::string::npos)
    out_gmv_fname = out_gmv_fname.substr(path_end + 1);
  out_gmv_fname += "_post_xmof_rank" + std::to_string(comm_rank) + ".gmv";
  Jali::MeshFactory mesh_factory(comm);
  mesh_factory.framework(Jali::MSTK);
  mesh_factory.included_entities({Jali::Entity_kind::EDGE, Jali::Entity_kind::FACE});
  std::shared_ptr<Jali::Mesh> mesh = mesh_factory(argv[2]);

  assert(mesh != nullptr);
  Wonton::Jali_Mesh_Wrapper mesh_wrapper(*mesh, true, false, false);
  int ncells = mesh_wrapper.num_entities(Tangram::Entity_kind::CELL);

  std::vector<int> cell_num_mats;
  std::vector<int> cell_mat_ids;
  std::vector<double> cell_mat_volfracs;
  std::vector<Tangram::Point2> cell_mat_centroids;
  read_mat_data(mesh_wrapper, in_data_fname, cell_num_mats, cell_mat_ids,
                cell_mat_volfracs, cell_mat_centroids);
  // Distance(angle) and volume tolerances
  std::vector<Tangram::IterativeMethodTolerances_t> ims_tols(2);
  ims_tols[0] = {.max_num_iter = 1000, .arg_eps = 1.0e-15, .fun_eps = 1.0e-15};
  ims_tols[1] = {.max_num_iter = 1000, .arg_eps = 1.0e-14, .fun_eps = 1.0e-14};

  Tangram::Driver<Tangram::XMOF2D_Wrapper, 2,
    Wonton::Jali_Mesh_Wrapper> xmof_driver(mesh_wrapper, ims_tols, true);

  xmof_driver.set_volume_fractions(cell_num_mats, cell_mat_ids, cell_mat_volfracs, cell_mat_centroids);
  xmof_driver.reconstruct();

  const std::vector<std::shared_ptr<Tangram::CellMatPoly<2>>>&
    cellmatpoly_list = xmof_driver.cell_matpoly_ptrs();

  std::vector<int> ncells_with_nmats(1, 0);
  std::vector<int> imaxcells;
  for (int icell = 0; icell < ncells; icell++) {
    const Tangram::CellMatPoly<2>* cell_mat_poly_ptr = cellmatpoly_list[icell].get();
    if (cell_mat_poly_ptr) {
      int cur_npolys = cell_mat_poly_ptr->num_matpolys();
      if (cur_npolys > ncells_with_nmats.size()) {
        ncells_with_nmats.resize(cur_npolys, 0);
        imaxcells.resize(1);
        imaxcells[0] = icell;
      }
      else if (cur_npolys == ncells_with_nmats.size())
        imaxcells.push_back(icell);

      ncells_with_nmats[cur_npolys - 1]++;
    }
  }

  int max_nmats = ncells_with_nmats.size();
  int top_max_nmats;
  MPI_Allreduce(&max_nmats, &top_max_nmats, 1, MPI_INT, MPI_MAX,
    MPI_COMM_WORLD);
  ncells_with_nmats.resize(top_max_nmats, 0);

  int* ncells_per_rank = NULL;
  int* nowned_cells_per_rank = NULL;
  if (comm_rank == 0) {
    ncells_per_rank = new int [world_size];
    nowned_cells_per_rank = new int [world_size];
  }
  int nowned_cells = mesh_wrapper.num_owned_cells();
  MPI_Gather(&ncells, 1, MPI_INT, ncells_per_rank, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Gather(&nowned_cells, 1, MPI_INT, nowned_cells_per_rank, 1, MPI_INT, 0, MPI_COMM_WORLD);
  std::vector<int*> cellmat_stats(top_max_nmats);
  for (int inm = 0; inm < top_max_nmats; inm++) {
    if (comm_rank == 0)
      cellmat_stats[inm] = new int [world_size];
    MPI_Gather(&ncells_with_nmats[inm], 1, MPI_INT, cellmat_stats[inm], 1, MPI_INT, 0,
      MPI_COMM_WORLD);
  }
  if (comm_rank == 0) {
    std::cout << "Number of ranks: " << world_size << std::endl;
    std::cout << "Total number of cells on each rank: ";
    for (int irank = 0; irank < world_size; irank++)
      std::cout << ncells_per_rank[irank] << " ";
    std::cout << std::endl;
    std::cout << "Number of owned cells on each rank: ";
    for (int irank = 0; irank < world_size; irank++)
      std::cout << nowned_cells_per_rank[irank] << " ";
    std::cout << std::endl;
    std::cout << "Number of single-material CellMatPoly's created on each rank: ";
    for (int irank = 0; irank < world_size; irank++)
      std::cout << cellmat_stats[0][irank] << " ";
    std::cout << std::endl;
    for (int inm = 0; inm < ncells_with_nmats.size() - 1; inm++) {
      std::cout << "Number of CellMatPoly's with " << inm + 2 << " materials created on each rank: ";
      for (int irank = 0; irank < world_size; irank++)
        std::cout << cellmat_stats[inm + 1][irank] << " ";
      std::cout << std::endl;
    }
    std::cout << std::endl;
    delete ncells_per_rank;
    delete nowned_cells_per_rank;
    for (int inm = 0; inm < top_max_nmats; inm++)
      delete cellmat_stats[inm];
  }
  if ( (world_size == 1) && (ncells_with_nmats.size() > 1) && (imaxcells.size() <= 25) ) {
    std::cout << "Cells with max number of material polygons are:" << std::endl;
    std::cout << "\t#CellID: [MaterialID, MaterialVolumeFraction] [..., ...]" << std::endl;
    for (int imc = 0; imc < imaxcells.size(); imc++) {
      double volume = mesh_wrapper.cell_volume(imaxcells[imc]);
      std::cout << "\t#" << imaxcells[imc] << ": ";
      for (int imp = 0; imp < ncells_with_nmats.size(); imp++) {
        const Tangram::CellMatPoly<2>* cell_mat_poly =
          cellmatpoly_list[imaxcells[imc]].get();
        std::cout << "[" << cell_mat_poly->matpoly_matid(imp) << ", " <<
          cell_mat_poly->matpoly_volume(imp)/volume << "]\t";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
  int nmats = *std::max_element(cell_mat_ids.begin(), cell_mat_ids.end()) + 1;
  Tangram::write_to_gmv(mesh_wrapper, nmats, cell_num_mats, cell_mat_ids,
                        cellmatpoly_list, out_gmv_fname);
  std::cout << "Resulting material polygons were written to " << out_gmv_fname << std::endl;
  MPI_Finalize();
}
