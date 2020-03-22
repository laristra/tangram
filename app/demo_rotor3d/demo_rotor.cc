/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#include <stdlib.h>
#include <ios>

#include <memory>
#include <vector>
#include <fstream>
#include <ostream>
#include <iterator>
#include <algorithm>
#include <string>
#include <sstream>
#include <cmath>
#include <limits>

#ifdef ENABLE_MPI 
  #include "mpi.h"
#endif

#include "wonton/mesh/simple/simple_mesh.h"
#include "wonton/mesh/simple/simple_mesh_wrapper.h"

#include "tangram/support/tangram.h"
#include "tangram/driver/driver.h"
#include "tangram/reconstruct/VOF.h"
#include "tangram/reconstruct/MOF.h"
#include "tangram/driver/write_to_gmv.h"

#include "tangram/utility/rpgtools/cuts.h"
#include "tangram/utility/rpgtools/primitives.h"
#include "tangram/utility/get_mat_sym_diff_vol.h"
#include "tangram/utility/rpgtools/examples/matdata_rotor3d.h"

/* Demo app for the 3D Rotor example.
   Uses SimpleMesh, VOF, and MOF.
   Generates mesh, computes material moments and reference poly's
   for the 3D Rotor example, performs interface reconstruction
   with both VOF and MOF methods, computes reconstruction errors
   based on volumes of symmetric differences with respect to the
   reference poly's, and outputs reconstruction results to gmv files. */

int main(int argc, char** argv) {
#ifdef ENABLE_MPI  
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  int world_size = 1;
  MPI_Comm_size(comm, &world_size);
  if (world_size > 1)
    throw std::runtime_error("This app is designed to run in serial!");
#endif  

  if (argc != 5) {
    std::ostringstream os;
    os << std::endl <<
      "Correct usage: rotor_demo <decompose_cells> [0|1]" <<
      " <nx> <ny> <nz> " << std::endl;
    throw std::runtime_error(os.str());
  }

  bool decompose_cells = atoi(argv[1]);
  std::string mesh_name;
  int nx = 0, ny = 0, nz = 0;

  nx = atoi(argv[2]); ny = atoi(argv[3]); nz = atoi(argv[4]);
  {
    std::ostringstream os;
    os << "box_mesh_" << nx << "x" << ny << "x" << nz;
    mesh_name = os.str();
  }

  std::string ref_gmv_fname = mesh_name + "_ref_matpolys.gmv";

  float seconds_taken = 0.0;
  struct timeval begin_timeval, end_timeval, diff_timeval;

  gettimeofday(&begin_timeval, 0);

  Wonton::Simple_Mesh mesh(0.25, 0.1, 0.25,
                           0.75, 0.9, 0.75,
                           nx, ny, nz);
  Wonton::Simple_Mesh_Wrapper mesh_wrapper(mesh);

  gettimeofday(&end_timeval, 0);
  timersub(&end_timeval, &begin_timeval, &diff_timeval);
  seconds_taken = diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;

  std::cout << "Time taken to generate mesh -> " << 
    seconds_taken << " (s)" << std::endl;

  int ncells = mesh_wrapper.num_owned_cells();

  std::vector<std::string> IR_names = {"VOF", "MOF"};

  // Volume and angle tolerances
  // The choice of the distance tolerance below ensures that all non-identical points
  // have at least one coordinate that is machine epsilon apart. Because distance
  // tolerance corresponds to a geometrical distance, point-equivalence checks have
  // no geometrical bias.  
  double dst_tol = sqrt(3)*std::numeric_limits<double>::epsilon();
  double vol_tol = std::numeric_limits<double>::epsilon();
  //double vol_tol = 1.0e-24;
  std::vector< Tangram::IterativeMethodTolerances_t> ims_tols(2);
  ims_tols[0] = {1000, dst_tol, vol_tol};
  ims_tols[1] = {100, 1.0e-15, 1.0e-12};

  std::vector<int> cell_num_mats;
  std::vector<int> cell_mat_ids;
  std::vector<double> cell_mat_volfracs;
  std::vector<Tangram::Point3> cell_mat_centroids;
  std::vector< std::vector< std::vector<r3d_poly> > > reference_mat_polys;

  std::vector<int> mesh_material_IDs;
  std::vector< std::string > mesh_material_names;

  gettimeofday(&begin_timeval, 0);

  rotor_material_moments(mesh_wrapper, mesh_material_IDs, mesh_material_names, 
    cell_num_mats, cell_mat_ids, cell_mat_volfracs, cell_mat_centroids, 
    vol_tol, dst_tol, decompose_cells, &reference_mat_polys);

  gettimeofday(&end_timeval, 0);
  timersub(&end_timeval, &begin_timeval, &diff_timeval);
  seconds_taken = diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;

  std::cout << "Time taken to generate reference material polyhedra and their moments -> " << 
    seconds_taken << " (s)" << std::endl;

  int nmesh_material_IDs = mesh_material_IDs.size();
  std::vector<int> offsets(ncells, 0);  
  for (int icell = 0; icell < ncells - 1; icell++)
    offsets[icell + 1] = offsets[icell] + cell_num_mats[icell];

  int max_cell_nref_mats = *std::max_element(cell_num_mats.begin(), cell_num_mats.end());
  int nmmcells = 0;
  // Number of cells containing a certain number of materials ([0] for the # of single-material,
  // [1] for the # of two-material, etc.), used for stats output.
  std::vector<int> ncells_with_xmats(max_cell_nref_mats, 0);
  // For every material in the mesh, its volume contained in multi-material cells only.
  // Used to compute per-material errors over all multi-material cells.
  std::vector<double> mmcells_material_volumes(nmesh_material_IDs, 0.0);
  // For every material in the mesh, its volume contained in multi-material cells with
  // a certain number of materials ([0] for volume contained in two-material cell, etc.).
  // Used to compute per-material errors over multi-material cells with a specific number
  // of contained materials.
  std::vector < std::vector<double> > xmat_cells_material_volumes(nmesh_material_IDs);
  for (int imat = 0; imat < nmesh_material_IDs; imat++)
      xmat_cells_material_volumes[imat].resize(max_cell_nref_mats - 1, 0.0);

  for (int icell = 0; icell < ncells; icell++) {
    int ncmats = cell_num_mats[icell];
    ncells_with_xmats[ncmats - 1]++;
    if (ncmats > 1) {
      double cell_volume = mesh_wrapper.cell_volume(icell);
      for (int icmat = 0; icmat < ncmats; icmat++) {
        int cur_mat_id = cell_mat_ids[offsets[icell] + icmat];
        int mesh_matid = std::distance(mesh_material_IDs.begin(),
          std::find(mesh_material_IDs.begin(), 
                    mesh_material_IDs.end(), cur_mat_id));

        double cell_mat_vol = cell_volume*cell_mat_volfracs[offsets[icell] + icmat];
        if (cell_mat_vol >= vol_tol) {
          mmcells_material_volumes[mesh_matid] += cell_mat_vol;
          xmat_cells_material_volumes[mesh_matid][ncmats - 2] += cell_mat_vol;
        }
      }
      nmmcells++;  
    }
  }

  std::cout << std::endl << "Stats for structured " <<
    nx << "x" << ny << "x" << nz << " computational mesh";

  std::cout << ":" << std::endl << "  " << ncells << " cells," <<
    std::endl << "  " << nmesh_material_IDs << " materials," <<
    std::endl << "  " << ncells_with_xmats[0] << " single-material cells," <<
    std::endl << "  " << nmmcells << " multi-material cells:" << std::endl;

  int count_cells_with_xmats = ncells_with_xmats.size();
  for (int inm = 0; inm < count_cells_with_xmats - 1; inm++)
    std::cout << "    Number of cells with " << inm + 2 << " materials -> " <<
      ncells_with_xmats[inm + 1] << std::endl;

  int nreconstructors = static_cast<int>(IR_names.size());
  std::vector< std::vector< std::shared_ptr< Tangram::CellMatPoly<3> > > >
    IR_cellmatpoly_list(nreconstructors);

  //Do we need VOF?
  int ivof = std::distance(IR_names.begin(),
                           std::find(IR_names.begin(), IR_names.end(), "VOF"));
  if (ivof != nreconstructors) {
    // Build the VOF driver
    std::cout << std::endl << IR_names[ivof] << " interface reconstruction method:" << std::endl;
    Tangram::Driver<Tangram::VOF, 3, Wonton::Simple_Mesh_Wrapper, 
                    Tangram::SplitR3D, Tangram::ClipR3D> 
      vof_driver(mesh_wrapper, ims_tols, !decompose_cells);

    vof_driver.set_volume_fractions(cell_num_mats, cell_mat_ids, 
                                  cell_mat_volfracs, cell_mat_centroids);
    vof_driver.reconstruct();    
    IR_cellmatpoly_list[ivof] = vof_driver.cell_matpoly_ptrs();
  }
  //Do we need MOF?
  int imof = std::distance(IR_names.begin(),
                           std::find(IR_names.begin(), IR_names.end(), "MOF"));
  if (imof != nreconstructors) {
    // Build the MOF driver
    std::cout << std::endl << IR_names[imof] << " interface reconstruction method:" << std::endl; 
    Tangram::Driver<Tangram::MOF, 3, Wonton::Simple_Mesh_Wrapper, 
                    Tangram::SplitR3D, Tangram::ClipR3D> 
      mof_driver(mesh_wrapper, ims_tols, !decompose_cells);

    mof_driver.set_volume_fractions(cell_num_mats, cell_mat_ids, 
                                  cell_mat_volfracs, cell_mat_centroids);
    mof_driver.reconstruct();    
    IR_cellmatpoly_list[imof] = mof_driver.cell_matpoly_ptrs();
  }

  for (int iIR = 0; iIR < nreconstructors; iIR++) {
    std::vector<double> total_mat_sym_diff_vol(nmesh_material_IDs, 0.0);

    gettimeofday(&begin_timeval, 0);

    std::vector< std::vector<double> > xmat_cells_mat_sym_diff_vol(nmesh_material_IDs);
    for (int imat = 0; imat < nmesh_material_IDs; imat++)
      xmat_cells_mat_sym_diff_vol[imat].resize(max_cell_nref_mats - 1, 0.0);

    for (int icell = 0; icell < ncells; icell++) {
      int ncmats = cell_num_mats[icell];
      if (ncmats == 1) {
        assert(IR_cellmatpoly_list[iIR][icell] == nullptr);
        continue;
      }

      std::vector<int> cell_ref_mat_ids(cell_mat_ids.begin() + offsets[icell], 
        cell_mat_ids.begin() + offsets[icell] + ncmats);

      std::vector<double> cell_ref_mat_vols(cell_mat_volfracs.begin() + offsets[icell], 
        cell_mat_volfracs.begin() + offsets[icell] + ncmats); 
      double cell_volume = mesh_wrapper.cell_volume(icell);  
      for (int icmat = 0; icmat < ncmats; icmat++)
        cell_ref_mat_vols[icmat] *= cell_volume;

      std::vector<double> cell_mat_sym_diff_vol;
      if (IR_cellmatpoly_list[iIR][icell] != nullptr)
        get_mat_sym_diff_vol(reference_mat_polys[icell], cell_ref_mat_ids, 
                             cell_ref_mat_vols, IR_cellmatpoly_list[iIR][icell], 
                             cell_mat_sym_diff_vol, !decompose_cells);

      for (int icmat = 0; icmat < ncmats; icmat++) {
        int material_id = cell_ref_mat_ids[icmat];
        int mesh_matid = std::distance(mesh_material_IDs.begin(),
          std::find(mesh_material_IDs.begin(), 
                    mesh_material_IDs.end(), material_id));
        
        if (IR_cellmatpoly_list[iIR][icell] != nullptr) {
          total_mat_sym_diff_vol[mesh_matid] += cell_mat_sym_diff_vol[icmat];
          xmat_cells_mat_sym_diff_vol[mesh_matid][ncmats - 2] += cell_mat_sym_diff_vol[icmat];
        }
        else if (cell_ref_mat_vols[icmat] > ims_tols[0].fun_eps) {
          total_mat_sym_diff_vol[mesh_matid] += cell_ref_mat_vols[icmat];
          xmat_cells_mat_sym_diff_vol[mesh_matid][ncmats - 2] += cell_ref_mat_vols[icmat];
        }
      }
    }

    gettimeofday(&end_timeval, 0);
    timersub(&end_timeval, &begin_timeval, &diff_timeval);
    seconds_taken = diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;

    std::cout.unsetf(std::ios_base::floatfield);
    std::cout << "Time taken to compute errors -> " <<
      seconds_taken << " (s)" << std::endl;

    std::cout << std::endl << IR_names[iIR] << " interface reconstruction method:" << std::endl <<
      std::endl;

    for (int imat = 0; imat < nmesh_material_IDs; imat++) {
      std::cout << "  " << mesh_material_names[imat] << ", material ID " << 
        mesh_material_IDs[imat] << " -> " << std::endl << 
        "    over all multi-material cells:" << std::endl << std::scientific <<
        "      Aggregate vol = " << mmcells_material_volumes[imat] << "," << std::endl <<
        "      aggregate sym.diff.vol = " << total_mat_sym_diff_vol[imat];
      if (total_mat_sym_diff_vol[imat] > std::numeric_limits<double>::epsilon())
        std::cout << "," << std::endl << "      relative sym.diff.vol = " << 
          total_mat_sym_diff_vol[imat]/mmcells_material_volumes[imat]; 
      std::cout << std::endl;

      for (int mc = 0; mc < count_cells_with_xmats - 1; mc++) {
        if (xmat_cells_material_volumes[imat][mc] != 0.0) {
          std::cout << "    over all " << mc + 2 << "-material cells:" << std::endl <<
          "      Aggregate vol = " << xmat_cells_material_volumes[imat][mc] << "," << std::endl <<
          "      aggregate sym.diff.vol = " << xmat_cells_mat_sym_diff_vol[imat][mc] << "," << std::endl <<
          "      relative sym.diff.vol = " <<
              xmat_cells_mat_sym_diff_vol[imat][mc]/xmat_cells_material_volumes[imat][mc] << std::endl;
        }
      }

      std::cout << std::endl;      
    }
    std::cout << std::endl;

    std::string out_gmv_fname = mesh_name;
    if (decompose_cells)
      out_gmv_fname += "_decomposed";
    out_gmv_fname += "_" + IR_names[iIR] + "_res_matpolys.gmv";

    gettimeofday(&begin_timeval, 0);

    write_to_gmv(mesh_wrapper, nmesh_material_IDs, cell_num_mats, cell_mat_ids,
                 IR_cellmatpoly_list[iIR], out_gmv_fname);

    gettimeofday(&end_timeval, 0);
    timersub(&end_timeval, &begin_timeval, &diff_timeval);
    seconds_taken = diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;

    std::cout.unsetf(std::ios_base::floatfield);
    std::cout << "Time taken to write CellMatPolys -> " <<
      seconds_taken << " (s)" << std::endl;

    IR_cellmatpoly_list[iIR].clear();
  }

  gettimeofday(&begin_timeval, 0);

  std::cout << "Writing out reference material polyhedra..." << std::endl;
  std::vector< std::shared_ptr< Tangram::CellMatPoly<3> > > ref_matpoly_list(ncells);
  for (int icell = 0; icell < ncells; icell++) {
    if (cell_num_mats[icell] == 1)
      continue;

    ref_matpoly_list[icell] = std::make_shared< Tangram::CellMatPoly<3> >(icell);

    for (int icmat = 0; icmat < cell_num_mats[icell]; icmat++) {
      int nmp = reference_mat_polys[icell][icmat].size();
      for (int imp = 0; imp < nmp; imp++) {
        std::vector< Tangram::MatPoly<3> > cur_matpoly;
        Tangram::r3dpoly_to_matpolys(reference_mat_polys[icell][icmat][imp], 
                                     cur_matpoly, vol_tol, dst_tol);
        assert(cur_matpoly.size() == 1);

        cur_matpoly[0].set_mat_id(cell_mat_ids[offsets[icell] + icmat]);
        ref_matpoly_list[icell]->add_matpoly(cur_matpoly[0]);
      }
    }
  }

  write_to_gmv(mesh_wrapper, nmesh_material_IDs, cell_num_mats, cell_mat_ids,
               ref_matpoly_list, ref_gmv_fname);  

  gettimeofday(&end_timeval, 0);
  timersub(&end_timeval, &begin_timeval, &diff_timeval);
  seconds_taken = diff_timeval.tv_sec + 1.0E-6*diff_timeval.tv_usec;

  std::cout << "Time taken to write reference material polyhedra -> " <<
    seconds_taken << " (s)" << std::endl;

#ifdef ENABLE_MPI 
  MPI_Finalize();
#endif

  return 0;
}
