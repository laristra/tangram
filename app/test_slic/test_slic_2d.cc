/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#include <cstdlib>

#include <memory>
#include <vector>
#include <fstream>
#include <ostream>
#include <iterator>
#include <algorithm>
#include <string>
#include <sstream>

#include "tangram/support/tangram.h"  // WONTON_ENABLE_MPI defined in this file

#ifdef WONTON_ENABLE_MPI
  #include "mpi.h"
#endif
#if defined(WONTON_ENABLE_Jali) && defined(WONTON_ENABLE_MPI)
  #include "Mesh.hh"
  #include "MeshFactory.hh"
  #include "wonton/mesh/jali/jali_mesh_wrapper.h"
#else
  #include "wonton/mesh/simple/simple_mesh.h"
  #include "wonton/mesh/simple/simple_mesh_wrapper.h"
#endif

// tangram includes
#include "tangram/driver/driver.h"
#include "tangram/reconstruct/SLIC.h"
#include "tangram/driver/write_to_gmv.h"
#include "tangram/utility/get_material_moments.h"
#include "tangram/utility/get_mat_sym_diff_vol.h"

// wonton includes
#include "wonton/support/Vector.h"

/* Test app for a 2D mesh and planar material interfaces.
   Uses SimpleMesh/Jali and SLIC.
   Generates (SimpleMesh) or reads mesh from file (Jali),
   computes material moments for a sequence of planar interfaces,
   performs interface reconstruction (SLIC), and outputs volumes of
   symmetric difference for every material in every
   multi-material cell */

#include <set>

const std::vector<int> mesh_materials = {5, 0, 3};
const std::vector< Tangram::Vector2 > material_interface_normals = {
  Tangram::Vector2(-0.5, 0.0), Tangram::Vector2(0, -0.5)
};
const std::vector< Tangram::Point2 > material_interface_points = {
  Tangram::Point2(0.5, 0.0), Tangram::Point2(1.0, 0.5)
};

int main(int argc, char** argv) {
#ifdef WONTON_ENABLE_MPI
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  int world_size = 1;
  MPI_Comm_size(comm, &world_size);
  if (world_size > 1)
    throw std::runtime_error("This app is designed to run in serial!");
#endif

  assert((material_interface_normals.size() == material_interface_points.size()) &&
         (mesh_materials.size() == material_interface_normals.size() + 1));

#if defined(WONTON_ENABLE_Jali) && defined(WONTON_ENABLE_MPI)
  if (argc != 3) {
    std::ostringstream os;
    os << std::endl <<
    "Correct usage: test_slic_2d" << " <decompose_cells> [0|1]" <<
      " <base_mesh_file>" << std::endl;
    throw std::runtime_error(os.str());
  }
#else
  if (argc != 4) {
    std::ostringstream os;
    os << std::endl <<
    "Correct usage: test_slic_2d" << " <decompose_cells> [0|1]" <<
      " <nx> <ny> " << std::endl;
    throw std::runtime_error(os.str());
  }
#endif

  bool decompose_cells = atoi(argv[1]);

  int nmesh_materials = static_cast<int>(mesh_materials.size());
  std::vector< Tangram::Plane_t<2> > material_interfaces(nmesh_materials - 1);
  for (int iplane = 0; iplane < nmesh_materials - 1; iplane++) {
    material_interfaces[iplane].normal = material_interface_normals[iplane];
    material_interfaces[iplane].normal.normalize();
    material_interfaces[iplane].dist2origin =
      -Wonton::dot(material_interface_points[iplane].asV(),
                    material_interfaces[iplane].normal);
  }

#if defined(WONTON_ENABLE_Jali) && defined(WONTON_ENABLE_MPI)
  std::string mesh_name = argv[2];
  mesh_name.resize(mesh_name.size() - 4);
#else
  int nx = atoi(argv[2]), ny = atoi(argv[3]);
  std::string mesh_name;
  {
    std::ostringstream os;
    os << "simple_mesh_" << nx << "x" << ny ;
    mesh_name = os.str();
  }
#endif

#ifndef NDEBUG
  std::string ref_gmv_fname = mesh_name + "_ref_matpolys.gmv";
  std::string out_gmv_fname = mesh_name + "_res_matpolys.gmv";
#endif

#if defined(WONTON_ENABLE_Jali) && defined(WONTON_ENABLE_MPI)
  Jali::MeshFactory mesh_factory(comm);
  mesh_factory.framework(Jali::MSTK);
  mesh_factory.included_entities(Jali::Entity_kind::EDGE);
  std::shared_ptr<Jali::Mesh> mesh = mesh_factory(argv[2]);

  assert(mesh != nullptr);
  Wonton::Jali_Mesh_Wrapper mesh_wrapper(*mesh, true, false, false);
#else
  Wonton::Simple_Mesh mesh(0.0, 0.0,
                           1.0, 1.0,
                           nx, ny );
  Wonton::Simple_Mesh_Wrapper mesh_wrapper(mesh);
#endif

  int ncells = mesh_wrapper.num_owned_cells();

  // Volume and distance tolerance
  double dst_tol = sqrt(2)*std::numeric_limits<double>::epsilon();
  double vol_tol = std::numeric_limits<double>::epsilon();
  std::vector< Tangram::IterativeMethodTolerances_t> ims_tols(1);
  ims_tols[0] = {1000, dst_tol, vol_tol}; 

  std::vector<int> cell_num_mats;
  std::vector<int> cell_mat_ids;
  std::vector<double> cell_mat_volfracs;
  std::vector<Tangram::Point2> cell_mat_centroids;
  std::vector< std::vector< std::vector<r2d_poly> > > reference_mat_polys;

#if defined(WONTON_ENABLE_Jali) && defined(WONTON_ENABLE_MPI)
  get_material_moments<Wonton::Jali_Mesh_Wrapper>(mesh_wrapper, material_interfaces,
    mesh_materials, cell_num_mats, cell_mat_ids, cell_mat_volfracs, cell_mat_centroids,
    vol_tol, dst_tol, decompose_cells, &reference_mat_polys);
#else
  get_material_moments<Wonton::Simple_Mesh_Wrapper>(mesh_wrapper, material_interfaces,
    mesh_materials, cell_num_mats, cell_mat_ids, cell_mat_volfracs, cell_mat_centroids,
    vol_tol, dst_tol, decompose_cells, &reference_mat_polys);
#endif

  std::vector<int> offsets(ncells, 0);
  for (int icell = 0; icell < ncells - 1; icell++)
    offsets[icell + 1] = offsets[icell] + cell_num_mats[icell];

  int nmmcells = 0;
  std::vector<double> mmcells_material_volumes(nmesh_materials, 0.0);
  for (int icell = 0; icell < ncells; icell++)
    if (cell_num_mats[icell] > 1) {
      double cell_volume = mesh_wrapper.cell_volume(icell);
      for (int icmat = 0; icmat < cell_num_mats[icell]; icmat++) {
        int cur_mat_id = cell_mat_ids[offsets[icell] + icmat];
        int mesh_matid = std::distance(mesh_materials.begin(),
          std::find(mesh_materials.begin(),
                   mesh_materials.end(), cur_mat_id));
        mmcells_material_volumes[mesh_matid] +=
          cell_volume*cell_mat_volfracs[offsets[icell] + icmat];
      }
      nmmcells++;
    }

#ifndef NDEBUG
  std::vector< std::shared_ptr< Tangram::CellMatPoly<2> > > ref_matpoly_list(ncells);
  for (int icell = 0; icell < ncells; icell++) {
    ref_matpoly_list[icell] = std::make_shared< Tangram::CellMatPoly<2> >(icell);

    for (int icmat = 0; icmat < cell_num_mats[icell]; icmat++) {
      int nmp = static_cast<int>(reference_mat_polys[icell][icmat].size());
      for (int imp = 0; imp < nmp; imp++) {
        r2d_print(&reference_mat_polys[icell][icmat][imp]);

        Tangram::MatPoly<2> cur_matpoly;
        Tangram::r2dpoly_to_matpoly(reference_mat_polys[icell][icmat][imp], 
                                    cur_matpoly, dst_tol);
        cur_matpoly.set_mat_id(cell_mat_ids[offsets[icell] + icmat]);

        std::vector<Tangram::Point2> mt_pts = cur_matpoly.points();
        int const nverts = mt_pts.size();
        std::cout<<"nverts = "<< nverts << std::endl;
        for (int i = 0; i < nverts; i++)
          std::cout<<"Pt["<<i<<"] = { "<<mt_pts[i][0]<<", "<<mt_pts[i][1]<<"}"<<std::endl;
        

        ref_matpoly_list[icell]->add_matpoly(cur_matpoly);
      }
    }
  }

  write_to_gmv(ref_matpoly_list, ref_gmv_fname);
#endif

  // Build the driver
#if defined(WONTON_ENABLE_Jali) && defined(WONTON_ENABLE_MPI)
  Tangram::Driver<Tangram::SLIC, 2, Wonton::Jali_Mesh_Wrapper,
                  Tangram::SplitR2D, Tangram::ClipR2D>
    slic_driver(mesh_wrapper, ims_tols, !decompose_cells);
#else
  Tangram::Driver<Tangram::SLIC, 2, Wonton::Simple_Mesh_Wrapper,
                  Tangram::SplitR2D, Tangram::ClipR2D>
    slic_driver(mesh_wrapper, ims_tols, !decompose_cells);
#endif

  slic_driver.set_volume_fractions(cell_num_mats, cell_mat_ids,
                                  cell_mat_volfracs, cell_mat_centroids);
  slic_driver.reconstruct();

  std::vector<std::shared_ptr<Tangram::CellMatPoly<2>>> cellmatpoly_list =
    slic_driver.cell_matpoly_ptrs();

  std::vector<double> total_mat_sym_diff_vol(nmesh_materials, 0.0);
  std::vector<double> max_mat_sym_diff_vol(nmesh_materials, 0.0);
  std::vector<int> max_mat_sym_diff_icell(nmesh_materials, -1);
  std::vector< std::vector<double> > mat_cells_sym_diff_vol(nmesh_materials);
  for (int imat = 0; imat < nmesh_materials; imat++)
    mat_cells_sym_diff_vol[imat].resize(ncells, -1.0);

  for (int icell = 0; icell < ncells; icell++) {
    int ncmats = cell_num_mats[icell];
    if (ncmats == 1) {
      assert(cellmatpoly_list[icell] == nullptr);
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
    get_mat_sym_diff_vol(reference_mat_polys[icell], cell_ref_mat_ids,
                         cell_ref_mat_vols, cellmatpoly_list[icell],
                         cell_mat_sym_diff_vol, !decompose_cells);

    for (int icmat = 0; icmat < ncmats; icmat++) {
      int material_id = cell_ref_mat_ids[icmat];
      int mesh_matid = std::distance(mesh_materials.begin(),
        std::find(mesh_materials.begin(),
                  mesh_materials.end(), material_id));

      total_mat_sym_diff_vol[mesh_matid] += cell_mat_sym_diff_vol[icmat];
      if (cell_mat_sym_diff_vol[icmat] > max_mat_sym_diff_vol[mesh_matid]) {
        max_mat_sym_diff_vol[mesh_matid] = cell_mat_sym_diff_vol[icmat];
        max_mat_sym_diff_icell[mesh_matid] = icell;
      }

      mat_cells_sym_diff_vol[mesh_matid][icell] = cell_mat_sym_diff_vol[icmat];
    }
  }

  std::string res_out_fname = "cell_sym_diff_2d_" + mesh_name;
  if (decompose_cells)
    res_out_fname += "_decomposed";
  res_out_fname += ".txt";

  std::ofstream fout(res_out_fname);
  fout << std::scientific;
  fout.precision(17);
  for (int imat = 0; imat < nmesh_materials; imat++)
    for (int icell = 0; icell < ncells; icell++)
      if (mat_cells_sym_diff_vol[imat][icell] != -1.0)
        fout << mesh_materials[imat]*ncells + icell << " " <<
          mat_cells_sym_diff_vol[imat][icell] << std::endl;

  fout.close();

std::cout << std::endl << "Stats for ";
#if defined(WONTON_ENABLE_Jali) && defined(WONTON_ENABLE_MPI)
  std::cout << "computational mesh " << mesh_name;
#else
  std::cout << "structured " << nx << "x" << ny <<
    " computational mesh";
#endif
  std::cout << ":" << std::endl << "  " << ncells << " cells," <<
    std::endl << "  " << nmesh_materials << " materials," <<
    std::endl << "  " << nmmcells << " multi-material cells." <<
    std::endl << std::endl <<
    "For each material over all multi-material cells:" << std::endl;

  for (int imat = 0; imat < nmesh_materials; imat++) {

    int imaxcell = 0;
    double max_sym_diff_mat_vol = 0.;

    if (max_mat_sym_diff_icell[imat] != -1) {
      imaxcell = max_mat_sym_diff_icell[imat];
      int cell_matid =
        std::distance(cell_mat_ids.begin() + offsets[imaxcell],
          std::find(cell_mat_ids.begin() + offsets[imaxcell],
                    cell_mat_ids.begin() + offsets[imaxcell] + cell_num_mats[imaxcell],
                    mesh_materials[imat]));

      max_sym_diff_mat_vol = cell_mat_volfracs[offsets[imaxcell] + cell_matid]*
        mesh_wrapper.cell_volume(imaxcell);
    }

    std::cout << "  Material ID " << mesh_materials[imat] << " -> " << std::endl <<
      "    Aggregate vol = " << mmcells_material_volumes[imat] << "," << std::endl <<
      "    aggregate sym.diff.vol = " << total_mat_sym_diff_vol[imat];

    if (total_mat_sym_diff_vol[imat] >= vol_tol)
      std::cout << "," << std::endl << "    relative sym.diff.vol = " <<
        total_mat_sym_diff_vol[imat]/mmcells_material_volumes[imat];

    std::cout << std::endl;

    if (max_mat_sym_diff_icell[imat] != -1)
      std::cout << "    Max sym.diff.vol in cell " << imaxcell << ":" << std::endl <<
      "      cell material vol = " << max_sym_diff_mat_vol << "," << std::endl <<
      "      sym.diff.vol = " << max_mat_sym_diff_vol[imat] << "," << std::endl <<
      "      relative sym.diff.vol = " <<
        max_mat_sym_diff_vol[imat]/max_sym_diff_mat_vol << std::endl;
  }

#ifndef NDEBUG
  //Create MatPoly's for single-material cells
  for (int icell = 0; icell < ncells; icell++) {
    if (cell_num_mats[icell] == 1) {
      assert(cellmatpoly_list[icell] == nullptr);
      std::shared_ptr< Tangram::CellMatPoly<2> >
        cmp_ptr(new Tangram::CellMatPoly<2>(icell));
      // Create a MatPoly with the same geometry as the current cell        
      Tangram::MatPoly<2> cell_matpoly;
      cell_get_matpoly(mesh_wrapper, icell, &cell_matpoly, dst_tol);
      cell_matpoly.set_mat_id(cell_mat_ids[offsets[icell]]);
      cmp_ptr->add_matpoly(cell_matpoly);
      cellmatpoly_list[icell] = cmp_ptr;
    }
  }

  write_to_gmv(cellmatpoly_list, out_gmv_fname);
#endif

#ifdef WONTON_ENABLE_MPI
  MPI_Finalize();
#endif

  return 0;
}
