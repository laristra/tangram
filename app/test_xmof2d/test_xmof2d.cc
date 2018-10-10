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

#ifdef ENABLE_MPI 
  #include "mpi.h"
#endif
#if ENABLE_JALI
  #include "Mesh.hh"
  #include "MeshFactory.hh"
  #include "tangram/wrappers/mesh/jali/jali_mesh_wrapper.h"
#else
  #include "tangram/simple_mesh/simple_mesh.h"
  #include "tangram/wrappers/mesh/simple_mesh/simple_mesh_wrapper.h"
#endif

#include "tangram/support/tangram.h"
#include "tangram/driver/driver.h"
#include "tangram/reconstruct/xmof2D_wrapper.h"
#include "tangram/driver/write_to_gmv.h"
#include "tangram/utility/get_material_moments.h"
#include "tangram/utility/get_mat_sym_diff_vol.h"

/* Test app for a 2D mesh and planar material interfaces.
   Uses SimpleMesh/Jali and XMOF2D.
   Generates (SimpleMesh) or reads mesh from file (Jali), 
   computes material moments for a sequence of planar interfaces,
   performs interface reconstruction (XMOF2D), and outputs volumes of
   symmetric difference for every material in every
   multi-material cell */

#include <set>

const std::vector<int> mesh_materials = {5, 0, 3}; 
const std::vector< Tangram::Vector2 > material_interface_normals = {
  Tangram::Vector2(0.5, 0.5), Tangram::Vector2(0.5, -0.375)
};
const std::vector< Tangram::Point2 > material_interface_points = {
  Tangram::Point2(0.5, 0.5), Tangram::Point2(0.5, 0.5)
};



int main(int argc, char** argv) {
#ifdef ENABLE_MPI  
  MPI_Init(&argc, &argv);
  MPI_Comm comm = MPI_COMM_WORLD;

  int world_size = 1;
  MPI_Comm_size(comm, &world_size);
  if (world_size > 1)
    throw std::runtime_error("This app is designed to run in serial!");
#endif  

  assert((material_interface_normals.size() == material_interface_points.size()) &&
         (mesh_materials.size() == material_interface_normals.size() + 1));

#ifdef ENABLE_JALI
  if (argc != 2) {
    std::ostringstream os;
    os << std::endl <<
    "Correct usage: test_xmof2d <base_mesh_file>" << std::endl;
    throw std::runtime_error(os.str());
  }
#else
  if (argc != 3) {
    std::ostringstream os;
    os << std::endl <<
    "Correct usage: test_xmof2d <nx> <ny>" << std::endl;
    throw std::runtime_error(os.str());
  }
#endif  

  bool decompose_cells = false;

  int nmesh_materials = static_cast<int>(mesh_materials.size()); 
  std::vector< Tangram::Plane_t<2> > material_interfaces(nmesh_materials - 1);
  for (int iline = 0; iline < nmesh_materials - 1; iline++) {
    material_interfaces[iline].normal = material_interface_normals[iline];
    material_interfaces[iline].normal.normalize();
    material_interfaces[iline].dist2origin = 
      -Tangram::dot(material_interface_points[iline].asV(),
                    material_interfaces[iline].normal);
  }

#ifdef ENABLE_JALI
  std::string mesh_name = argv[1];
  mesh_name.resize(mesh_name.size() - 4);
#else
  int nx = atoi(argv[1]), ny = atoi(argv[2]);
  std::string mesh_name;
  {
    std::ostringstream os;
    os << "simple_mesh_" << nx << "x" << ny;
    mesh_name = os.str();
  }
#endif  

#ifdef DEBUG
  std::string ref_gmv_fname = mesh_name + "_ref_matpolys.gmv";
  std::string out_gmv_fname = mesh_name + "_res_matpolys.gmv";
#endif

#ifdef ENABLE_JALI
  Jali::MeshFactory mesh_factory(comm);
  mesh_factory.framework(Jali::MSTK);
  mesh_factory.included_entities({Jali::Entity_kind::EDGE, Jali::Entity_kind::FACE});
  std::shared_ptr<Jali::Mesh> mesh = mesh_factory(argv[1]);
  
  assert(mesh != nullptr);
  Tangram::Jali_Mesh_Wrapper mesh_wrapper(*mesh, true, false, false);
#else
  Tangram::Simple_Mesh mesh(0.0, 0.0,
                            1.0, 1.0,
                            nx, ny);
  Tangram::Simple_Mesh_Wrapper mesh_wrapper(mesh);
#endif

  int ncells = mesh_wrapper.num_owned_cells();

  std::vector<int> cell_num_mats;
  std::vector<int> cell_mat_ids;
  std::vector<double> cell_mat_volfracs;
  std::vector<Tangram::Point2> cell_mat_centroids;
  std::vector< std::vector< std::vector<r2d_poly> > > reference_mat_polys;

#ifdef ENABLE_JALI
  get_material_moments<Tangram::Jali_Mesh_Wrapper>(mesh_wrapper, material_interfaces, 
    mesh_materials, cell_num_mats, cell_mat_ids, cell_mat_volfracs, cell_mat_centroids,
    reference_mat_polys, decompose_cells);
#else
  get_material_moments<Tangram::Simple_Mesh_Wrapper>(mesh_wrapper, material_interfaces, 
    mesh_materials, cell_num_mats, cell_mat_ids, cell_mat_volfracs, cell_mat_centroids,
    reference_mat_polys, decompose_cells);
#endif    

  std::vector<int> offsets(ncells, 0);  
  for (int icell = 0; icell < ncells - 1; icell++)
    offsets[icell + 1] = offsets[icell] + cell_num_mats[icell];

  std::vector<double> material_volumes(nmesh_materials, 0.0);
  for (int icell = 0; icell < ncells; icell++) {
    double cell_volume = mesh_wrapper.cell_volume(icell);
    for (int icmat = 0; icmat < cell_num_mats[icell]; icmat++) {
      int cur_mat_id = cell_mat_ids[offsets[icell] + icmat];
      int mesh_matid = std::distance(mesh_materials.begin(),
        std::find(mesh_materials.begin(), 
                  mesh_materials.end(), cur_mat_id));

      material_volumes[mesh_matid] += 
        cell_volume*cell_mat_volfracs[offsets[icell] + icmat];
    }
  }

#ifdef DEBUG
  std::vector< std::shared_ptr< Tangram::CellMatPoly<2> > > ref_matpoly_list(ncells);
  for (int icell = 0; icell < ncells; icell++) {
    ref_matpoly_list[icell] = std::make_shared< Tangram::CellMatPoly<2> >(icell);

    for (int icmat = 0; icmat < cell_num_mats[icell]; icmat++) {
      int nmp = static_cast<int>(reference_mat_polys[icell][icmat].size());
      for (int imp = 0; imp < nmp; imp++) {
        
        const r2d_poly& cur_r2d_poly = reference_mat_polys[icell][icmat][imp];
        std::vector<Tangram::Point2> poly_vrts;
        int ir2dvrt = 0;
        for (int ivrt = 0; ivrt < cur_r2d_poly.nverts; ivrt++) {
          poly_vrts.push_back(Tangram::Point2(cur_r2d_poly.verts[ir2dvrt].pos.xy[0],
                                              cur_r2d_poly.verts[ir2dvrt].pos.xy[1]));
          ir2dvrt = cur_r2d_poly.verts[ir2dvrt].pnbrs[0];                           
        }

        Tangram::MatPoly<2> cur_matpoly(cell_mat_ids[offsets[icell] + icmat]);
        cur_matpoly.initialize(poly_vrts);

        ref_matpoly_list[icell]->add_matpoly(cur_matpoly);
      }
    }
  }

  write_to_gmv(ref_matpoly_list, ref_gmv_fname);
#endif

  // Volume and angles tolerance
  std::vector<Tangram::IterativeMethodTolerances_t> ims_tols(2);
  ims_tols[0] = {.max_num_iter = 1000, .arg_eps = 1.0e-15, .fun_eps = 1.0e-15};
  ims_tols[1] = {.max_num_iter = 1000, .arg_eps = 1.0e-14, .fun_eps = 1.0e-15};

  // Build the driver
#ifdef ENABLE_JALI  
  Tangram::Driver<Tangram::XMOF2D_Wrapper, 2, Tangram::Jali_Mesh_Wrapper> 
    xmof_driver(mesh_wrapper, ims_tols, !decompose_cells);
#else
  Tangram::Driver<Tangram::XMOF2D_Wrapper, 2, Tangram::Simple_Mesh_Wrapper> 
    xmof_driver(mesh_wrapper, ims_tols, !decompose_cells);
#endif    

  xmof_driver.set_volume_fractions(cell_num_mats, cell_mat_ids, 
                                   cell_mat_volfracs, cell_mat_centroids);
  xmof_driver.reconstruct();    

  std::vector<std::shared_ptr<Tangram::CellMatPoly<2>>> cellmatpoly_list = 
    xmof_driver.cell_matpoly_ptrs();

  std::vector<double> total_mat_sym_diff_vol(nmesh_materials, 0.0);
  std::vector<double> max_mat_sym_diff_vol(nmesh_materials, 0.0);
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
                         cell_mat_sym_diff_vol, true);

    for (int icmat = 0; icmat < ncmats; icmat++) {
      int material_id = cell_ref_mat_ids[icmat];
      int mesh_matid = std::distance(mesh_materials.begin(),
        std::find(mesh_materials.begin(), 
                  mesh_materials.end(), material_id));
                  
      total_mat_sym_diff_vol[mesh_matid] += cell_mat_sym_diff_vol[icmat];
      if (cell_mat_sym_diff_vol[icmat] > max_mat_sym_diff_vol[mesh_matid])
        max_mat_sym_diff_vol[mesh_matid] = cell_mat_sym_diff_vol[icmat];

      mat_cells_sym_diff_vol[mesh_matid][icell] = cell_mat_sym_diff_vol[icmat];
    }
  }

  std::string res_out_fname = "cell_sym_diff_" + mesh_name;
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

  std::cout << std::endl << 
    "Symmetric difference over mesh: " << std::endl;
  for (int imat = 0; imat < nmesh_materials; imat++)
    std::cout << "Material #" << mesh_materials[imat] << " -> Total vol = " <<
      material_volumes[imat] << ", Total sym.diff.vol = " <<
      total_mat_sym_diff_vol[imat] << ", Max cell sym.diff.vol = " <<
      max_mat_sym_diff_vol[imat] << std::endl;

#ifdef DEBUG
  //Create MatPoly's for single-material cells
  for (int icell = 0; icell < ncells; icell++) {
    if (cell_num_mats[icell] == 1) {
      assert(cellmatpoly_list[icell] == nullptr);
      std::shared_ptr< Tangram::CellMatPoly<2> > 
        cmp_ptr(new Tangram::CellMatPoly<2>(icell));
      Tangram::MatPoly<2> cell_matpoly;
      cell_get_matpoly(mesh_wrapper, icell, &cell_matpoly);
      cell_matpoly.set_mat_id(cell_mat_ids[offsets[icell]]);
      cmp_ptr->add_matpoly(cell_matpoly);
      cellmatpoly_list[icell] = cmp_ptr;
    }
  }

  write_to_gmv(cellmatpoly_list, out_gmv_fname);
#endif

#ifdef ENABLE_MPI
  MPI_Finalize();
#endif  

  return 0;
}
