/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

// Generate volume fractions and centroids for materials in mesh cells
// given a configuration file with different geometric features.
//
//
// Assumption: 2D meshes are in X-Y plane only
//
// Author: Rao Garimella, rao@lanl.gov

#include <sys/time.h>
#include <string>
#include <memory>
#include <vector>
#include <iostream>

#include <mpi.h>


// Jali mesh infrastructure library
// See https://github.com/lanl/jali

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliStateVector.h"
#include "JaliState.h"

// tangram includes
#include "tangram/support/tangram.h"

// wonton includes
#include "wonton/mesh/jali/jali_mesh_wrapper.h"

#include "vfgen.h"


int main(int argc, char *argv[]) {
  int global_nmats;

  struct timeval begin_total {};
  struct timeval begin_core {};
  struct timeval end_total {};
  struct timeval end_core {};
  struct timeval diff_total {};
  struct timeval diff_core {};

  double ptol = 1.0e-12;

  MPI_Init(&argc, &argv);

  if (argc < 3) {
    std::cerr << "usage: " << argv[0] <<
        " meshfilename dim <featurefile.inp>\n";
    exit(-1);
  }

  gettimeofday(&begin_total, 0);

  std::string meshfilename(argv[1]);
  std::size_t pos = meshfilename.find_last_of(".");
  std::string basename(meshfilename.substr(0, pos));

  std::string featfilename;
#ifndef NDEBUG
  int probdim = std::stoi(argv[2]);
#endif
  if (argc > 3)
    featfilename = argv[3];
  else
    featfilename = basename + std::string(".inp");


  // Initialize the mesh from Jali

  std::shared_ptr<Jali::Mesh> mesh;
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities(Jali::Entity_kind::ALL_KIND);

  mf.partitioner(Jali::Partitioner_type::METIS);
  mesh = mf(meshfilename);

  int meshdim = mesh->space_dimension();
#ifndef NDEBUG
  assert(meshdim == probdim);
#endif

  int ncells = mesh->num_cells();

  Wonton::Jali_Mesh_Wrapper mesh_wrapper(*mesh);

  std::string outfilename(basename + std::string(".vf"));
  std::string boutfilename(basename + std::string(".bvf"));

  int ncells_nmats[MAXMATS] = {};

  if (meshdim == 2) {

    std::vector<FEATURE<2>> features;
    read_features<2>(featfilename, &features, &global_nmats);

    Wonton::vector<vfcen_t<2>> vfcen(ncells);
    InFeatureEvaluator<2> feature_functor(features);
    VolfracEvaluator<2, Wonton::Jali_Mesh_Wrapper>
      volfrac_calculator(mesh_wrapper, feature_functor, ptol, global_nmats);

    gettimeofday(&begin_core, 0);

    // Populate the vf array using the volfrac functor
    Wonton::vector<int> counter(ncells);
    for (int i = 0; i < ncells; i++) counter[i] = i;
    Wonton::transform(counter.begin(), counter.end(), vfcen.begin(),
                       volfrac_calculator);

    gettimeofday(&end_core, 0);

    writeAsciiFile<2>(outfilename, vfcen, global_nmats);
    writeBinaryFile<2>(boutfilename, vfcen);


    // Count the cells the with 1, 2 , 3 ... materials
    for (int i = 0; i < ncells; i++) {
      // When Thrust is turned on, vfcen is on the device and vfcen_i
      // is on the host. Thrust does not allow us to modify device
      // memory directly from host, so we have to use a const
      // reference or do an explicit copy

      vfcen_t<2> const& vfcen_i = vfcen[i];

      int nmats_i = vfcen_i.nmats;
      ncells_nmats[nmats_i]++;
    }

  }  else if (meshdim == 3) {

    std::vector<FEATURE<3>> features;
    read_features<3>(featfilename, &features, &global_nmats);

    Wonton::vector<vfcen_t<3>> vfcen(ncells);
    InFeatureEvaluator<3> feature_functor(features);
    VolfracEvaluator<3, Wonton::Jali_Mesh_Wrapper>
      volfrac_calculator(mesh_wrapper, feature_functor, ptol, global_nmats);

    gettimeofday(&begin_core, 0);

    // Populate the vf array using the volfrac functor
    Wonton::vector<int> counter(ncells);
    for (int i = 0; i < ncells; i++) counter[i] = i;
    Wonton::transform(counter.begin(), counter.end(), vfcen.begin(),
                       volfrac_calculator);

    gettimeofday(&end_core, 0);

    writeAsciiFile<3>(outfilename, vfcen, global_nmats);
    writeBinaryFile<3>(boutfilename, vfcen);

    // Count the cells the with 1, 2 , 3 ... materials
    for (int i = 0; i < ncells; i++) {
      // When Thrust is turned on, vfcen is on the device and vfcen_i
      // is on the host. Thrust does not allow us to modify device
      // memory directly from host, so we have to use a const
      // reference or do an explicit copy

      vfcen_t<3> const& vfcen_i = vfcen[i];

      int nmats_i = vfcen_i.nmats;
      ncells_nmats[nmats_i]++;
    }
  }

  int tot = 0;
  for (int i = 0; i < 5; i++) {
    std::cerr << "Number of cells with " << i << " materials " <<
        ncells_nmats[i] << "\n";
    tot += ncells_nmats[i];
  }
  if (tot < ncells)
    std::cerr <<
        "There may be cells with no materials or more than 5 materials\n";

  MPI_Finalize();

  gettimeofday(&end_total, 0);

  timersub(&end_core, &begin_core, &diff_core);
  timersub(&end_total, &begin_total, &diff_total);
  double core_seconds = diff_core.tv_sec + 1.0E-6*diff_core.tv_usec;
  double total_seconds = diff_total.tv_sec + 1.0E-6*diff_total.tv_usec;

  std::cerr << "Time taken to compute volume fractions and centroids: " <<
      core_seconds << "\n";
  std::cerr << "Total time taken by app: " << total_seconds << "\n";

  return 1;
}
