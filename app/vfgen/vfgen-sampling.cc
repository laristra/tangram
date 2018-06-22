// Generate volume fractions and centroids for materials in mesh cells
// given a configuration file with different geometric features.
//
//
// Assumption: 2D meshes are in X-Y plane only
//
// Author: Rao Garimella, rao@lanl.gov
// Copyright Los Alamos National Laboratory, 2016


int global_nmats;

#include <sys/time.h>

#include <mpi.h>


#include "tangram/support/tangram.h"


// Jali mesh infrastructure library
// See https://github.com/lanl/jali

#include "Mesh.hh"
#include "MeshFactory.hh"
#include "JaliStateVector.h"
#include "JaliState.h"

#include "tangram/support/tangram.h"
#include "tangram/support/Point.h"
#include "tangram/support/Vector.h"
#include "tangram/wrappers/mesh/jali/jali_mesh_wrapper.h"

#include "vfgen.h"


// Reads the feature file
template <int dim>
void read_features(std::string featfilename,
                   std::vector<FEATURE<dim>> *features) {

  features->clear();

  std::ifstream featfile;
  featfile.open(featfilename.c_str());
  if (!featfile.is_open()) {
    std::cerr << "Could not open file " << featfilename << "\n";
    exit(-98);
  }

  std::string temp_str;
  featfile >> temp_str;
  if (temp_str == "nmats")
    featfile >> global_nmats;
  else
    std::cerr << featfilename << ": first line must be 'nmats n' " <<
        "where n is number of materials\n";


  std::string inout_str, feat_str;
  int nfeat = 0;
  while (!featfile.eof()) {
    featfile >> feat_str;

    if (feat_str == "fill") {

      FEATURE<dim> this_feature;
      this_feature.type = FEATURETYPE::FILL;
      featfile >> this_feature.matid;
      features->push_back(this_feature);

    } else if (feat_str == "halfspace") {

      FEATURE<dim> this_feature;
      this_feature.type = FEATURETYPE::HALFSPACE;
      featfile >> this_feature.matid;
      for (int i = 0; i < dim; i++)
        featfile >> this_feature.plane_xyz[i];
      for (int i = 0; i < dim; i++)
        featfile >> this_feature.plane_normal[i];
      features->push_back(this_feature);

    }  else if (feat_str == "box") {

      FEATURE<dim> this_feature;
      this_feature.type = FEATURETYPE::BOX;
      featfile >> inout_str;
      if (inout_str == "in")
        this_feature.inout = (inout_str == "in") ? 1 : 0;
      featfile >> this_feature.matid;
      for (int i = 0; i < dim; i++)
        featfile >> this_feature.minxyz[i];
      for (int i = 0; i < dim; i++)
        featfile >> this_feature.maxxyz[i];
      features->push_back(this_feature);

    } else if (feat_str == "polytope" || feat_str == "polygon") {

      FEATURE<dim> this_feature;
      this_feature.type = FEATURETYPE::POLY;
      featfile >> inout_str;
      this_feature.inout = (inout_str == "in") ? 1 : 0;
      featfile >> this_feature.matid;

      if (dim == 3) {
        std::cerr << "Cannot handle 3D polytopes yet\n";
        // WHEN WE DO, WE HAVE TO CONVERT THE 3D POLYTOPES INTO POLYHEDRA
        // WITH A TRIANGULATED BOUNDARY
        exit(-1);
      } else {

        featfile >> this_feature.nppoly;
        for (int j = 0; j < this_feature.nppoly; j++) {
          double polyx, polyy;
          featfile >> this_feature.polyxyz[j][0];
          featfile >> this_feature.polyxyz[j][1];
        }
      }
      features->push_back(this_feature);

    } else if (feat_str == "sphere" || feat_str == "circle") {

      FEATURE<dim> this_feature;
      this_feature.type = FEATURETYPE::SPHERE;
      featfile >> inout_str;
      this_feature.inout = (inout_str == "in") ? 1 : 0;
      featfile >> this_feature.matid;
      for (int i = 0; i < dim; i++)
        featfile >> this_feature.cen[i];
      featfile >> this_feature.radius;
      features->push_back(this_feature);

    } else if (feat_str[0] == '#') {
      /* comment - do nothing */
      continue;
    } else if (feat_str == "end") {
      break;
    } else {
      std::cerr << "Unrecognized keyword: " << feat_str << "\n";
      continue;
    }
  }
  nfeat = features->size();

  featfile.close();
}  // read_features



int main(int argc, char *argv[]) {
  char mstk_fname[256], mname[256], feat_fname[256], out_fname[256];
  int i, j, t, ok, status, idx, idx2, nrv, len, flag, done, ir;
  int imat, ptin, inout, nppoly, mkid, rid, refmk, ptnum, nblock;
  double xmin, xmax, ymin, ymax, xlen, ylen;
  double XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX, XLEN, YLEN, ZLEN;
  double radius, radius2, area;
  double cen[3];
  char temp_str[256], feat_str[256], inout_str[8];

  struct timeval begin_total, begin_core, end_total, end_core,
      diff_total, diff_core;

  double ptol = 1.0e-12;
  double ptol2 = ptol*ptol;

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
  int probdim = std::stoi(argv[2]);
  if (argc > 3)
    featfilename = argv[3];
  else
    featfilename = basename + std::string(".inp");


  // Initialize the mesh from Jali

  std::shared_ptr<Jali::Mesh> mesh;
  Jali::MeshFactory mf(MPI_COMM_WORLD);
  mf.included_entities({Jali::Entity_kind::ALL_KIND});

  mf.partitioner(Jali::Partitioner_type::METIS);
  mesh = mf(meshfilename);

  int meshdim = mesh->space_dimension();
  assert(meshdim == probdim);
  
  int ncells = mesh->num_cells();

  Tangram::Jali_Mesh_Wrapper mesh_wrapper(*mesh);

  std::string outfilename(basename + std::string(".vf"));
  std::string boutfilename(basename + std::string(".bvf"));

  int ncells_nmats[MAXMATS] = {};

  if (meshdim == 2) {

    std::vector<FEATURE<2>> features;
    read_features<2>(featfilename, &features);

    Tangram::vector<vfcen_t<2>> vfcen(ncells);
    InFeatureEvaluator<2> feature_functor(features);
    VolfracEvaluator<2, Tangram::Jali_Mesh_Wrapper>
        volfrac_calculator(mesh_wrapper, feature_functor, ptol);

    gettimeofday(&begin_core, 0);

    // Populate the vf array using the volfrac functor
    Tangram::vector<int> counter(ncells);
    for (int i = 0; i < ncells; i++) counter[i] = i;
    Tangram::transform(counter.begin(), counter.end(), vfcen.begin(),
                       volfrac_calculator);

    gettimeofday(&end_core, 0);

    writeAsciiFile<2>(outfilename, vfcen);
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
    read_features<3>(featfilename, &features);

    Tangram::vector<vfcen_t<3>> vfcen(ncells);
    InFeatureEvaluator<3> feature_functor(features);
    VolfracEvaluator<3, Tangram::Jali_Mesh_Wrapper>
        volfrac_calculator(mesh_wrapper, feature_functor, ptol);

    gettimeofday(&begin_core, 0);

    // Populate the vf array using the volfrac functor
    Tangram::vector<int> counter(ncells);
    for (int i = 0; i < ncells; i++) counter[i] = i;
    Tangram::transform(counter.begin(), counter.end(), vfcen.begin(),
                       volfrac_calculator);

    gettimeofday(&end_core, 0);

    writeAsciiFile<3>(outfilename, vfcen);
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
