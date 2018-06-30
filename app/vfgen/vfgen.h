// Generate volume fractions for mesh cells given a configuration file
// with different geometric features.
//
//
// Assumption: 2D meshes are in X-Y plane only
//
// Author: Rao Garimella, rao@lanl.gov
// Copyright Los Alamos National Laboratory, 2016

#ifndef TGENVOLFRAC_H_
#define TGENVOLFRAC_H_


#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <string>

#include <omp.h>
#include <mpi.h>

#ifdef THRUST
#include "thrust/device_vector.h"
#include "thrust/transform.h"
#include "thrust/for_each.h"
#include "thrust/reduce.h"
#endif

#include "tangram/support/tangram.h"




#define MAXMATS 100



// Structure to hold volume fraction and centroid data for
// materials in a cell. Assume a maximum of MAXMATS materials per cell

template <int dim>
struct vfcen_t {
  int nmats = 0;
  int matids[MAXMATS] = {};  // initialize to zero
  double vf[MAXMATS] = {};  // initialize to zero
  Tangram::Point<dim> cen[MAXMATS];  // Constructor initializes to origin
};


#define NPARTICLES 1000


// Use Jordan curve algorithm to determine the number of crossings of
// the polygon a ray shot in the +X direction from the test point
// has. An odd number of crossings implies the point is in the polygon
// and an even number of crossings implies its outside
//
// We don't care how boundary points are classified and we don't
// consider tolerance
//
// @param ptest      Test point
// @param np         Number of points in polygon
// @param polypnts   Polygon points
// @return bool      true if point is inside

bool P_InPoly2D(Tangram::Point<2> ptest, int np,
                 Tangram::Point<2> *polypnts) {
  /* Basic test - will work for strictly interior and exterior points */

  double x = ptest[0];
  double y = ptest[1];

  int c = 0;
  for (int i = 0; i < np; i++) {
    int ip1 = (i+1)%np;
    double x_i = polypnts[i][0];
    double y_i = polypnts[i][1];
    double x_iplus1 = polypnts[ip1][0];
    double y_iplus1 = polypnts[ip1][1];

    if (((y_i > y && y_iplus1 <= y) || (y_iplus1 > y && y_i <= y)) &&
        (x <= (x_i + (y-y_i)*(x_iplus1-x_i)/(y_iplus1-y_i))))
      c = !c;
  }

  return (c == 1 ? true : false);
}

// Check if point is inside a polyhedron represented by a triangular
// facetization of its boundary. No tolerance is used for the in/out check
//
// @param ptest     Test point
// @param ntri      Number of triangles
// @param tripnts   Triangle points  (3*ntri)
// @return bool     if point is inside

bool P_InTriPoly3D(Tangram::Point<3> ptest,
                       int ntri, Tangram::Point<3> *tripnts) {
  for (int t = 0; t < ntri; t++) {
    Tangram::Point<3>& p0 = tripnts[3*t];
    Tangram::Point<3>& p1 = tripnts[3*t+1];
    Tangram::Point<3>& p2 = tripnts[3*t+2];
    Tangram::Vector<3> v0 = p0 - ptest;
    Tangram::Vector<3> v1 = p1 - ptest;
    Tangram::Vector<3> v2 = p2 - ptest;
    Tangram::Vector<3> vcross = cross(v0, v1);
    double vol = dot(vcross, v2);
    if (vol < 0.0) return false;
  }
  return true;
}


// Check if point is inside 2D polygon or 3D polyhedron. The 3D
// polyhedron must be specified with triangular facets in the
// boundary. Accordingly, the meaning of the points array is different
// for 2D and 3D. In 2D, it is just a list of points of the
// polygon. In 3D, it is a flat array of points of the triangle
// specified such that the normal of each triangle points out of the
// polyhedron.
//
// @param ptest        Test point
// @param npnts        Number of points sent in
// @param points       Points of polygon in 2D or points of triangles of
//                     of faceted polyhedron in 3D

template<int dim>
bool P_InPoly(Tangram::Point<dim> ptest, int npnts, Tangram::Point<dim> *points) {}

template<>
bool P_InPoly<2>(Tangram::Point<2> ptest, int npnts, Tangram::Point<2> *points) {
  return P_InPoly2D(ptest, npnts, points);
}

template<>
bool P_InPoly<3>(Tangram::Point<3> ptest, int npnts, Tangram::Point<3> *points) {
  int ntris = npnts/3;
  return P_InTriPoly3D(ptest, ntris, points);
}


// Geometric Features

enum class FEATURETYPE {FILL, HALFSPACE, BOX, POLY, SPHERE};


// Class to represent features. Didn't want to use inherited classes
// here because we want this on the GPU

template<int dim>
struct FEATURE {
  constexpr static int MAXPV2 = 100;  // Maximum number of points per face
  constexpr static int MAXPF3 = 100;  // Maximum number of faces per region

  FEATURETYPE type;
  int matid;            /* ID of material */
  int inout;            /* Do we consider the outside or inside */
                        /* For halfplane, points inside are in the back */

  /* points of polygons and polyhedra - for polygons, these are
   * expected to be in ccw manner */
  int nppoly;
  Tangram::Point<dim> polyxyz[MAXPV2];

  /* Additional info for 3D Polyhedron description - Not filled in for 2D */
  int nfpoly;           /* Number of polyhedron faces */
  int nfpnts[MAXPF3];    /* Number of points for each face */
  int fpnts[MAXPF3][MAXPV2]; /* face vertex indices listed so that
                                face normal points out of polyhedron */

  /* Info about computed triangulation of polyhedron - polytrixyz is a
   * flat array containing coordinates of the three vertices of each
   * triangle. It is derived from nfpoly, nfpnts and fpnts - has no
   * meaning and not filled in 2D */
  int ntris;
  Tangram::Point<dim> polytrixyz[MAXPF3*MAXPV2*3];


  /* Circle or Sphere */
  Tangram::Point<dim> cen;
  double radius;

  /* Half-space defined by plane */
  Tangram::Point<dim> plane_xyz;
  Tangram::Point<dim> plane_normal;

  /* Box */
  Tangram::Point<dim> minxyz;     /* Lower left corner of box */
  Tangram::Point<dim> maxxyz;     /* Upper right corner of box */
};



// Functor to evaluate which feature/material a point is in. Features
// are layered on top of another in the order they are listed

// Functor to evaluate which feature/material a point is in. Features
// are layered on top of another in the order they are listed


template <int dim>
struct InFeatureEvaluator {
  int nfeat_;

  std::vector<FEATURE<dim>> features_;

  explicit InFeatureEvaluator(std::vector<FEATURE<dim>> const& features) :
      features_(features), nfeat_(features.size()) {}

  int operator() (Tangram::Point<dim> const& ptxyz) {
    int pmatid = -1;
    int imat = 1;
    bool test_in = true;
    for (int j = nfeat_-1; j >= 0; j--) {
      FEATURE<dim> curfeat = features_[j];
      imat = curfeat.matid;
      test_in = (curfeat.inout == 1);

      bool ptin = true;
      if (curfeat.type == FEATURETYPE::FILL) {  /* Fill */
        test_in = true;  // Nothing to do really
      } else if (curfeat.type == FEATURETYPE::HALFSPACE) { /* Halfspace */
        // variable to store the dot product
        double dot_prod = 0.0;
        for (int i = 0; i < dim; i++)
          // Calcualte the dot product between vector
          // created by point on plane and normal vector
          dot_prod += (ptxyz[i] - curfeat.plane_xyz[i])*curfeat.plane_normal[i];

        // If dot product is negative, point is inside the halspace
        ptin = (dot_prod <= 0.0);
      } else if (curfeat.type == FEATURETYPE::BOX) { /* Box */
        for (int i = 0; i < dim; i++)
          ptin &= (curfeat.minxyz[i] < ptxyz[i] &&
                   curfeat.maxxyz[i] > ptxyz[i]);
      } else if (curfeat.type == FEATURETYPE::POLY) { /* Polyhedron */
        // Need to explicitly specify <dim> because if dim is interpreted as
        // a long int in the calling routine the template deduction fails
        // See https://stackoverflow.com/questions/32891471/template-argument-deduction-which-compiler-is-right-here
        if (dim == 2) {
          ptin = P_InPoly<dim>(ptxyz, curfeat.nppoly, &(curfeat.polyxyz[0]));
        } else {
          ptin = P_InPoly<dim>(ptxyz, curfeat.nppoly, &(curfeat.polytrixyz[0]));
        }
      } else if (curfeat.type == FEATURETYPE::SPHERE) { /* Circle */
        Tangram::Vector<dim> v = ptxyz-curfeat.cen;
        ptin = (v.norm(false) < curfeat.radius*curfeat.radius);  // norm does not compute square root by default
      }
      else {
        std::cerr << "Unknown feature type\n";
        continue;
      }
      // If point is inside material, save the id
      if (test_in == ptin) {
        pmatid = imat;
        break;
      }
    }
    return pmatid;
  }  // operator()


};  // struct InFeatureEvaluator




// General Volume Fraction Evaluator Class for cells
// Computes the volume fractions and centroids for all materials in a cell

template <int dim, class Mesh_Wrapper>
class VolfracEvaluator {
 public:
  // Constructor
  VolfracEvaluator(Mesh_Wrapper const& mesh,
                   InFeatureEvaluator<dim> feature_evaluator,
                   double ptol, int nmats) :
    mesh_(mesh), feature_evaluator_(feature_evaluator), ptol_(ptol), nmats_(nmats)
  {}

  vfcen_t<dim> operator()(int entity_ID) {
    std::cerr << "Operator not evaluated for dim " << dim << "\n";
    vfcen_t<dim> vfcen;
    return vfcen;
  }
 private:
  Mesh_Wrapper const& mesh_;
  InFeatureEvaluator<dim> feature_evaluator_;
  double ptol_;
  double nmats_;
};

// Specialization for dim = 2
template <class Mesh_Wrapper>
class VolfracEvaluator<2, Mesh_Wrapper> {
 public:
  // Constructor
  VolfracEvaluator(Mesh_Wrapper const& mesh,
                   InFeatureEvaluator<2> feature_evaluator,
                   double ptol, int nmats) :
    mesh_(mesh), feature_evaluator_(feature_evaluator), ptol_(ptol), nmats_(nmats)
  {}

  // Operator to calculate volume fractions and centroids of materials
  // in entity_ID

  vfcen_t<2> operator()(int cellID) {
    unsigned int seed = 0;
    double XMIN, XMAX, YMIN, YMAX;

    std::vector<Tangram::Point<2>> fxyz;
    mesh_.cell_get_coordinates(cellID, &fxyz);
    int nfv = fxyz.size();

    XMIN = YMIN = 1.0e+20;
    XMAX = YMAX = -XMIN;
    for (int j = 0; j < nfv; j++) {
      XMIN = std::min(XMIN, fxyz[j][0]);
      XMAX = std::max(XMAX, fxyz[j][0]);
      YMIN = std::min(YMIN, fxyz[j][1]);
      YMAX = std::max(YMAX, fxyz[j][1]);
    }
    double XLEN = XMAX-XMIN;
    double YLEN = YMAX-YMIN;

    int pmatid[NPARTICLES];
    for (int i = 0; i < NPARTICLES; i++)
      pmatid[i] = -1;

    srand(cellID);

    Tangram::Point<2> ptxyz[NPARTICLES];
    int np = 0;
    double xmult = XLEN/RAND_MAX;
    double ymult = YLEN/RAND_MAX;
    while (np < NPARTICLES) {
      ptxyz[np][0] = XMIN + rand_r(&seed)*xmult;
      ptxyz[np][1] = YMIN + rand_r(&seed)*ymult;

      // check if point is in cell. if its outside, don't increment np
      // - pnt will be overwritten

      if (P_InPoly2D(ptxyz[np], nfv, &(fxyz[0])))
        np++;
    }

    // Compute a material ID for each point based on their inclusion
    // in a feature

    Tangram::transform(ptxyz, ptxyz + np, pmatid, feature_evaluator_);

    // Tally up the particles to compute volume fractions

    vfcen_t<2> vfcen = {};  // Initialize to zero
    int npmat[MAXMATS] = {};
    for (int p = 0; p < np; p++) {
      bool found = false;
      int im = 0;
      for (im = 0; im < vfcen.nmats; im++)
        if (vfcen.matids[im] == pmatid[p]) {
          found = true;
          break;
        }
      // check if number of materials haven't exceed the global count
      if ((!found) && (vfcen.nmats < nmats_)) {
        im = vfcen.nmats++;
        vfcen.matids[im] = pmatid[p];
      }
      npmat[im]++;
      vfcen.cen[im] += ptxyz[p];
    }
    for (int im = 0; im < vfcen.nmats; im++) {
      vfcen.vf[im] = ((double) npmat[im])/np;
      vfcen.cen[im] /= npmat[im];
    }

    return vfcen;
  }

 private:
  Mesh_Wrapper const& mesh_;
  InFeatureEvaluator<2> feature_evaluator_;
  double ptol_;
  int nmats_;
};

template <class Mesh_Wrapper>
class VolfracEvaluator<3, Mesh_Wrapper> {
 public:
  // Constructor
  VolfracEvaluator(Mesh_Wrapper const& mesh,
                   InFeatureEvaluator<3> feature_evaluator,
                   double ptol, int nmats) :
    mesh_(mesh), feature_evaluator_(feature_evaluator), ptol_(ptol), nmats_(nmats)
  {}

  // Operator to compute volume fractions and centroids of materials
  // for this entity ID

  vfcen_t<3> operator()(int cellID) {
    unsigned int seed = 0;

    // Compute min-max bounds for cell

    double XMIN, XMAX, YMIN, YMAX, ZMIN, ZMAX;

    XMIN = YMIN = ZMIN = 1.0e+20;
    XMAX = YMAX = ZMAX = -XMIN;

    std::vector<Tangram::Point<3>> rxyz;
    mesh_.cell_get_coordinates(cellID, &rxyz);
    int nrv = rxyz.size();

    for (int i = 0; i < nrv; i++) {
      XMIN = std::min(XMIN, rxyz[i][0]);
      XMAX = std::max(XMAX, rxyz[i][0]);
      YMIN = std::min(YMIN, rxyz[i][1]);
      YMAX = std::max(YMAX, rxyz[i][1]);
      ZMIN = std::min(ZMIN, rxyz[i][2]);
      ZMAX = std::max(ZMAX, rxyz[i][2]);
    }

    double XLEN = XMAX-XMIN;
    double YLEN = YMAX-YMIN;
    double ZLEN = ZMAX-ZMIN;

    // Get a triangular facetization of cell boundary
    std::vector<std::vector<int>> tripnts;
    std::vector<Tangram::Point<3>> points;

    // ------ Will be replaced by call to mesh_.cell_get_facetization -------
    // ------ when PR-26 will get merged in ---------------------------------

    // mesh_.cell_get_facetization(cellID, &tripnts, &points);

    cell_get_facetization(cellID, &tripnts, &points);


    int ntris = tripnts.size();
    std::vector<Tangram::Point<3>> tripnts_flat(3*ntris);
    for (int t = 0; t < ntris; t++) {
      tripnts_flat[3*t] = points[tripnts[t][0]];
      tripnts_flat[3*t+1] = points[tripnts[t][1]];
      tripnts_flat[3*t+2] = points[tripnts[t][2]];
    }

    int pmatid[NPARTICLES];
    for (int i = 0; i < NPARTICLES; i++)
      pmatid[i] = -1;

    srand(cellID);

    /* Throw particles into cell and see which feature they lie in */

    Tangram::Point<3> ptxyz[NPARTICLES];
    int np = 0;
    double xmult = XLEN/RAND_MAX;
    double ymult = YLEN/RAND_MAX;
    double zmult = ZLEN/RAND_MAX;

    while (np < NPARTICLES) {
      ptxyz[np][0] = XMIN + rand_r(&seed)*xmult;
      ptxyz[np][1] = YMIN + rand_r(&seed)*ymult;
      ptxyz[np][2] = ZMIN + rand_r(&seed)*zmult;

      // Check if point is in cell - point can be outside of cell if
      // cell is not a coordinate aligned box if its outside, don't
      // increment np - pnt will be overwritten

      if (P_InTriPoly3D(ptxyz[np], ntris, &(tripnts_flat[0])))
        np++;
    }

    // Compute a material ID for each point based on their inclusion
    // in a feature

    Tangram::transform(ptxyz, ptxyz + NPARTICLES, pmatid, feature_evaluator_);

    // Tally up the particles to compute volume fractions

    vfcen_t<3> vfcen = {};  // Initialize to zero
    int npmat[MAXMATS] = {};
    for (int p = 0; p < np; p++) {
      bool found = false;
      int im = 0;
      for (im = 0; im < vfcen.nmats; im++)
        if (vfcen.matids[im] == pmatid[p]) {
          found = true;
          break;
        }
      // check if number of materials haven't exceed the global count
      if ((!found) && (vfcen.nmats < nmats_)) {
        im = vfcen.nmats++;
        vfcen.matids[im] = pmatid[p];
      }
      npmat[im]++;
      vfcen.cen[im] += ptxyz[p];
    }
    for (int im = 0; im < vfcen.nmats; im++) {
      vfcen.vf[im] = ((double) npmat[im])/np;
      vfcen.cen[im] /= npmat[im];
    }

    return vfcen;
  }

 private:
  Mesh_Wrapper const& mesh_;
  InFeatureEvaluator<3> feature_evaluator_;
  double ptol_;
  int nmats_;


  // Temporary copy of routine from AuxMeshTopology.h. Can be eliminated
  // when we merge in Tangram PR-26

  void cell_get_facetization(int const cellid,
                             std::vector<std::vector<int>> *facetpoints,
                             std::vector<Tangram::Point<3>> *points) const {
    facetpoints->clear();
    points->clear();

    std::vector<int> cnodes;
    mesh_.cell_get_nodes(cellid, &cnodes);
    int ncnodes = cnodes.size();

    points->resize(ncnodes);
    for (int n = 0; n < ncnodes; ++n)
      mesh_.node_get_coordinates(cnodes[n], &((*points)[n]));

    std::vector<int> cfaces, cfdirs;
    mesh_.cell_get_faces_and_dirs(cellid, &cfaces, &cfdirs);
    int ncfaces = cfaces.size();

    for (int f = 0; f < ncfaces; f++) {
      std::vector<int> fnodes;
      mesh_.face_get_nodes(cfaces[f], &fnodes);
      int nfnodes = fnodes.size();

      // Get the local indices (in the cell node list) of the face nodes
      std::vector<int> fnodes_local(nfnodes);
      for (int n = 0; n < nfnodes; n++) {
        bool found = false;
        int i = 0;
        while (!found && i < ncnodes) {
          found = (fnodes[n] == cnodes[i]);
          if (!found) i++;
        }
        assert(found);
        fnodes_local[n] = i;
      }

      if (nfnodes == 3) {  // Triangle; guaranteed to be planar - don't split
        if (cfdirs[f] != 1) {  // reverse direction of nodes
          int tmp = fnodes_local[0];
          fnodes_local[0] = fnodes_local[1];
          fnodes_local[1] = tmp;
        }
        facetpoints->emplace_back(fnodes_local);
      } else {
        // quad or more general polygonal face which could be
        // curved. Facetize/Triangulate it by connecting each edge to a
        // central point in the face. This central point is computed as
        // the geometric center of the nodes of the face

        // Add centroid of face a new point to the point list
        Tangram::Point<3> fcen;
        mesh_.face_centroid(cfaces[f], &fcen);
        points->push_back(fcen);
        int icen = points->size() - 1;

        // Add the triangular facets formed using edges of face and centroid
        std::vector<int> fctpnts(3);
        for (int n = 0; n < nfnodes; n++) {
          if (cfdirs[f] == 1) {
            fctpnts[0] = fnodes_local[n];
            fctpnts[1] = fnodes_local[(n+1)%nfnodes];
          } else {
            fctpnts[0] = fnodes_local[(n+1)%nfnodes];
            fctpnts[1] = fnodes_local[n];
          }
          fctpnts[2] = icen;
          facetpoints->push_back(fctpnts);
        }
      }
    }  // for (f...)
  }  // cell_get_facetization

};


// Write to ascii file - ascii file assumes a maximum of 10 materials
// It also puts entries for all materials for all cells

template<int dim>
void writeAsciiFile(std::string filename,
                    Tangram::vector<vfcen_t<dim>> vfcen,
                    int nmats_) {
  std::ofstream outfile;
  outfile.open(filename.c_str());
  if (!outfile.is_open()) {
    std::cerr << "Could not open output file " << filename << "\n";
    exit(-98);
  }

  int ncells = vfcen.size();

  outfile << dim << "\n";
  outfile << nmats_ << "\n";
  for (int i = 0; i < nmats_; i++)
    outfile << "mat" << i << "\n";

  /* Write out volume fractions and centroids to file */
  int ncells_nmats[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  for (int i = 0; i < ncells; i++) {
    double sum = 0.0;

    vfcen_t<dim> vfcen_i = vfcen[i];  // device to host copy, if Thrust is used

    int nmats_cell = vfcen_i.nmats;
    outfile << nmats_cell << "\n";
    for (int imat = 0; imat < nmats_cell; imat++)
      outfile << vfcen_i.matids[imat] << " ";
    outfile << "\n";
    for (int imat = 0; imat < nmats_cell; imat++)
      outfile << vfcen_i.vf[imat] << " ";
    outfile << "\n";
    for (int imat = 0; imat < nmats_cell; imat++) {
      for (int d = 0; d < dim; d++)
        outfile << vfcen_i.cen[imat][d] << " ";
      outfile << "  ";
    }
    outfile << "\n";
  }
  outfile.close();
}

// Write to binary file - ascii file assumes a maximum of 10 materials
// It also puts entries for all materials for all cells
//
// Format
// <int, dimension>
// <int, number of cells>
// for every cell:
//   <int, number of material in the cell>
//   for every cell's material:
//     <int, index of material>
// for every cell:
//   if cell is multi-material:
//     for every cell's material:
//       <double, material's volume fraction>
// for every cell:
//   if cell is multi-material:
//   for every cell's material:
//     <double, material centroids x-coordinate>
//     if 2D or 3D:
//       <double, material centroids y-coordinate>
//     if 3D:
//       <double, material centroids z-coordinate>

template<int dim>
void writeBinaryFile(std::string filename, Tangram::vector<vfcen_t<dim>> vfcen) {
  std::ofstream outfile;
  outfile.open(filename.c_str(), std::ofstream::out | std::ofstream::binary);
  if (!outfile.is_open()) {
    std::cerr << "Could not open output file " << filename << "\n";
    exit(-98);
  }

  int ncells = vfcen.size();

  int dimcpy = dim;
  outfile.write((char *) &dimcpy, sizeof(int));
  outfile.write((char *) &ncells, sizeof(int));

  /* Write out volume fractions and centroids to file */
  int ncells_nmats[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  for (int i = 0; i < ncells; i++) {
    vfcen_t<dim> const& vfcen_i = vfcen[i];
    int nmats_cell = vfcen_i.nmats;
    outfile.write((char *) &nmats_cell, sizeof(int));
    outfile.write((char *) &(vfcen_i.matids[0]), vfcen_i.nmats*sizeof(int));
  }
  for (int i = 0; i < ncells; i++) {
    vfcen_t<dim> const& vfcen_i = vfcen[i];
    if (vfcen_i.nmats > 1)
      outfile.write((char *) &(vfcen_i.vf[0]), vfcen_i.nmats*sizeof(double));
  }
  for (int i = 0; i < ncells; i++) {
    vfcen_t<dim> const& vfcen_i = vfcen[i];
    if (vfcen_i.nmats > 1)
      for (int imat = 0; imat < vfcen_i.nmats; imat++)
        outfile.write((char *) &(vfcen_i.cen[imat][0]), dim*sizeof(double));
  }
  outfile.close();
}



// Reading the feature file (feature files ends with keyword 'end')
// Feature file example (adjust coordinates for 2D as necessary)
// Example assumes mesh on a domain that goes from (0,0,0) to (1,1,1)
//-------------------------------------
// nmats N
// fill 0
//
// # SPHERE (in 2D use 'circle' and a 2D center)
// # 'in' means consider inside of sphere
// # Instead of 'in' one can use 'out'
// # The last number is material ID
// #
// sphere in 1
// 0.5 0.5 0.5
// 0.3
// #
// # PLANE
// # 'front' means consider the front of the plane
// # i.e. in direction of normal
// # Alternative to 'front' is of course 'back'
// # to mean behind the plane
// # Plane is specified by point and normal
// #
// plane front 2
// 0.5 0.5 0.5
// 1.0 1.0 2.0
//
// # BOX
// # Box is specified by lower left (bottom) corner
// # and upper right (top) corner
// #
// box in 3
// 0.1 0.1 0.1
// 0.2 0.2 0.3
//
// # POLYHEDRON
// # Specified by number of points
// # point coordinates,
// # number of faces
// # number of points and indices of points of each face
// # (use 0 based indexing)
// # FOR POLYGON, only number of points, coordinates is needed
// #
// poly in 4
// 5
// 0.7 0.7 0.2
// 0.8 0.7 0.2
// 0.8 0.8 0.2
// 0.7 0.8 0.2
// 0.75 0.75 0.5
// 5
// 4  0 3 2 1
// 3  0 1 4
// 3  1 2 4
// 3  2 3 4
// 3  3 0 4
// #
// end

template <int dim>
void read_features(std::string featfilename,
                   std::vector<FEATURE<dim>> *features,
                   int *nmats_) {
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
    featfile >> *nmats_;
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

      // Read points of the polygon or polyhedron
      featfile >> this_feature.nppoly;
      for (int j = 0; j < this_feature.nppoly; j++) {
        double polyx, polyy;
        featfile >> this_feature.polyxyz[j][0];
        featfile >> this_feature.polyxyz[j][1];
      }

      if (dim == 3) {
        // In 3D, we have to read the faces and the face points too
        featfile >> this_feature.nfpoly;
        for (int f = 0; f < this_feature.nfpoly; f++) {
          featfile >> this_feature.nfpnts[f];
          for (int p = 0; p < this_feature.nfpnts[f]; p++)
            featfile >> this_feature.fpnts[f][p];
        }

        // Now make a triangular faceted version of the polyhedron
        // Triangular facets assumed to point out of polyhedron
        this_feature.ntris = 0;
        int ntripnts = 0;
        for (int f = 0; f < this_feature.nfpoly; f++) {
          Tangram::Point<dim> fcen;   // will only come here when dim == 3
          int nfp = this_feature.nfpnts[f];
          for (int p = 0; p < nfp; p++)
            fcen += this_feature.polyxyz[this_feature.fpnts[f][p]];
          fcen /= this_feature.nfpnts[f];

          for (int p = 0; p < nfp; p++) {
            int pid0 = this_feature.fpnts[f][p];
            this_feature.polytrixyz[ntripnts] = this_feature.polyxyz[pid0];
            int pid1 = this_feature.fpnts[f][(p+1)%nfp];
            this_feature.polytrixyz[ntripnts+1] = this_feature.polyxyz[pid1];
            this_feature.polytrixyz[ntripnts+2] = fcen;
          }
          ntripnts += 3;
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



#endif
