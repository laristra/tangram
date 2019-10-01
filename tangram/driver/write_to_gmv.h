/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#ifndef TANGRAM_DRIVER_WRITE_TO_GMV_H_
#define TANGRAM_DRIVER_WRITE_TO_GMV_H_

#include <fstream>
#include <vector>
#include <array>
#include <algorithm>
#include <string>

#include "tangram/support/tangram.h"
#include "tangram/driver/CellMatPoly.h"

namespace Tangram {

template<class Mesh_Wrapper, int D>
void write_to_gmv(Mesh_Wrapper const& mesh,
                  int nmats,
                  std::vector<int> const& cell_num_mats,
                  std::vector<int> const& cell_mat_ids,
                  std::vector<std::shared_ptr<CellMatPoly<D>>> cellmatpoly_list,
                  std::string filename) {

  int nc = mesh.num_entities(Entity_kind::CELL, Entity_type::PARALLEL_OWNED);
  std::vector<int> cell_offsets(nc);
  cell_offsets[0] = 0;
  for (int c = 1; c < nc; c++)
    cell_offsets[c] = cell_offsets[c-1] + cell_num_mats[c-1];

  std::ofstream fout(filename);
  fout << std::scientific;
  fout.precision(17);

  fout << "gmvinput ascii" << std::endl;
  fout << "codename Tangram" << std::endl;
  fout << "simdate 01/01/01" << std::endl;

  // Make a list of unique points in the mesh and cell mat poly
  // structures combined

  int np = mesh.num_entities(Entity_kind::NODE, Entity_type::PARALLEL_OWNED);
  int nmatpnts = np;
  for (int c = 0; c < nc; c++) {
    CellMatPoly<D> *cellmatpoly = cellmatpoly_list[c].get();

    if (cellmatpoly) {  // Mixed cell
      int ncp = cellmatpoly->num_matvertices();
      for (int i = 0; i < ncp; i++)
        if (cellmatpoly->matvertex_parent_kind(i) != Entity_kind::NODE)
          nmatpnts++;
    }
  }

  std::vector<Point<D>> points(nmatpnts);

  for (int i = 0; i < np; i++)
    mesh.node_get_coordinates(i, &(points[i]));

  nmatpnts = np;  // reset nmatpnts so it can be used as a counter
  for (int c = 0; c < nc; c++) {
    CellMatPoly<D> *cellmatpoly = cellmatpoly_list[c].get();

    if (cellmatpoly) {  // Mixed cell
      int ncp = cellmatpoly->num_matvertices();
      for (int i = 0; i < ncp; i++) {
        if (cellmatpoly->matvertex_parent_kind(i) != Entity_kind::NODE) {
          points[nmatpnts] = cellmatpoly->matvertex_point(i);
          nmatpnts++;
        }
      }
    }
  }

  fout << "nodev " << nmatpnts << std::endl;
  for (int ip = 0; ip < nmatpnts; ip++) {
    fout << points[ip];
    if (D == 2)  // GMV requires us to output 3 coordinates
      fout << " " << 0.000000;
    fout << std::endl;
  }


  // Count the number of polygons we will write out

  int npoly = 0;
  for (int c = 0; c < nc; c++) {
    CellMatPoly<D> *cellmatpoly = cellmatpoly_list[c].get();
    if (cellmatpoly)
      npoly += cellmatpoly->num_matpolys();
    else
      npoly += 1;
  }

  fout << "cells " << npoly << std::endl;

  for (int c = 0; c < nc; c++) {
    CellMatPoly<D> *cellmatpoly = cellmatpoly_list[c].get();

    if (cellmatpoly) {

      // Write out the material polyhedra
      int nmp = cellmatpoly->num_matpolys();
      for (int i = 0; i < nmp; i++) {

        if (D == 1 || D == 2) {
          std::vector<int> mverts = cellmatpoly->matpoly_vertices(i);
          fout << "general 1 " << mverts.size() << " ";
          for (auto n : mverts) {
            if (cellmatpoly->matvertex_parent_kind(n) == Entity_kind::NODE)
              fout << cellmatpoly->matvertex_parent_id(n)+1 << " ";
            else {
              Point<D> const& pnt = cellmatpoly->matvertex_point(n);
              int np = mesh.num_entities(Entity_kind::NODE);
              for (int j = np; j < nmatpnts; j++) {
                if (points[j] == pnt) {
                  fout << j+1 << " ";
                  break;
                }
              }
            }
          }
          fout << std::endl;
        } else if (D == 3) {
          std::vector<int> const& mfaces = cellmatpoly->matpoly_faces(i);
          fout << "general " << mfaces.size() << std::endl;
          for (auto f : mfaces) {
            std::vector<int> const& mfverts = cellmatpoly->matface_vertices(f);
            fout << mfverts.size() << " ";
          }
          fout << std::endl;
          for (auto f : mfaces) {
            std::vector<int> const& mfverts = cellmatpoly->matface_vertices(f);
            int nfv = mfverts.size();
            int mp0, mp1;
            cellmatpoly->matface_matpolys(f, &mp0, &mp1);
            if (mp0 == i) {  // Natural order of vertices
              for (int j = 0; j < nfv; j++) {
                int n = mfverts[j];
                if (cellmatpoly->matvertex_parent_kind(n) == Entity_kind::NODE)
                  fout << cellmatpoly->matvertex_parent_id(n)+1 << " ";
                else {
                  Point<D> const& pnt = cellmatpoly->matvertex_point(n);
                  int np = mesh.num_entities(Entity_kind::NODE);
                  for (int j = np; j < nmatpnts; j++) {
                    if (points[j] == pnt) {
                      fout << j+1 << " ";
                      break;
                    }
                  }
                }
              }
              fout << std::endl;
            } else {  // Reverse order of vertices
              for (int j = 0; j < nfv; j++) {
                int n = mfverts[nfv-j-1];
                if (cellmatpoly->matvertex_parent_kind(n) == Entity_kind::NODE)
                  fout << cellmatpoly->matvertex_parent_id(n)+1 << " ";
                else {
                  Point<D> const& pnt = cellmatpoly->matvertex_point(n);
                  int np = mesh.num_entities(Entity_kind::NODE);
                  for (int j = np; j < nmatpnts; j++) {
                    if (points[j] == pnt) {
                      fout << j+1 << " ";
                      break;
                    }
                  }
                }
              }
              fout << std::endl;
            }
          }  // for (auto f : mfaces)
        }  // else if (D == 3)
      }  // for (i = 0; i < nmp; i++)
    } else {
      if (D == 1 || D == 2) {
        // Write out the cell
        std::vector<int> cverts;
        mesh.cell_get_nodes(c, &cverts);

        // write cell out as polygons (even if its a quad or tri)
        fout << "general 1 " << cverts.size() << " ";
        for (auto n : cverts)
          fout << n+1 << " ";
        fout << std::endl;
      } else if (D == 3) {
        std::vector<int> cfaces;
        std::vector<int> cfdirs;
        mesh.cell_get_faces_and_dirs(c, &cfaces, &cfdirs);
        fout << "general " << cfaces.size() << std::endl;
        for (auto f : cfaces) {
          std::vector<int> fverts;
          mesh.face_get_nodes(f, &fverts);
          fout << fverts.size() << " ";
        }
        fout << std::endl;

        int j = 0;
        for (auto f : cfaces) {
          std::vector<int> fverts;
          mesh.face_get_nodes(f, &fverts);
          int nfverts = fverts.size();
          if (cfdirs[j] == 1) {
            for (int k = 0; k < nfverts; k++)
              fout << fverts[k]+1 << " ";
            fout << std::endl;
          } else {
            for (int k = 0; k < nfverts; k++)
              fout << fverts[nfverts-k-1]+1 << " ";
            fout << std::endl;
          }
        }
      }
    }  // if (cellmatpoly) ... else
  }  // for (int c = 0; c < nc; c++)

  // Write out material ID of polygons

  fout << "material " << std::endl;
  fout << nmats << " 0" << std::endl;
  for (int i = 0; i < nmats; i++)
    fout << "mat" << i+1 << std::endl;

  for (int c = 0; c < nc; c++) {
    CellMatPoly<D> *cellmatpoly = cellmatpoly_list[c].get();
    if (cellmatpoly) {
      for (int i = 0; i < cellmatpoly->num_matpolys(); i++)
        fout << cellmatpoly->matpoly_matid(i)+1 << " ";
    } else {
      int ofst = cell_offsets[c];
      int matid = cell_mat_ids[ofst];
      fout << matid+1 << " ";
    }
  }
  fout << std::endl;

  fout << "endgmv" << std::endl;
}

//Simplified version: outputs CellMatPoly's list. If CellMatPoly's list
//was obtained from a reconstructor, it corresponds to writing out only
//multi-material cells. Does not require the base mesh and is therefore
//convenient for debugging and/or when a single CellMatPoly needs to be
//written out.
//Note that unlike in the full version, nodes and faces of CellMatPoly's
//are written as is and will be duplicated if shared by different CellMatPoly's.
template<int D>
void write_to_gmv(const std::vector<std::shared_ptr<CellMatPoly<D>>>& cellmatpoly_list,
                  const std::string& filename) {
  std::ofstream fout(filename);
  fout << std::scientific;
  fout.precision(17);

  fout << "gmvinput ascii" << std::endl;
  fout << "codename Tangram" << std::endl;
  fout << "simdate 01/01/01" << std::endl;

  int npoly = 0;
  int nmats = 0;
  std::vector<Point<D>> points;
  int const nb_cell_matpolys = cellmatpoly_list.size();

  for (int imp = 0; imp < nb_cell_matpolys; imp++) {
    CellMatPoly<D> *cellmatpoly = cellmatpoly_list[imp].get();

    if (cellmatpoly) {  // Mixed cell
      int ncp = cellmatpoly->num_matvertices();
      for (int i = 0; i < ncp; i++)
        points.push_back(cellmatpoly->matvertex_point(i));
      npoly += cellmatpoly->num_matpolys();
      const std::vector<int>& matpoly_matids = cellmatpoly->matpoly_matids();
      int max_matid = *std::max_element(matpoly_matids.begin(), matpoly_matids.end());
      if (max_matid >= nmats) nmats = max_matid + 1;
    }
  }

  int nmatpnts = static_cast<int>(points.size());
  fout << "nodev " << nmatpnts << std::endl;
  for (int ip = 0; ip < nmatpnts; ip++) {
    fout << points[ip];
    if (D == 2)  // GMV requires us to output 3 coordinates
      fout << " " << 0.000000;
    fout << std::endl;
  }

  fout << "cells " << npoly << std::endl;
  int pts_offset = 0;

  for (int imp = 0; imp < nb_cell_matpolys; imp++) {
    CellMatPoly<D> *cellmatpoly = cellmatpoly_list[imp].get();
    if (!cellmatpoly)
      continue;

    // Write out the material polyhedra
    int nmp = cellmatpoly->num_matpolys();
    for (int i = 0; i < nmp; i++) {
      if (D == 1 || D == 2) {
        std::vector<int> mverts = cellmatpoly->matpoly_vertices(i);
        fout << "general 1 " << mverts.size() << " ";
        for (auto n : mverts)
          fout << pts_offset + n + 1 << " ";
        fout << std::endl;
      }
      else if (D == 3) {
        std::vector<int> const& mfaces = cellmatpoly->matpoly_faces(i);
        fout << "general " << mfaces.size() << std::endl;
        for (auto f : mfaces) {
          std::vector<int> const& mfverts = cellmatpoly->matface_vertices(f);
          fout << mfverts.size() << " ";
        }
        fout << std::endl;
        for (auto f : mfaces) {
          std::vector<int> const& mfverts = cellmatpoly->matface_vertices(f);
          for (auto n : mfverts)
            fout << pts_offset + n + 1 << " ";
          fout << std::endl;
        }
      }  // else if (D == 3)
    }  // for (i = 0; i < nmp; i++)
    pts_offset += cellmatpoly->num_matvertices();
  }  // for (int imp = 0; imp < cellmatpoly_list.size(); imp++)

  // Write out material ID of polygons

  fout << "material " << std::endl;
  fout << nmats << " 0" << std::endl;
  for (int i = 0; i < nmats; i++)
    fout << "mat" << i + 1 << std::endl;

  for (int imp = 0; imp < nb_cell_matpolys; imp++) {
    CellMatPoly<D> *cellmatpoly = cellmatpoly_list[imp].get();
    if (cellmatpoly) {
      for (int i = 0; i < cellmatpoly->num_matpolys(); i++)
        fout << cellmatpoly->matpoly_matid(i) + 1 << " ";
    }
  }
  fout << std::endl;

  fout << "endgmv" << std::endl;
}

}  // namespace Tangram

#endif  // TANGRAM_DRIVER_WRITE_TO_GMV_H_
