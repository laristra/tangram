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

#ifndef TANGRAM_WRITE_TO_GMV_H_
#define TANGRAM_WRITE_TO_GMV_H_

#include <fstream>
#include <vector>
#include <array>
#include <algorithm>
#include <string>

#include "tangram/support/tangram.h"
#include "tangram/support/Point.h"

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
                if (approxEq(points[j], pnt, 1.0e-16)) {
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
            int dir, mp0, mp1;
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
                    if (approxEq(points[j], pnt, 1.0e-16)) {
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
                    if (approxEq(points[j], pnt, 1.0e-16)) {
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
        int nf;
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


}  // namespace Tangram

#endif  // TANGRAM_WRITE_TO_GMV_H_
