/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#ifndef UTILITY_COMPUTE_R2D_MOMENTS_H
#define UTILITY_COMPUTE_R2D_MOMENTS_H

#include <vector>

#include "wonton/support/CoordinateSystem.h"

extern "C" {
#include "r2d.h"
}

/* !
  @brief Compute coordinate system aware moments and place them in a C-array
*/

template<class CoordSys>
inline 
std::vector<double> compute_r2d_moments(r2d_poly* poly)
{
  // compute moments
  bool flag = std::is_same<CoordSys, Wonton::CylindricalAxisymmetricCoordinates>::value;
  int poly_order = (flag) ? 2 : 1;
  int nmoments = R2D_NUM_MOMENTS(poly_order);
  r2d_real r2d_moments[nmoments];

  r2d_reduce(poly, r2d_moments, poly_order);

  // shift moments
  std::vector<double> aux(nmoments);
  for (int i = 0; i < nmoments; ++i) aux[i] = r2d_moments[i];

  if (flag) {
    Wonton::CylindricalAxisymmetricCoordinates::shift_moments_list<2>(aux);
  }
  return aux;
}

#endif
