/*
This file is part of the Ristra tangram project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/tangram/blob/master/LICENSE
*/



#ifndef TANGRAM_LS_FITS_H_
#define TANGRAM_LS_FITS_H_

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#include "tangram/support/Point.h"
#include "tangram/support/Matrix.h"

namespace Tangram {


/*!
  @brief Compute least squares gradient from set of values
  @param[in] coords Vector of coordinates at which values are given
  @param[in] vals   Vector of values at said coordinates

  Compute a least squares gradient from a set of values. The first
  point is assumed to be the point where the gradient must be computed
  and the first value is assumed to the value at this reference point

  This operator does not know anything about a mesh.

*/

template<long D>
Vector<D> ls_gradient(std::vector<Point<D>> const & coords,
                      std::vector<double> const & vals) {

  Point<D> coord0 = coords[0];

  double val0 = vals[0];

  // There are nvals but the first is the reference point where we
  // are trying to compute the gradient; so the matrix sizes etc
  // will only be nvals-1

  int nvals = vals.size();

  // Each row of A contains the components of the vector from
  // coord0 to the candidate point being used in the Least Squares
  // approximation (X_i-X_0).

  Matrix A(nvals-1, D);
  for (int i = 0; i < nvals-1; ++i) {
    for (int j = 0; j < D; ++j)
      A[i][j] = coords[i+1][j]-coord0[j];
  }


  // A is a matrix of size nvals-1 by D (where D is the space
  // dimension). So transpose(A)*A is D by D

  Matrix AT = A.transpose();

  Matrix ATA = AT*A;

  // Each entry/row of F contains the difference between the
  // function value at the candidate point and the function value
  // at the point where we are computing (f-f_0)

  std::vector<double> F(nvals-1);
  for (int i = 0; i < nvals-1; ++i)
    F[i] = vals[i+1]-val0;

  // F is a vector of nvals. So transpose(A)*F is vector of D
  // (where D is the space dimension)

  Vector<D> ATF = Vector<D>(AT*F);

  // Inverse of ATA

  Matrix ATAinv = ATA.inverse();

  // Gradient of length D

  return ATAinv*ATF;
}

}  // namespace Tangram

#endif
