/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#ifndef TANGRAM_SUPPORT_RK4_H_
#define TANGRAM_SUPPORT_RK4_H_

#include <functional>

// tangram includes
#include "tangram/support/tangram.h"

// wonton includes
#include "wonton/support/wonton.h"

namespace Tangram {

/*!
  @brief Fourth order Runge–Kutta method for initial value problems
  @param[in] roc_f  Rate of change function for x, dx/dt=roc_f(t, x)
  @param[in] x_cur  Current value of x
  @param[in] t_cur Current time
  @param[in] dt Time step
  @return  Value of x after the time step
*/
template<int D>
Point<D> runge_kutta_4(const std::function<Vector<D>(double, const Point<D>&)>& roc_f,
                       const Point<D>& x_cur,
                       double t_cur,
                       double dt) {
  Vector<D> k1 = dt*roc_f(t_cur,          x_cur         );
  Vector<D> k2 = dt*roc_f(t_cur + 0.5*dt, x_cur + 0.5*k1);
  Vector<D> k3 = dt*roc_f(t_cur + 0.5*dt, x_cur + 0.5*k2);
  Vector<D> k4 = dt*roc_f(t_cur +     dt, x_cur +     k3);

  return x_cur + (k1 + 2*k2 + 2*k3 + k4)/6.0;
}

}  // namespace Tangram

#endif  // TANGRAM_SUPPORT_RK4_H_
