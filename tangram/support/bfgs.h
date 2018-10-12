/*
 This file is part of the Ristra tangram project.
 Please see the license file at the root of this repository, or at:
 https://github.com/laristra/tangram/blob/master/LICENSE
*/

#ifndef TANGRAM_SUPPORT_BFGS_H_
#define TANGRAM_SUPPORT_BFGS_H_

#include <functional>

// tangram includes
#include "tangram/support/tangram.h"

// wonton includes
#include "wonton/support/wonton.h"

namespace Tangram {

  using Wonton::pow2;

/*!
  @brief Computes the minimizer of a quadratic interpolant in the interval
  @param[in] arg0  End of the domain
  @param[in] arg1  Other end of the domain
  @param[in] fval0  Value of the function at arg0
  @param[in] fval1  Value of the function at arg1
  @param[in] df0  Value of the derivative at arg0
  @param[in] safeguards  Two-element array to keep the minimizer away
  from the boundaries of the domain
  @return Value of the minimizer
*/
double quad_interpolant_minimizer(const double arg0,
                                  const double arg1,
                                  const double fval0,
                                  const double fval1,
                                  const double df0,
                                  const double* safeguards) {
  double darg = arg1 - arg0;
  double denominator = fval1 - fval0 - df0*(darg);
  double arg_min;
  double bnd0 = arg0 + safeguards[0]*darg,
         bnd1 = arg0 + safeguards[1]*darg;

  if (std::fabs(denominator) <= std::numeric_limits<double>::epsilon())
    arg_min = (fval1 > fval0) ? bnd0 : bnd1;
  else {
    arg_min = arg0 - df0*pow2(darg)/(2*denominator);
    double d2arg0 = std::fabs(arg_min - arg0);
    if (d2arg0 < std::fabs(bnd0 - arg0)) arg_min = bnd0;
    else if (d2arg0 > std::fabs(bnd1 - arg0)) arg_min = bnd1;
  }

  return arg_min;
}

/*!
  @brief Forward finite difference approximation of a gradient
  @param[in] fun  Function for which we compute the gradient
  @param[in] arg_val  We approximate the gradient
  for this value of the argument
  @param[in] fun_val  Value of function at arg_val
  @param[out] approx_grad  Computed approximation of the gradient at arg_val
  @param[out] fdiff_h  Finite difference step
*/
template<int arg_dim>
void fwd_diff_grad(const std::function<double(const Vector<arg_dim>&)>& fun,
                   const Vector<arg_dim>& arg_val,
                   const double& fun_val,
                   Vector<arg_dim>& approx_grad,
                   double& fdiff_h) {
  fdiff_h = sqrt(std::numeric_limits<double>::epsilon());

  Vector<arg_dim> fwd_arg = arg_val;
  for (int i = 0; i < arg_dim; i++) {
    fwd_arg[i] += fdiff_h;
    double forward_fval = fun(fwd_arg);
    fwd_arg[i] = arg_val[i];

    approx_grad[i] = (forward_fval - fun_val)/fdiff_h;
  }
}

/*!
  @brief Linesearch algorithm: finds the argument satisfying
  the strong Wolfe conditions
  @param[in] obj_fun  Objective function
  @param[in] arg  Current argument
  @param[in] dir  Linesearch direction
  @param[in] initial_guess  Starting value of the argument
  @param[in] c1  Value of c1 parameter in Wolfe conditions
  @param[in] c2  Value of c2 parameter in Wolfe conditions
  @param[in] im_tols  Tolerances
  @param[in/out] cur_fval  Value of the function at arg [in]
  and arg + alpha*dir [out]
  @param[in/out] cur_grad  Value of the gradient at arg [in]
  and arg + alpha*dir [out]
  @return Value of the minimizer
*/
template<int arg_dim>
double linesearch(const std::function<double(const Vector<arg_dim>&)>& obj_fun,
                  const double& obj_fun_lbnd,
                  const Vector<arg_dim>& arg,
                  const Vector<arg_dim>& dir,
                  const double& initial_guess,
                  const double& c1,
                  const double& c2,
                  const IterativeMethodTolerances_t& im_tols,
                  double& cur_fval,
                  Vector<arg_dim>& cur_grad) {
  // Bracketing step, solution will be in (alpha_bnd[0], alpha_bnd[1])
  double fval0 = cur_fval;
  double df0 = dot(cur_grad, dir);
  assert(df0 < -std::numeric_limits<double>::epsilon());
  double fdiff_h;

  auto fval = [&obj_fun, &arg, &dir](double alpha) {
    return obj_fun(arg + alpha*dir);
  };

  // alpha_max is chosen so that it should violate the sufficient decrease condition
  // and therefore provides an acceptable upper bound for the interval
  double alpha_max = (obj_fun_lbnd - fval0)/(c1*df0);

  double alpha_bnd[2] = {0.0, initial_guess};
  double bnd_fval[2]; bnd_fval[0] = fval0;
  double df_lbnd = df0;
  // At the beginning of each loop, cur_fval and cur_grad are at alpha_bnd[0]
  for (int i = 0; i < im_tols.max_num_iter; i++) {
    bnd_fval[1] = fval(alpha_bnd[1]);

    if (bnd_fval[1] - obj_fun_lbnd <= im_tols.fun_eps) {
      cur_fval = bnd_fval[1];
      fwd_diff_grad<arg_dim>(obj_fun, arg + alpha_bnd[1]*dir, bnd_fval[1], cur_grad, fdiff_h);
      return alpha_bnd[1];
    }

    if ( (bnd_fval[1] > fval0 + c1*alpha_bnd[1]*df0) ||
         ((bnd_fval[1] >= bnd_fval[0]) && (i > 0)) )
      break;

    // alpha_bnd[1] will be the next alpha_bnd[0]
    cur_fval = bnd_fval[1];
    fwd_diff_grad<arg_dim>(obj_fun, arg + alpha_bnd[1]*dir, bnd_fval[1], cur_grad, fdiff_h);
    df_lbnd = dot(cur_grad, dir);

    if (std::fabs(df_lbnd) <= -c2*df0)
      return alpha_bnd[1];

    if (df_lbnd >= 0.0) {
        std::reverse(std::begin(alpha_bnd), std::end(alpha_bnd));
        std::reverse(std::begin(bnd_fval), std::end(bnd_fval));
        break;
    }

    double dalpha = alpha_bnd[1] - alpha_bnd[0];
    double tau = 9.0; // interval length scaling

    alpha_bnd[0] = alpha_bnd[1];
    bnd_fval[0] = bnd_fval[1];
    alpha_bnd[1] += tau*dalpha;

    if (alpha_bnd[1] > alpha_max - im_tols.arg_eps) {
      alpha_bnd[1] = alpha_max;
      bnd_fval[1] = obj_fun(alpha_max);
      break;
    }
  }

  // Sectioning/zooming step
  // At the beginning of each loop, cur_fval and cur_grad are at alpha_bnd[0]
  double safeguards[2] = {0.1, 0.5};
  for (int i = 0; i < im_tols.max_num_iter; i++) {
    if (std::fabs(alpha_bnd[1] - alpha_bnd[0]) < im_tols.arg_eps)
      return alpha_bnd[0];

    double cur_alpha = quad_interpolant_minimizer(alpha_bnd[0], alpha_bnd[1], bnd_fval[0],
                                                  bnd_fval[1], df_lbnd, safeguards);

    cur_fval = fval(cur_alpha);
    if (cur_fval - obj_fun_lbnd <= im_tols.fun_eps) {
      fwd_diff_grad<arg_dim>(obj_fun, arg + cur_alpha*dir, cur_fval, cur_grad, fdiff_h);
      return cur_alpha;
    }

    if ( (cur_fval > fval0 + c1*cur_alpha*df0) || (cur_fval >= bnd_fval[0]) ) {
      alpha_bnd[1] = cur_alpha;
      bnd_fval[1] = cur_fval;
      cur_fval = bnd_fval[0];
      continue;
    }

    fwd_diff_grad<arg_dim>(obj_fun, arg + cur_alpha*dir, cur_fval, cur_grad, fdiff_h);
    double cur_df = dot(cur_grad, dir);
    if (std::fabs(cur_df) <= -c2*df0)
      return cur_alpha;

    if (cur_df*(alpha_bnd[1] - alpha_bnd[0]) >= 0) {
      alpha_bnd[1] = alpha_bnd[0];
      bnd_fval[1] = bnd_fval[0];
    }

    alpha_bnd[0] = cur_alpha;
    bnd_fval[0] = cur_fval;
    df_lbnd = cur_df;
  }

  return alpha_bnd[0];
}

/*!
  @brief Linesearch algorithm, based on MWWP algorithm from
  "Global convergence of BFGS and PRP methods under
  a modified weak Wolfe–Powell line search"
  by Gonglin Yuana, Zengxin Weia, Xiwen Lu,
  Applied Mathematical Modelling 47 (2017) 811–825.

  @param[in] obj_fun  Objective function
  @param[in] arg  Current argument
  @param[in] dir  Linesearch direction
  @param[in] initial_guess  Starting value of the argument
  @param[in] delta  Value of delta parameter in MWWP conditions
  @param[in] delta1  Value of delta_1 parameter in MWWP conditions
  @param[in] im_tols  Tolerances
  @param[in/out] cur_fval  Value of the function at arg [in]
  and arg + alpha*dir [out]
  @param[in/out] cur_grad  Value of the gradient at arg [in]
  and arg + alpha*dir [out]

  @return Value of the minimizer
*/
template<int arg_dim>
double mwwp_linesearch(const std::function<double(const Vector<arg_dim>&)>& obj_fun,
                       const double& obj_fun_lbnd,
                       const Vector<arg_dim>& arg,
                       const Vector<arg_dim>& dir,
                       const double initial_guess,
                       const double delta,
                       const double delta1,
                       const IterativeMethodTolerances_t& im_tols,
                       double& cur_fval,
                       Vector<arg_dim>& cur_grad) {
  // Bracketing step, solution will be in (alpha_bnd[0], alpha_bnd[1])
  double phi0 = cur_fval;
  double dir_dot_dir = dir.norm(false);

  double dphi0 = dot(cur_grad, dir);
  assert(dphi0 < -std::numeric_limits<double>::epsilon());
  double fdiff_h;

  auto fval = [&obj_fun, &arg, &dir](double alpha) {
    return obj_fun(arg + alpha*dir);
  };

  auto phi = [&delta, &delta1, &phi0, &dphi0, &dir_dot_dir](double alpha, double f_alpha) {
    return f_alpha - phi0 - delta*dphi0*alpha -
      alpha*std::min(-delta1*dphi0, 0.5*delta*dir_dot_dir*alpha);
  };

  auto dphi = [&dir, &delta, &delta1, &dphi0, &dir_dot_dir](
    double alpha, Vector<arg_dim>& grad_alpha) {
    return dot(grad_alpha, dir) - delta*dphi0 -
      std::min(-delta1*dphi0, delta*dir_dot_dir*alpha);
  };

  double alpha_max = (obj_fun_lbnd - phi0)/((delta - delta1)*dphi0);

  double alpha_bnd[2] = {0.0, initial_guess};
  double bnd_phi[2]; bnd_phi[0] = phi0;
  double dphi_lbnd = dphi0;

  // At the beginning of each loop, cur_fval and cur_grad are at alpha_bnd[0]
  for (int i = 0; i < im_tols.max_num_iter; i++) {
    double ubnd_fval = fval(alpha_bnd[1]);
    if (ubnd_fval - obj_fun_lbnd <= im_tols.fun_eps) {
      cur_fval = ubnd_fval;
      fwd_diff_grad<arg_dim>(obj_fun, arg + alpha_bnd[1]*dir, ubnd_fval, cur_grad, fdiff_h);
      return alpha_bnd[1];
    }

    bnd_phi[1] = phi(alpha_bnd[1], ubnd_fval);
    if (bnd_phi[1] > 0.0)
      break;

    cur_fval = ubnd_fval;
    fwd_diff_grad<arg_dim>(obj_fun, arg + alpha_bnd[1]*dir, ubnd_fval, cur_grad, fdiff_h);
    dphi_lbnd = dphi(alpha_bnd[1], cur_grad);

    if (dphi_lbnd >= 0.0)
      return alpha_bnd[1];

    if (alpha_bnd[1] > alpha_max - im_tols.arg_eps)
      return alpha_bnd[1];


    double dalpha = alpha_bnd[1] - alpha_bnd[0];
    double tau = 9.0; // interval length scaling

    alpha_bnd[0] = alpha_bnd[1];
    bnd_phi[0] = bnd_phi[1];
    alpha_bnd[1] += tau*dalpha;

    if (alpha_bnd[1] > alpha_max - im_tols.arg_eps) {
      alpha_bnd[1] = alpha_max;
      double alpha_max_fval = obj_fun(alpha_max);
      bnd_phi[1] = phi(alpha_max, alpha_max_fval);
      break;
    }
  }

  // Sectioning/zooming step
  // At the beginning of each loop, cur_fval and cur_grad are at alpha_bnd[0]
  double lbnd_fval = cur_fval;
  double safeguards[2] = {0.1, 0.5};
  for (int i = 0; i < im_tols.max_num_iter; i++) {
    if (std::fabs(alpha_bnd[1] - alpha_bnd[0]) < im_tols.arg_eps)
      return alpha_bnd[0];

    double cur_alpha = quad_interpolant_minimizer(alpha_bnd[0], alpha_bnd[1], bnd_phi[0],
                                                  bnd_phi[1], dphi_lbnd, safeguards);

    cur_fval = fval(cur_alpha);
    if (cur_fval - obj_fun_lbnd <= im_tols.fun_eps) {
      fwd_diff_grad<arg_dim>(obj_fun, arg + cur_alpha*dir, cur_fval, cur_grad, fdiff_h);
      return cur_alpha;
    }

    double cur_phi = phi(cur_alpha, cur_fval);
    if ( (cur_phi > 0.0) || (cur_phi >= bnd_phi[0]) ) {
      alpha_bnd[1] = cur_alpha;
      bnd_phi[1] = cur_phi;
      cur_fval = lbnd_fval;
      continue;
    }

    fwd_diff_grad<arg_dim>(obj_fun, arg + cur_alpha*dir, cur_fval, cur_grad, fdiff_h);
    double cur_dphi = dphi(cur_alpha, cur_grad);
    if (cur_dphi >= 0.0)
      return cur_alpha;

    alpha_bnd[0] = cur_alpha;
    bnd_phi[0] = cur_phi;
    lbnd_fval = cur_fval;
    dphi_lbnd = cur_dphi;
  }

  return alpha_bnd[0];
}

/*!
  @brief BFGS algorithm: finds the minimizer of the objective function
  @param[in] obj_fun  Objective function
  @param[in] init_guess  Initial value of the argument
  @param[in] max_num_iter  Maximum number of iterations
  @param[in] im_tols  Tolerances
  @return Value of the minimizer
*/
template<int arg_dim>
Vector<arg_dim> bfgs(const std::function<double(const Vector<arg_dim>&)>& obj_fun,
                     const double& obj_fun_lbnd,
                     const Vector<arg_dim>& init_guess,
                     const IterativeMethodTolerances_t& im_tols) {
  const double linesearch_c1 = 1.0e-4;
  const double linesearch_c2 = 0.9;
  const double damping_scalar = 0.2;

  double fdiff_h;

  //Initial Hessian approximation: identity matrix
  Matrix B(arg_dim, arg_dim, 0.0);
  for (int i = 0; i < arg_dim; i++)
    B[i][i] = 1.0;

  Vector<arg_dim> cur_arg = init_guess, prev_grad, cur_grad, cur_dir, darg, dgrad;
  double dfval, cur_fval = obj_fun(cur_arg);
  fwd_diff_grad<arg_dim>(obj_fun, cur_arg, cur_fval, cur_grad, fdiff_h);

  for (int i = 0; i < im_tols.max_num_iter; i++) {
    // Forward differences are O(fdiff_h) accurate
    if (cur_grad.max_norm()*fdiff_h/std::max(1.0, std::fabs(cur_fval)) <=
        std::numeric_limits<double>::epsilon())
      break;

    Wonton::solve<arg_dim>(B, -cur_grad, cur_dir);

    double df0 = dot(cur_grad, cur_dir);
    // Function is flat in the linesearch direction
    if (df0 >= -std::numeric_limits<double>::epsilon())
      break;

    double initial_alpha = 1.0;
    if (i > 0) initial_alpha = std::min(1.0,
      -2*std::max(-dfval, 10*std::numeric_limits<double>::epsilon())/df0);

    dfval = -cur_fval;
    prev_grad = cur_grad;

    double alpha = linesearch<arg_dim>(obj_fun, obj_fun_lbnd,
                                       cur_arg, cur_dir, initial_alpha,
                                       linesearch_c1, linesearch_c2,
                                       im_tols, cur_fval, cur_grad);
    dfval += cur_fval;
    darg = alpha*cur_dir;

    if (darg.norm() < im_tols.arg_eps)
      break;

    cur_arg += darg;
    if (std::fabs(cur_fval) <= im_tols.fun_eps)
      return cur_arg;

    dgrad = cur_grad - prev_grad;

    // We update B if the change in gradient is above
    // forward differences error
    if (dgrad.max_norm() > fdiff_h) {
      Vector<arg_dim> B_darg = B*darg;
      if (B_darg.max_norm() <= fdiff_h)
        continue;

      double darg_B_darg = dot(darg, B_darg);
      double darg_dgrad = dot(darg, dgrad);

      // Perform damping if necessary
      if (darg_dgrad < damping_scalar*darg_B_darg) {
        double delta = (1.0 - damping_scalar)*darg_B_darg/(darg_B_darg - darg_dgrad);
        dgrad *= delta;
        dgrad += (1.0 - delta)*B_darg;
        darg_dgrad = dot(darg, dgrad);
      }

      if (darg_dgrad <= std::numeric_limits<double>::epsilon())
        continue;

      // Correction to the initial Hessian guess
      if (i == 0) {
        B *= dgrad.norm(false)/darg_dgrad;
        B_darg = B*darg;
        darg_B_darg = dot(darg, B_darg);
      }

      B += (1.0/darg_dgrad)*(dgrad*dgrad) - (1.0/darg_B_darg)*(B_darg*B_darg);
    }
  }

  return cur_arg;
}

/*!
  @brief BFGS algorithm: finds the minimizer of the objective function.
  Uses MWP linesearch.
  @param[in] obj_fun  Objective function
  @param[in] init_guess  Initial value of the argument
  @param[in] max_num_iter  Maximum number of iterations
  @param[in] im_tols  Tolerances
  @return Value of the minimizer
*/
template<int arg_dim>
Vector<arg_dim> bfgs_mwwp(const std::function<double(const Vector<arg_dim>&)>& obj_fun,
                          const double& obj_fun_lbnd,
                          const Vector<arg_dim>& init_guess,
                          const IterativeMethodTolerances_t& im_tols) {
  const double mwp_delta = 1.0e-2;
  const double mwp_delta1 = 5.0e-3;
  const double damping_scalar = 0.2;

  double fdiff_h;

  //Initial Hessian approximation: identity matrix
  Matrix B(arg_dim, arg_dim, 0.0);
  for (int i = 0; i < arg_dim; i++)
    B[i][i] = 1.0;

  Vector<arg_dim> cur_arg = init_guess, cur_grad, cur_dir, darg, dgrad;
  double dfval, cur_fval = obj_fun(cur_arg);
  fwd_diff_grad<arg_dim>(obj_fun, cur_arg, cur_fval, cur_grad, fdiff_h);

  for (int i = 0; i < im_tols.max_num_iter; i++) {
    // Forward differences are O(fdiff_h) accurate
    if (cur_grad.max_norm()*fdiff_h/std::max(1.0, std::fabs(cur_fval)) <=
        std::numeric_limits<double>::epsilon())
      break;

    Wonton::solve<arg_dim>(B, -cur_grad, cur_dir);

    double df0 = dot(cur_grad, cur_dir);
    // Function is flat in the linesearch direction
    if (df0 >= -std::numeric_limits<double>::epsilon())
      break;

    double initial_alpha = 1.0;
    if (i > 0) initial_alpha = std::min(1.0,
      -2*std::max(-dfval, 10*std::numeric_limits<double>::epsilon())/df0);

    dfval = -cur_fval;
    dgrad = -cur_grad;

    double alpha = mwwp_linesearch<arg_dim>(obj_fun, obj_fun_lbnd,
                                            cur_arg, cur_dir, initial_alpha,
                                            mwp_delta, mwp_delta1,
                                            im_tols, cur_fval, cur_grad);
    dfval += cur_fval;
    darg = alpha*cur_dir;

    if (darg.norm() < im_tols.arg_eps)
      break;

    cur_arg += darg;
    if (std::fabs(cur_fval) <= im_tols.fun_eps)
      return cur_arg;

    dgrad += cur_grad;

    // We update B if the change in gradient is above
    // forward differences error
    if (dgrad.max_norm() > fdiff_h) {
      Vector<arg_dim> B_darg = B*darg;
      if (B_darg.max_norm() <= fdiff_h)
        continue;

      double darg_B_darg = dot(darg, B_darg);
      double darg_dgrad = dot(darg, dgrad);

      // Perform damping if necessary
      if (darg_dgrad < damping_scalar*darg_B_darg) {
        double delta = (1.0 - damping_scalar)*darg_B_darg/(darg_B_darg - darg_dgrad);
        dgrad *= delta;
        dgrad += (1.0 - delta)*B_darg;
        darg_dgrad = dot(darg, dgrad);
      }

      if (darg_dgrad <= std::numeric_limits<double>::epsilon())
        continue;

      // Correction to the initial Hessian guess
      if (i == 0) {
        B *= dgrad.norm(false)/darg_dgrad;
        B_darg = B*darg;
        darg_B_darg = dot(darg, B_darg);
      }

      B += (1.0/darg_dgrad)*(dgrad*dgrad) - (1.0/darg_B_darg)*(B_darg*B_darg);
    }
  }

  return cur_arg;
}

/*!
  @brief D-BFGS algorithm: finds the minimizer of the objective function
  Based on "Improved Damped Quasi-Newton Methods for Unconstrained Optimization"
  by Mehiddin Al-Baali and Lucio Grandinetti,
  Pacific Journal of Optimization (to appear).

  @param[in] obj_fun  Objective function
  @param[in] init_guess  Initial value of the argument
  @param[in] max_num_iter  Maximum number of iterations
  @param[in] im_tols  Tolerances
  @return Value of the minimizer
*/
template<int arg_dim>
Vector<arg_dim> dbfgs(const std::function<double(const Vector<arg_dim>&)>& obj_fun,
                      const double& obj_fun_lbnd,
                      const Vector<arg_dim>& init_guess,
                      const IterativeMethodTolerances_t& im_tols) {
  const double linesearch_c1 = 1.0e-4;
  const double linesearch_c2 = 0.9;

  double fdiff_h;

  //Initial Hessian approximation: identity matrix
  Matrix B(arg_dim, arg_dim, 0.0);
  for (int i = 0; i < arg_dim; i++)
    B[i][i] = 1.0;

  Vector<arg_dim> cur_arg = init_guess, prev_grad, cur_grad, cur_dir, darg, dgrad;
  double dfval, cur_fval = obj_fun(cur_arg);
  fwd_diff_grad<arg_dim>(obj_fun, cur_arg, cur_fval, cur_grad, fdiff_h);

  for (int i = 0; i < im_tols.max_num_iter; i++) {
    // Forward differences are O(fdiff_h) accurate
    if (cur_grad.max_norm()*fdiff_h/std::max(1.0, std::fabs(cur_fval)) <=
        std::numeric_limits<double>::epsilon())
      break;

    Wonton::solve<arg_dim>(B, -cur_grad, cur_dir);

    double df0 = dot(cur_grad, cur_dir);
    // Function is flat in the linesearch direction
    if (df0 >= -std::numeric_limits<double>::epsilon())
      break;

    double initial_alpha = 1.0;
    if (i > 0) initial_alpha = std::min(1.0,
      -2*std::max(-dfval, 10*std::numeric_limits<double>::epsilon())/df0);

    dfval = -cur_fval;
    prev_grad = cur_grad;

    double alpha = linesearch<arg_dim>(obj_fun, obj_fun_lbnd,
                                       cur_arg, cur_dir, initial_alpha,
                                       linesearch_c1, linesearch_c2,
                                       im_tols, cur_fval, cur_grad);
    dfval += cur_fval;
    darg = alpha*cur_dir;
    if (darg.norm() < im_tols.arg_eps)
      break;

    cur_arg += darg;
    if (std::fabs(cur_fval) <= im_tols.fun_eps)
      return cur_arg;

    dgrad = cur_grad - prev_grad;

    // We update B if the change in gradient is above
    // forward differences error
    if (dgrad.max_norm() > fdiff_h) {
      Vector<arg_dim> B_darg = B*darg;
      if (B_darg.max_norm() <= fdiff_h)
        continue;

      double darg_B_darg = dot(darg, B_darg);
      double darg_dgrad = dot(darg, dgrad);

      double b_rec = darg_dgrad/darg_B_darg;
      Vector<arg_dim> H_dgrad;
      Wonton::solve<arg_dim>(B, dgrad, H_dgrad);
      double h_rec = darg_dgrad/dot(dgrad, H_dgrad);

      double l = std::min(b_rec, b_rec*h_rec);
      double m = std::max(b_rec, b_rec*h_rec);

      double sigma2 = std::max(1.0 - 1.0/alpha, 0.5);
      double sigma3 = exp(1);

      double phi = 1.0;
      if (l < 1.0 - sigma2)
        phi = sigma2/(1.0 - l);
      else if (m > 1 + sigma3)
        phi = sigma3/(m - 1.0);
      else {
        double a = 1.0/(b_rec*h_rec);
        if ( (l >= 1.0 - sigma2) && (m <= 1.0 + sigma3) && (a > sigma3) )
          phi = sqrt(sigma3/a);
      }

      // Perform damping if necessary
      if (phi != 1.0) {
        dgrad *= phi;
        dgrad += (1.0 - phi)*B_darg;
        darg_dgrad = dot(darg, dgrad);
      }

      if (darg_dgrad <= std::numeric_limits<double>::epsilon())
        continue;

      // Correction to the initial Hessian guess
      if (i == 0) {
        B *= dgrad.norm(false)/darg_dgrad;
        B_darg = B*darg;
        darg_B_darg = dot(darg, B_darg);
      }

      B += (1.0/darg_dgrad)*(dgrad*dgrad) - (1.0/darg_B_darg)*(B_darg*B_darg);
    }
  }

  return cur_arg;
}

}  // namespace Tangram

#endif  // TANGRAM_SUPPORT_BFGS_H_
