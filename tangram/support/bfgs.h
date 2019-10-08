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
inline
double quad_interpolant_minimizer(const double arg0,
                                  const double arg1,
                                  const double fval0,
                                  const double fval1,
                                  const double df0,
                                  const double* safeguards) {
  double darg = arg1 - arg0;
  double denominator = fval1 - fval0 - df0*(darg);
  std::vector<double> denominator_terms = {fval0, fval1, df0*(darg)};
  double max_denominator_term = *std::max(denominator_terms.begin(), denominator_terms.end());
  double arg_min;
  double bnd0 = arg0 + safeguards[0]*darg,
         bnd1 = arg0 + safeguards[1]*darg;

  if (std::fabs(denominator/max_denominator_term) <= 
      std::numeric_limits<double>::epsilon())
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
  @brief Central finite difference approximation of a gradient
  @param[in] fun  Function for which we compute the gradient
  @param[in] arg_val  We approximate the gradient
  for this value of the argument
  @param[out] approx_grad  Computed approximation of the gradient at arg_val
  @param[out] fdiff_h  Finite difference step
*/
template<int arg_dim>
void cen_diff_grad(const std::function<double(const Vector<arg_dim>&)>& fun,
                   const Vector<arg_dim>& arg_val,
                   Vector<arg_dim>& approx_grad,
                   double& fdiff_h) {
  fdiff_h = pow(std::numeric_limits<double>::epsilon(), 1.0/3.0);

  Vector<arg_dim> cur_arg = arg_val;
  for (int i = 0; i < arg_dim; i++) {
    cur_arg[i] += fdiff_h;
    double forward_fval = fun(cur_arg);
    cur_arg[i] = arg_val[i] - fdiff_h;
    double backward_fval = fun(cur_arg);
    cur_arg[i] = arg_val[i];

    approx_grad[i] = (forward_fval - backward_fval)/(2*fdiff_h);
  }
}

/*!
  @brief Finite difference approximation of a gradient. Switches from forward to
  central finite differences if the value of the gradient is not significant
  relative to the method's error.
  @param[in] fun  Function for which we compute the gradient
  @param[in] arg_val  We approximate the gradient
  for this value of the argument
  @param[in] fun_val  Value of function at arg_val
  @param[in] rel_err_bnd  If the estimate of relative error in gradient is above this,
  we switch from forward to central finite differences
  @param[out] approx_grad  Computed approximation of the gradient at arg_val
*/
template<int arg_dim>
void finite_diff_grad(const std::function<double(const Vector<arg_dim>&)>& fun,
                      const Vector<arg_dim>& arg_val,
                      const double& fun_val,
                      const double rel_err_bnd,
                      Vector<arg_dim>& approx_grad) {
  double fdiff_h;
  fwd_diff_grad(fun, arg_val, fun_val, approx_grad, fdiff_h);

  double max_norm = approx_grad.max_norm();
  if (fdiff_h*max_norm/fun_val < std::numeric_limits<double>::epsilon()) {
    approx_grad.zero();
    return;
  }
  if (fdiff_h/max_norm <= rel_err_bnd) return;

  //Fallback to central differences
  cen_diff_grad(fun, arg_val, approx_grad, fdiff_h);
  if (2*fdiff_h*max_norm/fun_val < std::numeric_limits<double>::epsilon())
    approx_grad.zero();
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
  @param[in] rel_err_bnd  If the estimate of relative error in gradient is above this,
  we switch from forward to central finite differences  
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
                  const double rel_err_bnd,
                  double& cur_fval,
                  Vector<arg_dim>& cur_grad) {
  // Bracketing step, solution will be in (alpha_bnd[0], alpha_bnd[1])
  double fval0 = cur_fval;
  double df0 = dot(cur_grad, dir);
  assert(df0 < 0.0);
  double alpha_tol = im_tols.arg_eps/dir.norm();

  auto fval = [&obj_fun, &arg, &dir](double alpha) {
    return obj_fun(arg + alpha*dir);
  };

  // alpha_max is chosen so that it should violate the sufficient decrease condition
  // and therefore provides an acceptable upper bound for the interval
  double alpha_max = (obj_fun_lbnd - fval0)/(c1*df0);

  double alpha_bnd[2] = {0.0, initial_guess};
  double bnd_fval[2]; bnd_fval[0] = fval0;
  Vector<arg_dim>& lbnd_grad = cur_grad;
  double df_lbnd = df0;
  // At the beginning of each loop, cur_fval and cur_grad are at alpha_bnd[0]
  for (int i = 0; i < im_tols.max_num_iter; i++) {
    double tau = 9.0; // interval length scaling

    bnd_fval[1] = fval(alpha_bnd[1]);

    if (bnd_fval[1] - obj_fun_lbnd < im_tols.fun_eps) {
      cur_fval = bnd_fval[1];
      return alpha_bnd[1];
    }

    if ( (bnd_fval[1] > fval0 + c1*alpha_bnd[1]*df0) ||
         ((bnd_fval[1] >= bnd_fval[0]) && (i > 0)) )
      break;

    cur_fval = bnd_fval[1];
    finite_diff_grad<arg_dim>(obj_fun, arg + alpha_bnd[1]*dir, bnd_fval[1], 
                              rel_err_bnd, cur_grad);
    // alpha_bnd[1] will be the next alpha_bnd[0]
    df_lbnd = dot(cur_grad, dir);

    if (std::fabs(df_lbnd) <= -c2*df0)
      return alpha_bnd[1];

    lbnd_grad = cur_grad;
    if (df_lbnd >= 0.0) {
        std::reverse(std::begin(alpha_bnd), std::end(alpha_bnd));
        std::reverse(std::begin(bnd_fval), std::end(bnd_fval));
        break;
    }

    double dalpha = alpha_bnd[1] - alpha_bnd[0];

    alpha_bnd[0] = alpha_bnd[1];
    bnd_fval[0] = bnd_fval[1];
    alpha_bnd[1] += tau*dalpha;

    if (alpha_bnd[1] > alpha_max - alpha_tol) {
      alpha_bnd[1] = alpha_max;
      bnd_fval[1] = obj_fun(alpha_max);
      break;
    }
  }

  // Sectioning/zooming step
  // At the beginning of each loop, cur_fval and cur_grad are at alpha_bnd[0]
  double safeguards[2] = {0.1, 0.5};
  for (int i = 0; i < im_tols.max_num_iter; i++) {
    if (std::fabs(alpha_bnd[1] - alpha_bnd[0]) < alpha_tol) {
      cur_grad = lbnd_grad;
      return alpha_bnd[0];
    }

    double cur_alpha = quad_interpolant_minimizer(alpha_bnd[0], alpha_bnd[1], bnd_fval[0],
                                                  bnd_fval[1], df_lbnd, safeguards);
    cur_fval = fval(cur_alpha);
    if (cur_fval - obj_fun_lbnd < im_tols.fun_eps)
      return cur_alpha;

    if ( (cur_fval > fval0 + c1*cur_alpha*df0) || (cur_fval >= bnd_fval[0]) ) {
      alpha_bnd[1] = cur_alpha;
      bnd_fval[1] = cur_fval;
      continue;
    }

    finite_diff_grad<arg_dim>(obj_fun, arg + cur_alpha*dir, cur_fval,
                              rel_err_bnd, cur_grad);

    double cur_df = dot(cur_grad, dir);
    if (std::fabs(cur_df) <= -c2*df0)
      return cur_alpha;

    if (cur_df*(alpha_bnd[1] - alpha_bnd[0]) >= 0) {
      alpha_bnd[1] = alpha_bnd[0];
      bnd_fval[1] = bnd_fval[0];
    }

    alpha_bnd[0] = cur_alpha;
    bnd_fval[0] = cur_fval;
    lbnd_grad = cur_grad;
    df_lbnd = cur_df;
  }

  cur_fval = bnd_fval[0];
  cur_grad = lbnd_grad;
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
  const double grad_rel_err_bnd = 0.1;

  //Initial Hessian approximation: identity matrix
  Matrix B(arg_dim, arg_dim, 0.0);
  for (int i = 0; i < arg_dim; i++)
    B[i][i] = 1.0;
  Matrix B_old;

  Vector<arg_dim> cur_arg = init_guess, prev_grad, cur_grad, cur_dir, darg, dgrad;
  double dfval, cur_fval = obj_fun(cur_arg);
  finite_diff_grad<arg_dim>(obj_fun, cur_arg, cur_fval, 
                            grad_rel_err_bnd, cur_grad);

  for (int i = 0; i < im_tols.max_num_iter; i++) {
    // Check if gradient is nonzero, we're finished otherwise
    if (cur_grad.max_norm() == 0.0)
      return cur_arg;

    Wonton::solve<arg_dim>(B, -cur_grad, cur_dir);
    double df0 = dot(cur_grad, cur_dir);

    // Function should be decreasing in the linesearch direction
    if (df0 >= 0.0) {
      if (i > 0) {
        B = B_old;
        Wonton::solve<arg_dim>(B, -cur_grad, cur_dir);
        df0 = dot(cur_grad, cur_dir);
      }
      // Can't decrease the function value any further, finishing
      if (df0 >= 0.0)
        return cur_arg;
    }

    double initial_alpha = 1.0;
    if (i > 0) initial_alpha = std::min(1.0,
      -2*std::max(-dfval, 10*std::numeric_limits<double>::epsilon())/df0);

    dfval = -cur_fval;
    prev_grad = cur_grad;

    double alpha = linesearch<arg_dim>(obj_fun, obj_fun_lbnd,
                                       cur_arg, cur_dir, initial_alpha,
                                       linesearch_c1, linesearch_c2,
                                       im_tols, grad_rel_err_bnd,
                                       cur_fval, cur_grad);
    darg = alpha*cur_dir;
    // Check if the change in argument was above tolerance, finish otherwise
    if (darg.norm() < im_tols.arg_eps) return cur_arg;

    dfval += cur_fval;
    // Check if the change in function was numerically significant, finish otherwise
    if (std::fabs(dfval/cur_fval) <= std::numeric_limits<double>::epsilon())
      return cur_arg;

    cur_arg += darg;
    // Check if the value of the function is above tolerance, finish otherwise
    if (std::fabs(cur_fval) < im_tols.fun_eps) return cur_arg;

    dgrad = cur_grad - prev_grad;

    bool valid_dgrad = false;
    for (int idim = 0; idim < arg_dim; idim++)
      if (dgrad[idim]/cur_grad[idim] > std::numeric_limits<double>::epsilon()) {
        valid_dgrad = true;
        break;
      }

    // We update B if the change in gradient is not too small
    // relative to the value of the gradient
    if (valid_dgrad) {
      B_old = B;

      Vector<arg_dim> B_darg = B*darg;
      double darg_B_darg = dot(darg, B_darg);

      double darg_dgrad = dot(darg, dgrad);
      // Perform damping if necessary
      if (darg_dgrad < damping_scalar*darg_B_darg) {
        double delta = (1.0 - damping_scalar)*darg_B_darg/(darg_B_darg - darg_dgrad);
        dgrad *= delta;
        dgrad += (1.0 - delta)*B_darg;
        darg_dgrad = dot(darg, dgrad);
      }

      if (darg_dgrad <= im_tols.arg_eps*std::numeric_limits<double>::epsilon()) continue;

      // Correction to the initial Hessian guess
      if (i == 0) {
        B *= dgrad.norm(false)/darg_dgrad;
        B_darg = B*darg;
        darg_B_darg = dot(darg, B_darg);
      }

      B += (1.0/darg_dgrad)*(dgrad*dgrad) - (1.0/darg_B_darg)*(B_darg*B_darg);
    }
  }

  // Stopping conditions were NOT reached in the given number of iterations,
  // return the last value of the argument
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
  const double grad_rel_err_bnd = 0.1;

  //Initial Hessian approximation: identity matrix
  Matrix B(arg_dim, arg_dim, 0.0);
  for (int i = 0; i < arg_dim; i++)
    B[i][i] = 1.0;
  Matrix B_old;

  Vector<arg_dim> cur_arg = init_guess, prev_grad, cur_grad, cur_dir, darg, dgrad;
  double dfval, cur_fval = obj_fun(cur_arg);
  finite_diff_grad<arg_dim>(obj_fun, cur_arg, cur_fval, 
                            grad_rel_err_bnd, cur_grad);

  for (int i = 0; i < im_tols.max_num_iter; i++) {
    // Check if gradient is nonzero, we're finished otherwise
    if (cur_grad.max_norm() == 0.0)
      return cur_arg;

    Wonton::solve<arg_dim>(B, -cur_grad, cur_dir);
    double df0 = dot(cur_grad, cur_dir);

    // Function should be decreasing in the linesearch direction
    if (df0 >= 0.0) {
      if (i > 0) {
        B = B_old;
        Wonton::solve<arg_dim>(B, -cur_grad, cur_dir);
        df0 = dot(cur_grad, cur_dir);
      }
      // Can't decrease the function value any further, finishing
      if (df0 >= 0.0)
        return cur_arg;
    }

    double initial_alpha = 1.0;
    if (i > 0) initial_alpha = std::min(1.0,
      -2*std::max(-dfval, 10*std::numeric_limits<double>::epsilon())/df0);

    dfval = -cur_fval;
    prev_grad = cur_grad;

    double alpha = linesearch<arg_dim>(obj_fun, obj_fun_lbnd,
                                       cur_arg, cur_dir, initial_alpha,
                                       linesearch_c1, linesearch_c2,
                                       im_tols, grad_rel_err_bnd,
                                       cur_fval, cur_grad);
    darg = alpha*cur_dir;
    // Check if the change in argument was above tolerance, finish otherwise
    if (darg.norm() < im_tols.arg_eps) return cur_arg;

    dfval += cur_fval;
    // Check if the change in function was numerically significant, finish otherwise
    if (std::fabs(dfval/cur_fval) <= std::numeric_limits<double>::epsilon())
      return cur_arg;

    cur_arg += darg;
    // Check if the value of the function is above tolerance, finish otherwise
    if (std::fabs(cur_fval) < im_tols.fun_eps) return cur_arg;

    dgrad = cur_grad - prev_grad;

    bool valid_dgrad = false;
    for (int idim = 0; idim < arg_dim; idim++)
      if (dgrad[idim]/cur_grad[idim] > std::numeric_limits<double>::epsilon()) {
        valid_dgrad = true;
        break;
      }

    // We update B if the change in gradient is not too small
    // relative to the value of the gradient
    if (valid_dgrad) {
      B_old = B;

      Vector<arg_dim> B_darg = B*darg;
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

      if (darg_dgrad <= im_tols.arg_eps*std::numeric_limits<double>::epsilon()) continue;

      // Correction to the initial Hessian guess
      if (i == 0) {
        B *= dgrad.norm(false)/darg_dgrad;
        B_darg = B*darg;
        darg_B_darg = dot(darg, B_darg);
      }

      B += (1.0/darg_dgrad)*(dgrad*dgrad) - (1.0/darg_B_darg)*(B_darg*B_darg);
    }
  }

  // Stopping conditions were NOT reached in the given number of iterations,
  // return the last value of the argument
  return cur_arg;
}

}  // namespace Tangram

#endif  // TANGRAM_SUPPORT_BFGS_H_
