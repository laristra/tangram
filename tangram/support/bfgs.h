/*
This file is part of the Ristra tangram project.
Please see the license file at the root of this repository, or at:
    https://github.com/laristra/tangram/blob/master/LICENSE
*/

#ifndef TANGRAM_BFGS_H_
#define TANGRAM_BFGS_H_

#include <functional>
#include "tangram/support/tangram.h"
#include "tangram/support/Matrix.h"

namespace Tangram {


/*!
  @brief Computes the minimizer of a quadratic interpolant in the interval
  [ maxarg*safeguards[0], maxarg*safeguards[1] ]
  @param[in] maxarg  End of the interval for argument's value
  @param[in] fval0  Value of the function at 0
  @param[in] fval_maxarg  Value of the function for maxarg
  @param[in] df0  Value of the derivative at 0
  @param[in] safeguards  Two-element array to keep the minimizer away
  from the ends of [0, argmax] interval
  @return Value of the minimizer
*/
double quad_interpolant_minimizer(const double maxarg,
                                  const double fval0,
                                  const double fval_maxarg,
                                  const double df0,
                                  const double* safeguards) {
  double numerator = df0*maxarg;
  double denominator = -2.0*(fval_maxarg - fval0 - numerator);
  double argmin, lbnd = maxarg*safeguards[0];

  if (std::fabs(denominator) < std::numeric_limits<double>::epsilon())
    argmin = lbnd;
  else {
    argmin = numerator / denominator;
    if (argmin < lbnd) argmin = lbnd;
    double ubnd = maxarg*safeguards[1];
    if (argmin > ubnd) argmin = ubnd;
  }

  return argmin;
}

/*!
  @brief Linesearch algorithm: finds the argument satisfying Wolfe conditions
  @param[in] obj_fun  Objective function
  @param[in] of_deriv  Derivative of the objective function
  @param[in] init_guess  Initial value of the argument
  @param[in] c1  Value of c1 parameter in Wolfe conditions
  @param[in] c2  Value of c2 parameter in Wolfe conditions
  @param[in] max_num_iter  Maximum number of iterations
  @param[in] arg_eps  Argument tolerance
  @return Value of the minimizer
*/
double linesearch(const std::function<double(const double)>& obj_fun,
                  const std::function<double(const double)>& of_deriv,
                  const double init_guess,
                  const double c1,
                  const double c2,
                  const int max_num_iter,
                  const double arg_eps) {
  //Bracketing step, solution will be in (arg_bnd[0], arg_bnd[1])
  double fval0 = obj_fun(0.0);
  double df0 = of_deriv(0.0);
  
  double arg_max = (std::fabs(df0) > std::numeric_limits<double>::epsilon()) ?
    -fval0/(c1*df0) : DBL_MAX;
  double arg_bnd[2] = {0.0, init_guess};
  double bnd_fval[2]; bnd_fval[0] = fval0;
  double df_lbnd = df0;
  for (int i = 0; i < max_num_iter; i++) {
    bnd_fval[1] = obj_fun(arg_bnd[1]);

    if ( (bnd_fval[1] > fval0 + c1*df0*arg_bnd[1]) || 
         ((bnd_fval[1] >= bnd_fval[0]) && (i > 0)) )
      break;

    double df_ubnd = of_deriv(arg_bnd[1]);
    if (std::fabs(df_ubnd) <= -c2*df0) 
      return arg_bnd[1];

    df_lbnd = df_ubnd;

    if (df_ubnd >= 0.0) {
        std::reverse(std::begin(arg_bnd), std::end(arg_bnd));
        std::reverse(std::begin(bnd_fval), std::end(bnd_fval));
        break;
    }
    
    arg_bnd[0] = arg_bnd[1];  
    bnd_fval[0] = bnd_fval[1];
    arg_bnd[1] *= 2.0;
    if (arg_bnd[1] > arg_max) {
      arg_bnd[1] = arg_max + std::numeric_limits<double>::epsilon();
      bnd_fval[1] = obj_fun(arg_bnd[1]);
      break;
    }
  }

  //Sectioning/zooming step
  double safeguards[2] = {0.2, 0.5};
  for (int i = 0; i < max_num_iter; i++) {
    if (std::fabs(arg_bnd[1] - arg_bnd[0]) < arg_eps)
      return 0.5*(arg_bnd[0] + arg_bnd[1]);

    double cur_arg = quad_interpolant_minimizer(arg_bnd[1] - arg_bnd[0], bnd_fval[0], 
                                                bnd_fval[1], df_lbnd, safeguards) + 
                     arg_bnd[0];

    double cur_fval = obj_fun(cur_arg);
    if ( (cur_fval > fval0 + c1*df0*cur_arg) || (cur_fval >= bnd_fval[0]) ) {
      arg_bnd[1] = cur_arg;
      bnd_fval[1] = cur_fval;
      continue;
    }

    double cur_df = of_deriv(cur_fval);
    if (std::fabs(cur_df) <= -c2*df0)
      return cur_arg;

    if (cur_df*(arg_bnd[1] - arg_bnd[0]) >= 0) {
      arg_bnd[1] = arg_bnd[0];
      bnd_fval[1] = bnd_fval[0];
    }

    arg_bnd[0] = cur_arg;
    bnd_fval[0] = cur_fval;
    df_lbnd = cur_df;
  }

  return 0.5*(arg_bnd[0] + arg_bnd[1]);
}

/*!
  @brief Central finite difference approximation of a gradient
  @param[in] fun  Function for which we compute the gradient
  @param[in] arg_val  We approximate the gradient 
  for this value of the argument
  @param[in] fdiff_h  Finite difference step
  @param[out] approx_grad  Computed approximation of the gradient at arg_val
*/
template<int arg_dim>
void cen_diff_grad(const std::function<double(const Vector<arg_dim>&)>& fun,
                   const Vector<arg_dim>& arg_val,
                   const double fdiff_h,
                   Vector<arg_dim>& approx_grad) {
  Vector<arg_dim> cur_arg = arg_val;
  for (int i = 0; i < arg_dim; i++) {
    cur_arg[i] += fdiff_h;
    double forward_fval = fun(cur_arg);
    cur_arg[i] = arg_val[i] - fdiff_h;
    double backward_fval = fun(cur_arg);
    cur_arg[i] = arg_val[i];

    approx_grad[i] = 0.5*(forward_fval - backward_fval)/fdiff_h;
  }
}

/*!
  @brief BFGS algorithm: finds the minimizer of the objective function
  @param[in] obj_fun  Objective function
  @param[in] init_guess  Initial value of the argument
  @param[in] max_num_iter  Maximum number of iterations
  @param[in] arg_eps  Argument tolerance
  @param[in] grad_eps  Gradient tolerance
  @return Value of the minimizer
*/
template<int arg_dim>
Vector<arg_dim> bfgs(const std::function<double(const Vector<arg_dim>&)>& obj_fun,
                     const Vector<arg_dim>& init_guess,
                     const int max_num_iter,
                     const double arg_eps,
                     const double grad_eps) {
  const double linesearch_c1 = 1.0e-4;
  const double linesearch_c2 = 0.9;
  const double linesearch_init_guess = 1.0;
  const double fdiff_h = pow(std::numeric_limits<double>::epsilon(), 1.0/3.0);
  const double damping_scalar = 0.2;

  //Initial Hessian approximation: identity matrix
  Matrix B_inv(arg_dim, arg_dim, 0.0);
  for (int i = 0; i < arg_dim; i++)
    B_inv[i][i] = 1.0;

  Vector<arg_dim> cur_arg = init_guess, cur_grad, cur_dir, darg, dgrad;
  cen_diff_grad<arg_dim>(obj_fun, cur_arg, fdiff_h, cur_grad);

  for (int i = 0; i < max_num_iter; i++) {    
    if (cur_grad.norm() < grad_eps)
      break;

    cur_dir = -(B_inv*cur_grad);
    
    std::function<double(const double)> linesearch_obj_fun = 
      [&obj_fun, &cur_arg, &cur_dir]
      (const double stepsize)->double {
      return obj_fun(cur_arg + stepsize*cur_dir);
    };

    std::function<double(const double)> linesearch_of_deriv = 
      [&obj_fun, &cur_arg, &linesearch_obj_fun, &fdiff_h]
      (const double stepsize)->double {
      double forward_fval = linesearch_obj_fun(stepsize + fdiff_h);
      double backward_fval = linesearch_obj_fun(stepsize - fdiff_h);
      return 0.5*(forward_fval - backward_fval)/fdiff_h;
    };

    double alpha = linesearch(linesearch_obj_fun, linesearch_of_deriv,
                              linesearch_init_guess, linesearch_c1, 
                              linesearch_c2, max_num_iter, arg_eps);
    darg = alpha*cur_dir;
    if (darg.norm() < arg_eps)
      break;

    cur_arg += darg;
    dgrad = -cur_grad;
    cen_diff_grad<arg_dim>(obj_fun, cur_arg, fdiff_h, cur_grad);
    dgrad += cur_grad;

    // We only update Hessian if gradient has changed
    if (dgrad.norm() > grad_eps) {
      double darg_dgrad = dot(darg, dgrad);
      // Correction to the initial Hessian guess
      if (i == 0)
        B_inv *= darg_dgrad/dgrad.norm(false);

      Vector<arg_dim> Binv_dgrad = B_inv*dgrad;
      double dgrad_Binv_dgrad = dot(dgrad, Binv_dgrad);

      // Perform damping if necessary
      if ((darg_dgrad < damping_scalar*dgrad_Binv_dgrad) && 
          (std::fabs(dgrad_Binv_dgrad - darg_dgrad) > 
           std::numeric_limits<double>::epsilon())) {
        double delta = (1.0 - damping_scalar)*dgrad_Binv_dgrad/
                       (dgrad_Binv_dgrad - darg_dgrad);
        darg *= delta;
        darg += (1.0 - delta)*Binv_dgrad;
        darg_dgrad = dot(darg, dgrad);
      }
      
      B_inv += (1.0/darg_dgrad) * (
        ( (darg_dgrad + dgrad_Binv_dgrad)/darg_dgrad )*(darg*darg) -
        Binv_dgrad*darg - darg*Binv_dgrad 
        );
    }  
  }

  return cur_arg;
}

}  // namespace Tangram

#endif