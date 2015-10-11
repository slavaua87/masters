
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;
using namespace std;

vec forwsub_c(mat A, vec b) {
  vec x = zeros(3);
  for (int i = 0; i < 3; ++i) {
    x[i] = (b[i] - as_scalar(A.row(i) * x)) / A(i, i);
  }
  return x;
}

double normal_logcopula_c(vec cop_value, mat sqrt_rho) {
  vec x(3), marginals(3);
  unsigned i;
  
  for (i = 0; i < 3; ++i) {
    x[i] = R::qnorm(cop_value[i], 0, 1, true, false);
    marginals[i] = R::dnorm(x[i], 0, 1, false);
  }
  uvec test = marginals == 0.0;
  if (any(test)) {
    marginals.elem(find(test)).fill(DBL_EPSILON);
  //numeric_limits<double>::min()
  }
  vec sqrt_rho_x = forwsub_c(trans(sqrt_rho), x);
  vec rho_x = sum(pow(sqrt_rho_x, 2));
  vec rho_diag = sqrt_rho.diag();
  test = rho_diag == 0.0;
  if (any(test)) {
    rho_diag.elem(find(test)).fill(DBL_EPSILON);

  }
  double logdens = as_scalar((-1.5 * log(2 * M_PI) - sum(log(rho_diag)) - 
                              0.5 * rho_x) - sum(log(marginals)));  
  return logdens;
}

double calc_logdensity_c(double delta, double beta, double t_nd,
                         double nu, double eta, double shape1, double shape2,
                         double shape, double scale, mat sqrt_rho,
                         string model) {
//  # Calculates probability density for independent and normal
//  # Input is a numeric scalars delta, beta, t_nd, 
//  # numeric data.frame params, character scalar model
//  # Output is a numeric scalar dens
//  # Notes: params vector has the following order - 1:alpha, 2:nu, 3:eta,
//  # 4:lambda, 5:gamma, 6:chi, 7:phi, 8:rho_db, 9:rho_dt, 10:rho_bt, 11:omega

  vec marginal_dens(3);
  marginal_dens[0] = R::dnorm(delta, nu, eta, false);
  marginal_dens[1] = R::dbeta(beta, shape1, shape2, false);
  marginal_dens[2] = R::dgamma(t_nd, shape, scale, false);
  
  uvec test_bounds = marginal_dens == 0.0;
  if (any(test_bounds)) {
    marginal_dens.elem(find(test_bounds)).fill(DBL_EPSILON);

  }
  bool test_inf = marginal_dens.is_finite();
  if (test_inf) {
    marginal_dens.elem(find_nonfinite(marginal_dens)).fill(
      numeric_limits<double>::max());
  }
  double marginal_logdens = sum(log(marginal_dens));

  if (model == "independent") 
    return marginal_logdens;
  
  vec cop_value(3);
  cop_value[0] = R::pnorm(delta, nu, eta, true, false);
  cop_value[1] = R::pbeta(beta, shape1, shape2, true, false);
  cop_value[2] = R::pgamma(t_nd, shape, scale, true, false);
  
  test_bounds = cop_value == 0.0;
  if (any(test_bounds)) {
    cop_value.elem(find(test_bounds)).fill(DBL_EPSILON);
  }
  test_bounds = cop_value == 1.0;
  if (any(test_bounds)) {
    cop_value.elem(find(test_bounds)).fill(1 - DBL_EPSILON);
  }
  
  double joint_logdens = normal_logcopula_c(cop_value, sqrt_rho);

  return marginal_logdens + joint_logdens;
}

double dwiener_cpp(double q, double alpha, double tau,
                   double beta, double delta, double choice) {
  double kl, ks, ans;
  int k,K;
  double err = 1e-6;

  // choice is 1 for lower, 2 for upper boundary
  if (choice == 2.0) {
    beta = 1 - beta;
    delta = -delta;
  }
  
  q = q - tau; // remove non-decision time from q
  q = q / pow(alpha, 2); // convert t to normalized time tt

  // calculate number of terms needed for large t
  if (M_PI * q * err < 1) { // if error threshold is set low enough
      kl = sqrt(-2 * log(M_PI * q * err) / (pow(M_PI, 2) * q)); // bound
      kl = kl > 1 / (M_PI * sqrt(q)) ? kl : 1 / (M_PI * sqrt(q)); // ensure boundary conditions met
  }
  else { // if error threshold set too high
      kl = 1 / (M_PI * sqrt(q)); // set to boundary condition
  }

  // calculate number of terms needed for small t
  if ((2 * sqrt(2 * M_PI * q) * err) < 1) { // if error threshold is set low enough
      ks = 2 + sqrt(-2 * q * log(2 * sqrt(2 * M_PI * q) * err)); // bound
      ks = ks > sqrt(q) + 1 ? ks : sqrt(q) + 1; // ensure boundary conditions are met
  }
  else { // if error threshold was set too high
      ks = 2; // minimal kappa for that case
  }

  // compute density: f(tt|0,1,beta)
  ans = 0; //initialize density
  if (ks < kl) { // if small t is better (i.e., lambda<0)
      K = ceil(ks); // round to smallest integer meeting error
      for (k = -floor((K - 1) / 2); k <= ceil((K - 1) / 2); k++) { // loop over k
          ans = ans + (beta + 2 * k) * exp(-(pow((beta + 2 * k), 2)) / 2 / q); // increment sum
      }
     // ans = log(ans) - 0.5 * (log(2) + log(M_PI) + log(pow(q, 3))); // add constant term
      ans = ans / sqrt(2 * M_PI * pow(q, 3));
  }
  else { // if large t is better...
      K = ceil(kl); // round to smallest integer meeting error
      for (k = 1; k <= K; k++) {
          ans = ans + k * exp(-(pow(k, 2)) * (pow(M_PI, 2)) * q / 2) * 
            sin(k * M_PI * beta); // increment sum
      }
      //ans = log(ans) + log(M_PI); // add constant term
      ans = ans * M_PI;
  }

// convert to f(t|v,a,w) and return result
//  ans = log(ans + ((-delta * alpha * beta -(pow(delta, 2)) *
//        (q * pow(alpha, 2)) / 2) - log(pow(alpha,2))));
  ans = ans * exp(-delta * alpha * beta - (pow(delta, 2)) * 
    (q * pow(alpha, 2)) / 2) / (pow(alpha, 2));
  if (ans < numeric_limits<double>::min() || ans == datum::inf || isnan(ans)) {
    ans = DBL_EPSILON;
  }
  return ans;
}

double calc_density_integrand_c(const double * x, double rt, double choice,
                                double sigma, double alpha, double nu,
                                double eta, double shape1, double shape2, 
                                double shape, double scale, mat sqrt_rho,
                                string model) {
//  # Purpose: Calculates integrand for a rt density for a given a 
//  # choice and model
//  # Input: numeric vector dimensions, numeric scalars rt, choice, sigma
//  # numeric vector params, character scalar model
//  # Output: numeric scalar integrand
//  # Notes: params vector has the following order - 1:alpha, 2:nu, 3:eta,
//  # 4:lambda, 5:gamma, 6:chi, 8:phi, 9:rho_db, 10:rho_dt, 11:rho_bt, 12:omega
  
  
  double delta_trans = *x;
  double beta = *(x + 1);
  double t_nd = *(x + 2);
  double delta = delta_trans / (1 - pow(delta_trans, 2));
    
  double behav_dens = dwiener_cpp(rt, alpha / sigma, t_nd, 
                                  beta, delta / sigma, choice);
  double param_logdens = calc_logdensity_c(delta, beta, t_nd,     
                                           nu, eta, shape1, shape2,
                                           shape, scale, sqrt_rho, model);
  
  double jacobian = (1 + pow(delta_trans, 2)) / pow(1 - pow(delta_trans, 2), 2);
  
  double integrand = exp(log(behav_dens) + param_logdens + log(jacobian));
  
  if (isnan(integrand) || integrand == datum::inf) 
    integrand = DBL_EPSILON;
      
  return integrand;
}








