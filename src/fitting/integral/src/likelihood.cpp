#include <limits>
#include <RcppArmadillo.h>
#include <omp.h>
#include "integrand.h"
#include "cubature.h" 

using namespace Rcpp;
using namespace std;
using namespace arma;

struct params {
  double rt;
  double choice;
  double sigma;
  double alpha;
  double nu;
  double eta;
  double shape1;
  double shape2;  
  double shape; 
  double scale; 
  mat sqrt_rho;
  string model;
  unsigned * counts;
};

double calc_prior_cpp(vec theta, string model, double omega = 1) {
  // Purpose: evaluates prior density of unknown parameters
  
  double prior_density = 0;
  unsigned i;
  for (i = 0; i < 20; ++i) 
    prior_density += R::dnorm(theta[i], 0, 100, true);
  if (model == "normal") 
    prior_density += log(pow((1 - pow(theta[20], 2) - pow(theta[21], 2) -
                             pow(theta[22], 2) + 2 * theta[20] * theta[21] * 
                             theta[22]), omega - 1.0));
  return prior_density;
}

int fwrap_v(unsigned ndim, size_t npts, const double * x, void * fdata,         
            unsigned fdim, double * fval) {
  // Purpose: wraps integrand to fit with hcubature algorithm

  params mydata = *((params *) fdata);
  vec res(fval, npts, false);
  unsigned i;
  for (i = 0; i < npts; ++i) {
    double y[3]; 
    y[0] = x[i * 3];
    y[1] = x[i * 3 + 1];
    y[2] = x[i * 3 + 2];
    
    res[i] = calc_density_integrand_c(y, mydata.rt,     
                      			          mydata.choice,
                    				      	  mydata.sigma,
                    				      	  mydata.alpha,
                    				      	  mydata.nu,
                    				      	  mydata.eta,
                    				      	  mydata.shape1,
                    				      	  mydata.shape2,
                    				      	  mydata.shape,
                    				      	  mydata.scale,
                    				      	  mydata.sqrt_rho,
                    				      	  mydata.model);
  }
  *mydata.counts += npts;
  return 0;
}

double integrate_density( double rt,  double choice,  double sigma,
                          double alpha,  double nu,  double eta,
                          double shape1,  double shape2,  double shape,
                          double scale,  mat sqrt_rho,  string model,
                          double tol,  size_t maxEval) {
  // Purpose: calculates tripple integral of behavioral density
  //          with respect to density of parameters
  
  double integral;
  double error;
  unsigned counts;
  params mydata;

  double * integral_pt = &integral;
  double * error_pt = &error;
  
  mydata.rt = rt;     
  mydata.choice = choice;
  mydata.sigma = sigma;
  mydata.alpha = alpha;
  mydata.nu = nu;
  mydata.eta = eta;
  mydata.shape1 = shape1;
  mydata.shape2 = shape2;
  mydata.shape = shape;
  mydata.scale = scale;
  mydata.sqrt_rho = sqrt_rho;
  mydata.model = model;
  mydata.counts = &counts;

  const double xmin[3] = {-1 + DBL_EPSILON,
    	                    numeric_limits<double>::min(),
			                    numeric_limits<double>::min()},
               xmax[3] = {1 - DBL_EPSILON,
			                    1 - DBL_EPSILON,
			                    rt};

  do {
    counts = 0;
    hcubature_v(1, fwrap_v, &mydata, 3, xmin, xmax,
                maxEval, 0, tol, ERROR_L1, integral_pt, error_pt);
    tol *= 0.1;
  } while (counts < 2e3 || tol < 1e-7);
  
  return integral;
}

vec calc_likelihood_cpp(mat data_mat, string model, size_t thread_n,
                        size_t chunk_n, double tol, size_t maxEvals) {
  // Purpose: evalutes likelihood function of response times and responses
  
  mat sqrt_rho;
  if (model == "normal") {
    double rho_db = data_mat(0, 9);
    double rho_dt = data_mat(0, 10);
    double rho_bt = data_mat(0, 11);
    sqrt_rho << 1 << rho_db << rho_dt << endr
             << rho_db << 1 << rho_bt << endr
	           << rho_dt << rho_bt << 1 << endr;
    sqrt_rho = chol(sqrt_rho);
 }

  unsigned i, n;
  n = data_mat.n_rows;
  vec behav_density(n);
 
  vec rt(data_mat.colptr(0), n, false);
  vec choice(data_mat.colptr(1), n, false);
  vec alpha(data_mat.colptr(2), n, false);
  vec nu(data_mat.colptr(3), n, false);
  vec eta(data_mat.colptr(4), n, false); 
  vec shape1(data_mat.colptr(5), n, false); 
  vec shape2(data_mat.colptr(6), n, false);
  vec shape(data_mat.colptr(7), n, false);
  vec scale(data_mat.colptr(8), n, false);

  #pragma omp parallel for \
  schedule(dynamic, chunk_n) \
  num_threads(thread_n)

  for (i = 0; i < n; ++i) {
    behav_density[i] = integrate_density(rt[i], choice[i], 0.1, 
                                         alpha[i], nu[i], eta[i], 
    		                                 shape1[i], shape2[i], shape[i],
                                         scale[i], sqrt_rho, "normal",
				                                 tol, maxEvals);  
  } 

  uvec zero_test = behav_density <= 0;
  if (any(zero_test))
    behav_density.elem(find(zero_test)).fill(DBL_EPSILON);
  return behav_density; 
}

double joint_logdensity_cpp(mat data_mat, vec theta_prop, string model,
                            unsigned thread_n, unsigned chunk_n, 
                            double tol, unsigned maxEvals) {
  // Purpose: evalutes joint density of behavior and parameters
  
  double logprior = calc_prior_cpp(theta_prop, model);
  vec likelihood = calc_likelihood_cpp(data_mat, model, thread_n,
                                       chunk_n, tol, maxEvals);
  unsigned i, n;
  n = likelihood.n_elem;
  double loglikelihood = 0;

  for (i = 0; i < n; ++i) {
    loglikelihood += log(likelihood[i]);
  }
  return logprior + loglikelihood; 
}

RCPP_MODULE(integral) {
  Rcpp::function("joint_logdensity_cpp",
                 &joint_logdensity_cpp,
                 "calculates joint density");
}



