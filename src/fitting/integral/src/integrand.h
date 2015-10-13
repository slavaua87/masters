#ifndef INTEGRAND_H
#define INTEGRAND_H

#include <RcppArmadillo.h>
#include <string>
using namespace Rcpp;
using namespace std;
using namespace arma;

// forward substitution 
vec forwsub_c(mat A, vec b);

// normal copula density
double normal_logcopula_c(NumericVector cop_value, mat sqrt_rho);

// parameter density
double calc_logdensity_c(double delta, double beta, double t_nd,
                         double nu, double eta, double shape1, double shape2,
                         double shape, double scale, mat sqrt_rho,
                         string model);

// wiener density
double dwiener_cpp(double q, double alpha, double tau,
                   double beta, double delta, string response);

// integrand
double calc_density_integrand_c(const double * x, double rt, double choice,
                                double sigma, double alpha, double nu,
                                double eta, double shape1, double shape2, 
                                double shape, double scale, mat sqrt_rho,
                                string model);

#endif
