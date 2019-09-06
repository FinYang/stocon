#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export(name=".test")]]
double beta_function(NumericVector para_beta, double w){
  // NumericVector beta = clone(para_beta);
  // double out = beta[1] +  beta[2] * w + beta[3] * w * w;
  double out = w * w;
  return out;
}

/***R
.test(1:3, 100)
*/
