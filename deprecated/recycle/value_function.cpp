#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

// [[Rcpp::export(name=".test")]]
double beta_function(NumericVector para_beta, double w){
  NumericVector beta = clone(para_beta);
  double out = beta[0] +  beta[1] * w + beta[2] * w * w;
  return out;
}

double c_function(NumericVector para_c, double w){
  return w*exp(para_c[0]+para_c[1]*w)/(1+exp(para_c[0]+para_c[1]*w));
}


// [[Rcpp::export(name=".varmean")]]
List value_varmean(NumericMatrix para,
                   int t,
                   NumericMatrix Rr,
                   double Rf,
                   int M,
                   int Tn,
                   NumericMatrix W,
                   double discount,
                   double lambda){
  NumericMatrix para_beta = para(_,Rcpp::Range(0,2));
  NumericMatrix para_c = para(_,Rcpp::Range(3, 4));
  NumericVector value(M);
  double final_value = 0;
  // NumericMatrix R(M,1);
  // double BETA = 0;
  // double C = 0;
  NumericMatrix BETA(M,Tn);
  NumericMatrix C(M,Tn);

  NumericMatrix w = clone(W);
  if(t < Tn){
    // for(int s=t+1; s < (Tn+1); s++){
    for(int s=t+1; s < (Tn+1); s++){
      // NumericVector R = Rr[s-1];
      for(int i=0; i<M; i++){
        BETA(i,s-1) = beta_function(para_beta(s-1,_), w(i,s-1));
        C(i,s-1) = c_function(para_c(s-1,_), w(i,s-1));
        w(i, s) = (w(i, s-1)-BETA(i,s-1))*Rf + BETA(i,s-1)*Rr(i,s) - C(i,s-1);
        value[i] += pow(discount, s-t) * (pow(C(i,s-1),2)-2*lambda*C(i,s-1));
      }
    }
  }
  for(int i=0; i<M;i++){
    value[i] += pow(discount, Tn-t)*(pow(W(i,Tn), 2)-2*lambda*W(i,Tn));
    final_value += value[i];
  }




  return List::create(_["W"] = w, _["value"] = final_value/M, BETA, C, para_beta);
}







/***R

varmean_cpp <- function(update_par = NULL, i = 0, para, Rr, Rf,
                        M = NROW(Rr[[1]]), Tn = length(Rr), W, discount=1/1.01,
                        lambda = 1/2, returnW = FALSE){
  if(!is.null(update_par)) para[i, ] <- update_par
  Rr <- do.call(cbind, Rr)
  out <- .varmean(para = para, t=i, Rr = Rr, Rf = Rf,
                  M = M, Tn = Tn,W = W, discount = discount, lambda = lambda)
  if(returnW) return(out$W)
  return(out$value)
  # return(out)
}
#
# Tn <- 5
# N <- 5
# M <- 10
#
# Rf <- 1.01
# discount <- 1/1.05
# W <- matrix(nrow = M, ncol = Tn+1)
# W[,1] <- 1000
# pm <- NULL
# lambda = 1/2
# discount=1/1.01
# # Rr <- simulate_riskys()
# # Rr <- step2(Rr)
# mean_Rr <- 1.05
# sd_Rr <- 0.03
# Rr <- replicate(Tn,  matrix(rnorm(M, mean_Rr, sd_Rr), ncol = 1), simplify = F)
#
# para <- matrix(0,nrow = Tn, ncol = 5)
#
# value_varmean(para = para, Rr = Rt, Rf = Rf, W = W,
#               M = M, Tn = Tn, discount = discount, lambda = lambda,
#               returnW = T)
# varmean_cpp(para = para, Rr = Rr, Rf = Rf, W = W,
#             M = M, Tn = Tn, discount = discount, lambda = lambda,
#             returnW = T)
# #
# bench::mark(value_varmean(para = para, Rr = Rr, Rf = Rf, W = W,
#                           M = M, Tn = Tn, discount = discount, lambda = lambda,
#                           returnW = TRUE),
#             varmean_cpp(para = para, Rr = Rr, Rf = Rf, W = W,
#                         M = M, Tn = Tn, discount = discount, lambda = lambda,
#                         returnW = T), check = F) %>% View
#
#
# value_varmean(update_par = c(0,2,1.3, 2, 0.4), para = para,i=Tn-1, Rr = Rr, Rf = Rf, W = W,
#               M = M, Tn = Tn, discount = discount, lambda = lambda,
#               returnW = F)

*/

