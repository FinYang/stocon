library(tidyverse)
library(glmnet)
library(boot)
library(quadprog)


# big loop
# for(bl in 1:1){

compar <- function(defaultsetting = T){
  v <- numeric(4)
  rt <- sim_simple(Tn=10, N=5, M=1e4, defaultsetting)

  qp_weights <- qp_weights(rt)
  qp_result <- optim_dynam(Rt = rt, weights = qp_weights)
  v[1] <- qp_result$v

  qp_l_weights <- lasso_weights(rt, qp_lasso = TRUE, qp_weights = qp_weights)
  qp_l_result <- optim_dynam(Rt = rt, weights = qp_l_weights)
  v[2] <- qp_l_result$v

  qp_cv_weights <- qp_weights(rt, constr = "cv")
  qp_cv_result <- optim_dynam(Rt = rt, weights = qp_cv_weights)
  v[3] <- qp_cv_result$v

  qp_nc_weights <- qp_weights(rt, constr = NULL)
  qp_nc_result <- optim_dynam(Rt = rt, weights = qp_nc_weights)
  v[4] <- qp_nc_result$v
  v
}

rt <- sim_simple(Tn=10, N=5, M=1e4, defaultsetting)

path <- function(M, rt){
  rt
  qp_weights <- qp_weights(rt)
  qp_result <- optim_dynam(Rt = rt, weights = qp_weights)
  v[1] <- qp_result$v

  qp_l_weights <- lasso_weights(rt, qp_lasso = TRUE, qp_weights = qp_weights)
  qp_l_result <- optim_dynam(Rt = rt, weights = qp_l_weights)
  v[2] <- qp_l_result$v

  qp_cv_weights <- qp_weights(rt, constr = "cv")
  qp_cv_result <- optim_dynam(Rt = rt, weights = qp_cv_weights)
  v[3] <- qp_cv_result$v

  qp_nc_weights <- qp_weights(rt, constr = NULL)
  qp_nc_result <- optim_dynam(Rt = rt, weights = qp_nc_weights)
  v[4] <- qp_nc_result$v
  v
}


