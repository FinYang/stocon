library(tidyverse)
library(glmnet)
library(boot)
library(quadprog)


# big loop
# for(bl in 1:1){

rt <- sim_simple(Tn=10, N=5, M=1e4)
weights <- lasso_weights(rt, qp_lasso = TRUE)
theta_v <- optim_dynam(Rt = rt, weights = weights)







# }
