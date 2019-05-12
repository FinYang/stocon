library(tidyverse)
set.seed(22222)
Rt <- sim_simple(Tn=10, N=5, M=1e5+1)
rt_t <- sim_simple(Tn=10, N=5, M=1e4)
path <- function(M, Rt, rt_t){
  Rt <- lapply(Rt, function(rt, M) rt[-((M+1):(1e5+1)),], M=M)
  # qp_weights <- qp_weights(rt)
  # qp_result <- optim_dynam(Rt = rt, weights = qp_weights)

  # v[1] <- qp_result$v
  qp_l_weights <- lasso_weights(Rt, qp_lasso = TRUE)
  qp_l_theta <- optim_dynam_raw(Rt = Rt, weights = qp_l_weights)
  value_func(updatetheta =  qp_l_theta[1,], i = 1, theta = qp_l_theta, Rt = rt_t,  qp_l_weights,
             W = NULL, ValueOnly = F, M=1e4, Tn=10)$Value %>% mean()

  # v[2] <- qp_l_result$v
  # qp_cv_weights <- qp_weights(rt, constr = "cv")
  # qp_cv_result <- optim_dynam(Rt = rt, weights = qp_cv_weights)
  # v[3] <- qp_cv_result$v
  # qp_nc_weights <- qp_weights(rt, constr = NULL)
  # qp_nc_result <- optim_dynam(Rt = rt, weights = qp_nc_weights)
  # v[4] <- qp_nc_result$v
  # v
}
# value_path_e55 <- sapply(c(seq(1e4+5e3, 1e5, 5e3)), path, Rt = Rt, rt_t=rt_t)
value_path_u01 <- sapply(c(seq(100, 1e4, 100)), path, Rt = Rt, rt_t=rt_t)
# qplot(x=c(seq(100, 1e4, 100), seq(1e4+5e3, 1e5, 5e3)),
#       y = c(value_path, value_path_e55), geom = "line")
qplot(x=1/c(seq(100, 1e4, 100)),
      y = value_path_u01, geom = "line") + xlab("1/N")

write_csv(as.data.frame(value_path_u01), "path_raw_u01.csv")
write_csv(as.data.frame(value_path), "path_raw.csv")

# trail and error ---------------------------------------------------------


set.seed(22222)
rt <- sim_simple(Tn=10, N=5, M=1e4)
qp_l_weights <- lasso_weights(rt, qp_lasso = TRUE)
qp_l_result <- optim_dynam(Rt = rt, weights = qp_l_weights)

l_test <- optim_dynam_raw(Rt = rt, weights = qp_l_weights)



# plot --------------------------------------------------------------------

rt <- sim_simple(Tn=10, N=5, M=1e5)
path <- function(M, rt){
  rt <- lapply(rt, function(rt, M) rt[-((M+1):1e5),], M=M)
  qp_l_weights <- lasso_weights(rt, qp_lasso = TRUE)
  qp_l_result <- optim_dynam_raw(Rt = rt, weights = qp_l_weights)


}
value_path <- sapply(c(seq(100, 1e4, 100)), path, rt=rt)








Value <- numeric(M)
for(t in i:Tn){
  consu[t,] <- W[t,]*consu_function(theta[t,1]+theta[t,2]*W[t,])
  W[t+1,] <- (W[t,]-consu[t,])*(1+ t(Rt[[t]] %*% weights_list[[t]]) )
  # W[t+1,] <- (W[t,]-consu[t,])*exp(t(Rt[[t]] %*% weights_list[[t]]) )
  Value <- Value + utility(consu[t,])*beta^(t-1)
  # print(c(consu[t,],Value, W[t,], exp(t(Rt[[t]] %*% weights_list[[t]]) )))
  print(c(consu_function(theta[t,1]+theta[t,2]*W[t,]),consu[t,],Value, W[t,], (1+ t(Rt[[t]] %*% weights_list[[t]]) )))
}
