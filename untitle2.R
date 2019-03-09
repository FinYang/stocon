library(tidyverse)
library(gridExtra)

set.seed(222222)
rt_w <- sim_simple(M=1e4+1)
rt_o <- sim_simple(M=1e4+1,for_weights = F)
rt_t <- sim_simple(for_weights = F)
path <- function(M, rt_w, rt_o, rt_t){
  rt_w <- lapply(rt_w, function(rt, M) rt[-((M+1):(1e4+1)),], M=M)
  rt_o <- lapply(rt_o, function(rt, M) rt[-((M+1):(1e4+1)),], M=M)
  # qp_weights <- qp_weights(rt)
  # qp_result <- optim_dynam(Rt = rt, weights = qp_weights)

  # v[1] <- qp_result$v

  # qp_l_weights <- lasso_weights(rt_w, qp_lasso = TRUE)
  qp_l_weights <- qp_weights(rt_w)
  weights_list <- qp_l_weights %>% data.frame() %>% as.list()
  sapply(1:10, function(t) exp(rt_o[[t]]) %*% weights_list[[t]]) %>% apply(2,mean)
}

  qp_l_theta <- optim_dynam_raw(Rt = rt_o, weights = qp_l_weights)
  V <- value_func(updatetheta =  qp_l_theta[1,], i = 1, theta = qp_l_theta, Rt = rt_t,  qp_l_weights,
             W = NULL, ValueOnly = F, M=1e4, Tn=10)
  mV <- V$Value %>% mean()
  W <- V$W
  return(list(qp_l_weights, qp_l_theta, mV, W))

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
value_path_3s10 <- lapply(c(100,seq(1000, 1e4, 1000)), path, rt_w = rt_w, rt_o = rt_o, rt_t=rt_t)
# qplot(x=c(seq(100, 1e4, 100), seq(1e4+5e3, 1e5, 5e3)),
#       y = c(value_path, value_path_e55), geom = "line")

value_path <- sapply(value_path_3s10, function(x) x[[3]])
p1 <- qplot(x=c(100,seq(1000, 1e4, 1000)),
            y = value_path, geom = "line") + xlab("N")
p2<- qplot(x=1/(c(100,seq(1000, 1e4, 1000))^(1)),
      y = value_path, geom = "line") + xlab("1/N")
grid.arrange(p1,p2)


weights_list <- qp_l_weights %>% data.frame() %>% as.list()
sapply(1:10, function(t) exp(rt_o[[t]]) %*% weights_list[[t]]) %>% View()
sapply(1:10, function(t) exp(rt_o[[t]]) %*% weights_list[[t]]) %>% apply(2, mean) ->x
qplot(x=c(100,seq(1000, 1e4, 1000)), y=x, geom="line")
#
erw <- lapply(c(100,seq(200, 1e4, 100)), path, rt_w = rt_w, rt_o = rt_o, rt_t=rt_t)
qplot(x=c(100,seq(200, 1e4, 100)), y=sapply(erw, function(x) x[7]), geom="line")

sapply(1:10, function(t) boot::inv.logit(qp_l_theta[t,1]+qp_l_theta[t,2]*V$W[t,])) %>% View


write_csv(as.data.frame(value_path_3s), "path_3s.csv")
