library(tidyverse)
Tn <- 10
N <- 5
M <- 1e4


Rt <- sim_simple(Tn=10, N=5, M=1e4) %>%
  lapply(function(x) x+1) # the return for each individual asset
Rr <- sapply(Rt, rowMeans)[,-11]# %>% t()
Rf <- 1.01
W <- matrix(nrow = M, ncol = Tn+1)
W[,1] <- 1000


value_varmean <- function(update_par = NULL, i = 1, para, Rr, Rf,
                       M = NROW(Rr), Tn = NCOL(Rr), W, discount=1/1.05,
                       lambda = 1/2, returnW = FALSE){
  if(!is.null(update_par)) para[i, ] <- update_par
  beta <- para[ ,1]
  C <- para[ ,2]
  Value <- numeric(M)
  for(t in i:(Tn)){
    W[,t+1] <- (W[,t]-beta[t])*Rf + beta[t]*Rr[,t] - C[t]
    Value <- Value + discount^(t-i+1)*(C[t]^2-2*lambda*C[t])
  }
  Value <- Value + discount^(Tn-i+1)*(W[,Tn+1]^2-2*lambda*W[,Tn+1])
  Value <- mean(Value)
  if(returnW) return(W)
  return(Value)
}





dytim <- function(Rr, Rf, valuefunction = value_varmean, M = NROW(Rr),
                  Tn = NCOL(Rr),
                  ini_W = 1000, discount=1/1.05,
                  lambda = 1/2){
  para <- matrix(0,nrow = Tn, ncol = 2)
  W <- matrix(nrow = M, ncol = Tn+1)
  W[,1] <- ini_W
  W <- valuefunction(para = para, Rr = Rr, Rf = Rf, W = W,
                     M = M, Tn = Tn, discount = discount, lambda = lambda,
                     returnW = TRUE)
  for(t in Tn:1){
    para[t,] <- optim(par = para[t,], fn = valuefunction, i=t, para = para,
                     Rr = Rr, Rf = Rf, W = W,
                     M = M, Tn = Tn, discount = discount, lambda = lambda)$par
  }
  return(para)
}

dytim(Rr, Rf)
