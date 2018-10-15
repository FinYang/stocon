library(tidyverse)
library(glmnet)
library(boot)
library(quadprog)



set.seed(2222)
time <- seq(from=0, to=3, by=1) # time span
t_incre <- 1
N <- 5
M<-1e4

theta  <- rep(0, (length(time)-1)*2) %>% matrix(ncol = 2)

# big loop
# for(bl in 1:1){
  
    # Rt <- vector("list", T) 
  rt <- matrix(nrow=M, ncol=N) 
  Rt <- lapply(seq_len(length(time)), function(X) rt)
  # return at time 0
  Rt[[1]][ ,1] <- runif(M, min=0.06, max=0.08)
  Rt[[1]][ ,2] <- runif(M, min=0.01, max=0.03)
  Rt[[1]][ ,3] <- runif(M, min=0.07, max=0.09)
  Rt[[1]][ ,4] <- runif(M, min=0.02, max=0.04)
  Rt[[1]][ ,5] <- runif(M, min=0.11, max=0.13)
  
  
  # simple simulation
  mu <- c(0.001, 0.0003, 0.0012, 0.0004, 0.0015)
  vol <- c(0.01, 0.06, 0.03, 0.07, 0.005)
  
  # covariance matrix
  rho_do <- function(i,j, par=0.9){
    exp(-par*(i-j)) %>% return()
  }
  rho_m <- sapply(1:N, function(i) {
    sapply(1:N, function(j) rho_do(i=i, j=j))
  })
  rho_m[lower.tri(rho_m, TRUE)] <- 0
  rho_m <- rho_m + t(rho_m) + diag(rep(1, N))
  
  varcov <- vol %*% t(vol) * rho_m
  A <- chol(varcov)
  # -
  
  for(t in 2:(length(time))){
    Z <- matrix(rnorm(M*N), ncol=N)
    Rt[[t]] <- Rt[[t-1]] + t(replicate(M,mu))*t_incre + sqrt(t_incre) * Z %*% A
  }
  
  
  # simulation CIR ------------------------------------------------------
  # 
  # 
  # a <- c(rep(0.09, 5))
  # sigma <- c(0.024, 0.026, 0.027, 0.06, 0.08)
  # b <- c(0.07, 0.08, 0.09, 0.01, 0.02)
  # 
  # d <- 4*a*b/sigma^2
  # c <- sigma^2*(1-exp(-a*t_incre))/(4*a)
  # 
  # Loop over time and assets -----------------------------------------------
  # 
  # for(n in 1:N) {
  #   if(d[n]>1){
  #     for(t in 2:length(time)){
  #       # Z<-randn(M,1)
  #       Z <- rnorm(M)
  #       # X<-chi2rnd(d-1,M,1)
  #       X <- rchisq(M, d-1) 
  #       lambda <- Rt[[t-1]][,n]*exp(-a[n]*t_incre)/c[n]
  #       Rt[[t]][,n]<-c[n]*((Z+sqrt(lambda))^2+X)
  #     }
  #   }else{
  #     for(t in 2:length(time)){
  #       lambda <- Rt[[t-1]][,n]*exp(-a[n]*t_incre)/c[n]
  #       # P<-poissrnd(lambda/2,M,1)
  #       P <- rpois(M, lambda/2)
  #       # X=chi2rnd(d+2*N)
  #       X <- rchisq(M, d[n]+2*P) 
  #       Rt[[t]][,n]<-c*X
  #     } 
  #   }
  # }
  
  # define weights ----
  #weights <- matrix(c(0.2, 0.3, 0.5, 0, 0))
  #####
  
  # no-short-sale portfolio quadratic programming
  
  qp_weights_do <- function(Rt){
    Dmat <- cov(Rt)
    dvec <- rep(0,5)
    l1norm_A <- as.matrix(expand.grid(rep(list(c(-1,1)),5)))
    Amat <- rbind(rep(-1,5), l1norm_A)
    bvec <- rep(-1,NROW(l1norm_A)+1)
    qp <- quadprog::solve.QP(Dmat, dvec, t(Amat), bvec, meq=1)
    qp_weights <- qp$solution
    return(qp_weights)
  }
  qp_weights <- mapply(qp_weights_do, Rt)
  Rt_noshortsale <- mapply(function(Rt,w) Rt %*% w, 
                           Rt = Rt, 
                           w = as.data.frame(qp_weights))
  
  # lasso_weights ----
  
  ns_Rt <- mapply(function(ns, Rt) cbind(ns, Rt), ns = as.data.frame(Rt_noshortsale), Rt = Rt, SIMPLIFY = F)
  lasso_data <- lapply(ns_Rt, function(Rt) cbind(Rt[,1], Rt[,2:NCOL(Rt)]-Rt[,1]))
  lasso.weights.m <- mapply(function(data) 
    cv.glmnet(x=as.matrix(data[,2:NCOL(data)]), 
              y=as.matrix(data[,1]), 
              alpha = 1, 
              intercept = TRUE,
              lambda = exp(seq(log(0.00001), log(3), length.out=200))) %>% 
      list(), 
    data = lasso_data)
  lasso.weights <- mapply(coef, lasso.weights.m)
  lasso.weights <- do.call(cbind, lasso.weights) %>% as.matrix() %>% .[-1,]
  ns.weights <- apply(lasso.weights, 2, function(x) 1-sum(x))
  weights <- matrix(rep(ns.weights,N), byrow = T, nrow = N)*qp_weights + lasso.weights
  # lasso.weights <- rbind(1-colSums(lasso.weights), lasso.weights)
  colnames(lasso.weights) <- time 
  rownames(lasso.weights) <- paste("Asset",1:N, sep = "_")
  
  
  
  # Compute Rtw ----
  weights_list <- lasso.weights %>% data.frame() %>% as.list()
  Rtw <- mapply(function(r,w) r %*% as.matrix(w), r=Rt, w=weights_list)
  
  
  # Rtw <-t(weights) %*% t(Rt_bar) 
  # wealth ----
  W <- matrix(nrow = M, ncol = length(time))
  W[,1] <- 1000 # initial wealth currently set to be the same 
  # discount factor ----
  beta <- 1/1.05 #discount factor currently set to be the same 
  
  # utility_function ----
  # u <- function(x, alpha = 0.5){
  #   -exp(-alpha*x) %>% return()
  # }
  # 
  #power utility
  u <- function(x, lambda = 0.9){
    x^(1-lambda)/(1-lambda) %>% return()
  }
  consu <- matrix(nrow = M, ncol = length(time)-1)
  # Value
  
  
  for(i in 1:(length(time)-1)){
    consu[,i] <- W[,i]*inv.logit(theta[i,1]+theta[i,2]*Rtw[,i]) 
    W[,i+1] <- (1+Rtw[,i])*(W[,i]-consu[,i])
  }
  # consumption has been set to logit
  
  Vt <- matrix(nrow=M, ncol=length(time)) 
  
  Vt[,length(time)] <- beta^(length(time)-1) * u(W[,length(time)])
  v_cmu <- Vt[,length(time)] 
  for(i in (length(time)-1):1){
    
    V <- function(thet, v_cmu){
      consum <- W[,i]*inv.logit(thet[1]+thet[2]*Rtw[,i])
      (mean(beta^i*u(consum)) + mean(v_cmu)) %>% return()
    }
    theta[i,] <- optim(c(0,0), V, v_cmu=v_cmu)$par
    for(l in i:(length(time)-1)){
      consu[,l] <- W[,l]*inv.logit(theta[l,1]+theta[l,2]*Rtw[,l])
      W[,l+1] <- (1+Rtw[,l])*(W[,l]-consu[,l])
      Vt[,l] <- beta^(l-1)*u(consu[,l])
    }
    Vt[,length(time)] <- beta^(length(time)-1) * u(W[,length(time)])
    
    v_cum <- rowSums(Vt[,i:length(time)])
  }
  
  
  
  
# }
