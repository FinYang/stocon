library(tidyverse)
library(glmnet)
library(boot)



set.seed(2222)
time <- seq(from=0, to=5, by=1) # time span
t_incre <- 1
N <- 5
M<-1e4

theta  <- matrix(rep(0, (length(time)-1)*2), ncol = 2)

# big loop
for(bl in 1:3){
  
  # new simulation CIR ------------------------------------------------------
  

  a <- c(rep(0.09, 5))
  sigma <- c(0.024, 0.026, 0.027, 0.003, 0.007)
  b <- c(0.07, 0.08, 0.09, 0.01, 0.02)
  

  
  d <- 4*a*b/sigma^2
  c <- sigma^2*(1-exp(-a*t_incre))/(4*a)
  
  

  
  # Rt <- vector("list", T) 
  rt <- matrix(nrow=M, ncol=N) 
  Rt <- lapply(seq_len(length(time)), function(X) rt)
  # return at time 0
  Rt[[1]][ ,1] <- runif(M, min=0.06, max=0.08)
  Rt[[1]][ ,2] <- runif(M, min=0.07, max=0.09)
  Rt[[1]][ ,3] <- runif(M, min=0.08, max=0.1)
  Rt[[1]][ ,4] <- runif(M, min=0.05, max=0.15)
  Rt[[1]][ ,5] <- runif(M, min=0.01, max=0.03)
  
  # Loop over time and assets -----------------------------------------------
  
  for(n in 1:N) {
    if(d[n]>1){
      for(t in 2:length(time)){
        # Z<-randn(M,1)
        Z <- rnorm(M)
        # X<-chi2rnd(d-1,M,1)
        X <- rchisq(M, d-1) 
        lambda <- Rt[[t-1]][,n]*exp(-a[n]*t_incre)/c[n]
        Rt[[t]][,n]<-c[n]*((Z+sqrt(lambda))^2+X)
      }
    }else{
      for(t in 2:length(time)){
        lambda <- Rt[[t-1]][,n]*exp(-a[n]*t_incre)/c[n]
        # P<-poissrnd(lambda/2,M,1)
        P <- rpois(M, lambda/2)
        # X=chi2rnd(d+2*N)
        X <- rchisq(M, d[n]+2*P) 
        Rt[[t]][,n]<-c*X
      } 
    }
  }
  
  # define weights ----
  #weights <- matrix(c(0.2, 0.3, 0.5, 0, 0))
  #####
  
  # lasso_weights ----
  
  lasso_data <- lapply(Rt, function(Rt) cbind(Rt[,1], Rt[,2:N]-Rt[,1]))
  lasso.weights.m <- mapply(function(data) 
    cv.glmnet(x=as.matrix(data[,2:N]), as.matrix(data[,1]), 
              alpha = 1) %>% list(), data = lasso_data)
  lasso.weights <- mapply(coef, lasso.weights.m)
  lasso.weights <- do.call(cbind, lasso.weights) %>% as.matrix() %>% .[-1,]
  lasso.weights <- rbind(1-colSums(lasso.weights), lasso.weights)
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
  # currently just identity
  u <- function(x, alpha = 0.5){
    -exp(-alpha*x) %>% return()
  }
  
  #power utility
  # u <- function(x, lambda = 0.9){
  #   x^(1-lambda)/(1-lambda) %>% return()
  # }
  consu <- matrix(nrow = M, ncol = length(time)-1)
  # Value
  
  
  for(i in 1:(length(time)-1)){
    consu[,i] <- W[,i]*inv.logit(theta[i,1]+theta[i,2]*Rtw[,i]) # consumption has been set to logit
    W[,i+1] <- (1+Rtw[,i])*(W[,i]-consu[,i])
  }
  
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
  
  
  
  
}
