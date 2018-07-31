library(tidyverse)
library(glmnet)



set.seed(2222)
# new simulation CIR ------------------------------------------------------

N <- 5
a <- c(rep(0.09, 5))
sigma <- c(0.024, 0.026, 0.027, 0.003, 0.007)
b <- c(0.07, 0.08, 0.09, 0.01, 0.02)

time <- seq(from=0, to=5, by=1) # time span
t_incre <- 1

d <- 4*a*b/sigma^2
c <- sigma^2*(1-exp(-a*t_incre))/(4*a)


M<-1e4

# Rt <- vector("list", T) 
rt <- matrix(nrow=M, ncol=N) 
Rt <- lapply(seq_len(length(time)), function(X) rt)
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
weights <- matrix(c(0.2, 0.3, 0.5, 0, 0))
#####


# getting Rtw using simulated Rt and defined weights
Rtw <- lapply(Rt, function(x) x %*% weights) %>% lapply(as.vector)
# Using lasso to get estimated weights. 
#currently drop the intercept after estimation
lasso.weights.m <- mapply(function(y, x) cv.glmnet(x, y, alpha = 1) %>% list(), x=Rt,y=Rtw)
lasso.weights <- mapply(coef, lasso.weights.m)
lasso.weights <- do.call(cbind, lasso.weights) %>% as.matrix() %>% .[-1,]
colnames(lasso.weights) <- time #currently just drop the intercept after estimation
# Rtw <-t(weights) %*% t(Rt_bar) 
Rtw <- mapply(function(list, newx) predict(list, newx), lasso.weights.m, Rt)
Rtw <- colMeans(Rtw) 
# Rtw used in dynimic estimation are the average of the fitted Rtw from lasso models
W <- numeric(length(time))
W[1] <- 1000 # initial wealth
beta <- 1/1.05 #discount factor
# utility function, currently just identity
u <- function(W){ 
  return(W)
}

c <- numeric(length(time)-1)
# Value
V <- function(theta, decision = length(time)-1){
  for(i in 1:(length(time)-1)){
    c[i] <- W[i]*inv.logit(theta[1]+theta[2]*Rtw[i]) # consumption has been set to logit
    W[i+1] <- Rtw[i]*(W[i]-c[i])
  }
  Vt <- numeric(length(decision:(length(time)-1)))
  Vt[length(Vt)] <- beta^(length(time)-1) * u(W[length(W)])
  for(i in (length(decision:(length(time)-1))-1):1){
    t <- (decision:(length(time)-1))[i]
    Vt[i] <- beta^t*u(c[t+1])+Vt[i+1]
  }
  # for(i in 1:(length(decision:(length(time)-1))-1)){
  #   t <- (decision:(length(time)-1))[i]
  #   Vt[i] <- beta^t*u(c[t+1])
  # }
  sum(Vt) %>% return()
}

optim(c(10, 3), V)

