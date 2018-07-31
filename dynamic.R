library(tidyverse)
library(glmnet)



set.seed(2222)
# old simulating Rt -----------------------------------------------------------

time <- seq(from=0, to=0.5, by=0.1) # time span
# Setting parameters:
# n: number of assets
# M: number of simulation
# alpha, b, sigma, d, C (capital) are used in CIR simulation
N <- 5
M <- 1e4
alpha <- 0.5
b <- 0.4
sigma <- 0.3

#Rt <- vector("list", leMgth(time))
rt <- matrix(nrow=M, ncol=N) 
Rt <- lapply(seq_len(length(time)), function(X) rt)


# setting the initial return at time 0
Rt[[1]][ ,1] <- runif(M, min=0.2, max=0.8)
Rt[[1]][ ,2] <- runif(M, min=0.2, max=0.8)
Rt[[1]][ ,3] <- runif(M, min=0.2, max=0.8)
Rt[[1]][ ,4] <- runif(M, min=0.2, max=0.8)
Rt[[1]][ ,5] <- runif(M, min=0.2, max=0.8)

# simulatin of return using CIR model 
d <- 4*b*alpha/sigma^2
if(d>1){
  C <- (sigma^2*(1-exp(-alpha*(time[seq(time)+1]-time[seq(time)])))/(4*alpha)) %>% na.omit()
  for(m in seq(N)){
    for(i in 1:(length(time)-1)){
      lambda <- Rt[[i]][,m]*(exp(-alpha*(time[i+1]-time[i])))/C[i]
      Z <- rnorm(M)
      X <- rchisq(M, d-1)
      Rt[[i+1]][, m] <- C[i]*((Z+sqrt(lambda))^2+X)
    }
    
  }
}
#Rt_bar <- sapply(Rt, colMeans) %>% t


# new simulation CIR ------------------------------------------------------

N <- 5
a <- c(rep(0.25, 5))
sigma <- c(rep(0.25, 5))
b <- c(rep(100, 5))
x_0 <- c(rep(100, 5))
T <- 2
d <- 4*a*b/sigma^2
c <- sigma^2*(1-exp(-a*T))/(4*a)
lambda <- x_0*exp(-a*T)/c

M<-1e4

# Rt <- vector("list", T) 
rt <- matrix(nrow=M, ncol=N) 
Rt <- lapply(seq_len(T), function(X) rt)

for(n in 1:N) {
  if(d[n]>1){
    for(t in 1:T){
      # Z<-randn(M,1)
      Z <- rnorm(M)
      # X<-chi2rnd(d-1,M,1)
      X <- rchisq(M, d-1)
      Rt[[t]][,n]<-c*((Z+sqrt(lambda[n]))^2+X)
    }
  }else{
    for(t in 1:T){
      
      # P<-poissrnd(lambda/2,M,1)
      P <- rpois(M, lambda/2)
      # X=chi2rnd(d+2*N)
      X<-rchisq(1, d+2*P)
      Rt[[t]][,n]<-c*X
    } 
  }
}
#define weights
weights <- matrix(c(0.2, 0.3, 0.5, 0, 0))
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

