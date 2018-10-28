library(tidyverse)
library(glmnet)
library(boot)
library(quadprog)



set.seed(2222)
t_incre <- 1
time <- seq(from=0, to=10, by=t_incre) # time span
Tn <- length(time)-1
N <- 5
M<-1e4



# big loop
# for(bl in 1:1){

rt <- sim_simple()
weights <- lasso_weights(rt)

weights_list <- weights %>% data.frame() %>% as.list()
# define weights ----
#weights <- matrix(c(0.2, 0.3, 0.5, 0, 0))


# wealth ----
W <- matrix(nrow = Tn+1, ncol = M)
W[1,] <- 1000 # initial wealth currently set to be the same
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


theta  <- numeric(Tn*2) %>% matrix(ncol = 2)

# consumption has been set to logit
ValueFunction <- function(updatetheta, i, theta,
                          Rt, weights_list, M, Tn, W, ValueOnly = FALSE){
  theta[i,] <- updatetheta
  Value <- numeric(M)
  consu <- matrix(nrow = Tn, ncol = M)

  for(t in i:Tn){
    consu[t,] <- W[t,]*inv.logit(theta[t,1]+theta[t,2]*W[t,])
    W[t+1,] <- (W[t,]-consu[t,])*(1+ t(Rt[[t]] %*% weights_list[[t]]) )
    Value <- Value + u(consu[t,])*beta^(t-1)
  }
  Value <- sum(Value + u(W[NROW(W),])*beta^Tn)
  if(ValueOnly) return(Value)
  list(W=W,Value=Value) %>% return()
}

W_V <- ValueFunction(0, 1, theta, Rt, weights_list, M, Tn, W)
W <- W_V$W

for(t in Tn:1){
  # theta(t,:)=fminunc(@(updatetheta)ValueFunction(updatetheta,t, theta, R,weight, M, 10, W),theta(t,:));
  theta[t,] <- optim(theta[t,], ValueFunction, i=t, theta = theta,
                     Rt = Rt, weights_list = weights_list, M=M, Tn=Tn, W=W,
                     ValueOnly = TRUE)$par
  # update?
  W <- ValueFunction(theta[t,], i=t, theta = theta,
                     Rt = Rt, weights_list = weights_list, M=M, Tn=Tn, W=W, ValueOnly = F)$W
}



# old V----
# Vt <- matrix(nrow=M, ncol=Tn+1)
#
# Vt[,Tn+1] <- beta^Tn * u(W[,Tn+1])
# v_cmu <- Vt[,Tn+1]
# for(i in Tn:1){
#
#   V <- function(thet, v_cmu){
#     consum <- W[,i]*inv.logit(thet[1]+thet[2]*Rtw[,i])
#     (mean(beta^i*u(consum)) + mean(v_cmu)) %>% return()
#   }
#   theta[i,] <- optim(c(0,0), V, v_cmu=v_cmu)$par
#   for(l in i:Tn){
#     consu[,l] <- W[,l]*inv.logit(theta[l,1]+theta[l,2]*Rtw[,l])
#     W[,l+1] <- (1+Rtw[,l])*(W[,l]-consu[,l])
#     Vt[,l] <- beta^(l-1)*u(consu[,l])
#   }
#   Vt[,Tn+1] <- beta^Tn * u(W[,Tn+1])
#
#   v_cum <- rowSums(Vt[,i:(Tn+1)])
# }




# }
