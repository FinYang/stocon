library(tidyverse)


Tn <- 10
N <- 5
M <- 1e4

step1 <- function(){
  rt <- matrix(nrow=M, ncol=N)
  Rt <- lapply(seq_len(Tn), function(X) rt)
  mu <- c(0.01, 0.003, 0.012, 0.004, 0.015)
  vol <- c(0.01, 0.06, 0.03, 0.07, 0.005)
  for(i in 1:N)
    Rt[[1]][,i] <- rnorm(M, mean = mu[i], sd = vol[i])
  rho_do <- function(i,j, par=0.2){
    exp(-par*(i-j))
  }
  rho_m <- sapply(1:N, function(i) {
    sapply(1:N, function(j) rho_do(i=i, j=j))
  })
  rho_m[lower.tri(rho_m, TRUE)] <- 0
  rho_m <- rho_m + t(rho_m) + diag(rep(1, N))

  varcov <- vol %*% t(vol) * rho_m
  A <- chol(varcov)
  # -

  for(t in 2:Tn){
    Z <- matrix(rnorm(M*N), ncol=N)
    Rt[[t]] <- Rt[[t-1]] + t(replicate(M,mu)) +  Z %*% A
  }

  Rt <- Rt %>%
    lapply(function(x) x+1) # simulated return for each resky assets
  return(Rt)
}



## EM
# step2 <- function(Rt){
# lm_weight <- function(Rt){
#   reg_data <- cbind(Rt[,1], Rt[,1] - Rt[,2:NCOL(Rt)])
#   model <- lm(reg_data[,1]~reg_data[,-1])
#   weights <- c(1-sum(coef(model)[-1]), coef(model)[-1])
#   names(weights) <- NULL
#   return(weights)
# }
# weights <- sapply(Rt, lm_weight)
#
#
# Rr <- mapply(function(Rt, omega) Rt %*% omega, Rt = Rt, omega = as.data.frame(weights))
# return(Rr)
# }

Rf <- 1.01
W <- matrix(nrow = M, ncol = Tn+1)
W[,1] <- 1000


value_varmean <- function(update_par = NULL, i = 1, para, Rt, Rf,
                          M = NROW(Rt[[1]]), Tn = length(Rt), W, discount=1/1.05,
                          lambda = 1/2, returnW = FALSE){
  if(!is.null(update_par)) para[i, ] <- update_par
  beta <- para[ ,1]
  C <- para[ ,2]
  omega <- para[,-1:-2]
  Value <- numeric(M)
  for(t in i:(Tn)){
    W[,t+1] <- (W[,t]-beta[t])*Rf + beta[t]*(Rt[[t]] %*% omega[t,]) - C[t]
    Value <- Value + discount^(t-i+1)*(C[t]^2-2*lambda*C[t])
  }
  Value <- Value + discount^(Tn-i+1)*(W[,Tn+1]^2-2*lambda*W[,Tn+1])
  Value <- mean(Value)
  if(returnW) return(W)
  return(Value)
}


dytim <- function(Rt, Rf, valuefunction = value_varmean, M = NROW(Rt[[1]]),
                  Tn = length(Rt),
                  ini_W = 1000, discount=1/1.05,
                  lambda = 1/2, para = NULL){
  if(is.null(para)){
    para <- matrix(0,nrow = Tn, ncol = 2+NCOL(Rt[[1]]))
    para[,-1:-2] <- 1/NCOL(Rt[[1]])
  }
  W <- matrix(nrow = M, ncol = Tn+1)
  W[,1] <- ini_W
  W <- valuefunction(para = para, Rt = Rt, Rf = Rf, W = W,
                     M = M, Tn = Tn, discount = discount, lambda = lambda,
                     returnW = TRUE)
  for(t in Tn:1){
    para[t,] <- optim(par = para[t,], fn = valuefunction, i=t, para = para,
                      Rt = Rt, Rf = Rf, W = W,
                      M = M, Tn = Tn, discount = discount, lambda = lambda)$par
  }
  colnames(para) <- c("beta", "C", paste0("omega", 1:NCOL(Rt[[1]])))
  value <- valuefunction(para = para, Rt = Rt, Rf = Rf, W = W,
                         M = M, Tn = Tn, discount = discount, lambda = lambda)
  return(list(para, value))
}

pm <- NULL

Rt <- step1()
# Rr <- step2(Rt)
pm[[1]] <- dytim(Rt, Rf)
for(it in 2:100){
  Rt <- step1()
  # Rr <- step2(Rt)
  pm[[it]] <- dytim(Rt, Rf, para = pm[[it-1]][[1]])
}
v <- sapply(pm, function(x) x[[2]])
qplot(y= v, x=seq_along(v), geom = "line")
