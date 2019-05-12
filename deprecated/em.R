library(tidyverse)


Tn <- 5
N <- 5
M <- 1e4

Rf <- 1.01
discount <- 1/1.05
W <- matrix(nrow = M, ncol = Tn+1)
W[,1] <- 1000

value_varmean <- function(update_par = NULL, i = 1, para, Rt, Rf,
                          M = NROW(Rt[[1]]), Tn = length(Rt), W, discount=1/1.01,
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
                  ini_W = 1000, discount=1/1.01,
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

Rt <- sim_simple(Tn = Tn, N=N, M=M)
# Rr <- step2(Rt)
pm[[1]] <- dytim(Rt, Rf, discount = discount)
for(it in 2:100){
  Rt <- simulate_riskys()
  # Rr <- step2(Rt)
  pm[[it]] <- dytim(Rt, Rf, para = pm[[it-1]][[1]], discount = discount)
}
v <- sapply(pm, function(x) x[[2]])
saveRDS(pm, "em_estimation.rds")
qplot(y= v, x=seq_along(v), geom = "line")
plotly::ggplotly()
# + geom_line(aes(y=v, x=x, color = "red"), data = data.frame(v=v100, x=seq_along(v100)))
