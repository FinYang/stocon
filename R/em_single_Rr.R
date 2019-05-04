library(tidyverse)


Tn <- 5
N <- 5
M <- 1e4

Rf <- 1.01
discount <- 1/1.05
W <- matrix(nrow = M, ncol = Tn+1)
W[,1] <- 1000

# para_beta <- numeric(3)
beta_function <- function(para_beta, w){
  para_beta[[1]] + para_beta[[2]]*w + para_beta[[3]] * w^2
}

# para_c <- numeric(2)
c_function <- function(para_c, w){
  w*exp(para_c[[1]]+para_c[[2]]*w)/(1+exp(para_c[[1]]+para_c[[2]]*w))
}

value_varmean <- function(update_par = NULL, i = 1, para, Rr, Rf,
                          M = NROW(Rr[[1]]), Tn = length(Rr), W, discount=1/1.01,
                          lambda = 1/2, returnW = FALSE){
  if(!is.null(update_par)) para[i, ] <- update_par
  para_beta <- para[ ,1:3]
  para_c <- para[ ,4:5]
  Value <- numeric(M)
  # t-1 : Tn-1
  for(t in i:(Tn)){
    W[,t+1] <- (W[,t]-beta_function(para_beta = para_beta[t, ], w = W[,t]))*Rf +
                  beta_function(para_beta = para_beta[t, ], w = W[,t])*Rr[[t]] - c_function(para_c = para_c[t, ], w = W[,t])
    Value <- Value +
      discount^(t-i+1)*(c_function(para_c = para_c[t, ], w = W[,t])^2-2*lambda*c_function(para_c = para_c[t, ], w = W[,t]))
  }
  Value <- Value + discount^(Tn-i+1)*(W[,Tn+1]^2-2*lambda*W[,Tn+1])
  Value <- mean(Value)
  if(returnW) return(W)
  return(Value)
}


dytim <- function(Rr, Rf, valuefunction = value_varmean, M = NROW(Rr[[1]]),
                  Tn = length(Rr),
                  ini_W = 1000, discount=1/1.01,
                  lambda = 1/2, para = NULL){
  if(is.null(para)){
    para <- matrix(0,nrow = Tn, ncol = 5)
    # para[,-1:-2] <- 1/NCOL(Rr[[1]])
  }
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
  colnames(para) <- c(paste0("beta_", 1:3), paste0("c_", 1:2))
  value <- valuefunction(para = para, Rr = Rr, Rf = Rf, W = W,
                         M = M, Tn = Tn, discount = discount, lambda = lambda)
  return(list(para, value))
}

pm <- NULL

# Rr <- simulate_riskys()
# Rr <- step2(Rr)
mean_Rr <- 1.05
sd_Rr <- 0.03
Rr <- replicate(Tn,  matrix(rnorm(M, mean_Rr, sd_Rr), ncol = 1), simplify = F)

pm[[1]] <- dytim(Rr, Rf, discount = discount)
pb <- txtProgressBar(min = 2, max = 100, style = 3)
for(it in 2:100){
  Rr <- replicate(Tn,  matrix(rnorm(M, mean_Rr, sd_Rr), ncol = 1), simplify = F)
  # Rr <- step2(Rr)
  pm[[it]] <-
    try(dytim(Rr, Rf, para = pm[[it-1]][[1]], discount = discount))
  setTxtProgressBar(pb, it)
}
v <- sapply(pm, function(x) x[[2]])
# saveRDS(pm, "em_estimation.rds")
qplot(y= v, x=seq_along(v), geom = "line")
plotly::ggplotly()
# + geom_line(aes(y=v, x=x, color = "red"), data = data.frame(v=v100, x=seq_along(v100)))


value_varmean(para = pm[[100]][[1]], Rr = Rr, Rf = Rf, W = W,
  M = M, Tn = Tn, discount = discount, lambda = lambda)
