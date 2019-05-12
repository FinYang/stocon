
library(tidyverse)
library(furrr)


our_method <- function(mean_Rr, sd_Rr){
  discount = 1/1.05
  lambda = 1/2
  Tn <- 5
  N <- 5
  M <- 1e4
  Rf <- 1.01
  W <- matrix(nrow = M, ncol = Tn+1)
  W[,1] <- 1000
  mean_J1 <- mean_Rr - Rf
  mean_J <- replicate(Tn,list(matrix(c(mean_J1, -1), 2)))

  mean_JJ <- matrix(c((mean_J1^2+sd_Rr^2), -(mean_Rr-Rf), -(mean_Rr-Rf), 1), 2)
  mean_JJ <- replicate(Tn,list(mean_JJ))

  # HDGF
  # H  0 : Tn-1 = Tn
  # D  0 : Tn   = Tn+1
  # G  0 : Tn   = Tn+1
  # F  0 : Tn   = Tn+1
  Ht <- vector("list", Tn)
  Dt <- vector("list", Tn+1)
  Gt <- vector("list", Tn+1)
  Ft <- vector("list", Tn+1)
  Dt[[Tn+1]] <- 1
  Gt[[Tn+1]] <- 0
  Ft[[Tn+1]] <- 0
  for(i in Tn:1){
    Ht[[i]] <- matrix(c(0,0,0,1), 2) + Dt[[i+1]]*mean_JJ[[i]]
    mul <- (1-Dt[[i+1]]*t(mean_J[[i]]) %*% solve(Ht[[i]]) %*% mean_J[[i]])
    Dt[[i]] <- c(discount*Rf^2*Dt[[i+1]] * mul)
    Gt[[i]] <- c((discount*Rf*Gt[[i+1]] + discount*Rf*(Rf -2)*lambda*Dt[[i+1]]) * mul )
    Ft[[i]] <- c((discount*(Rf-2)^2*lambda^2*Dt[[i+1]] + 2*discount*(Rf-2)*lambda*Gt[[i+1]] + discount*Ft[[i+1]]) * mul)
  }

  Wt <- W[1,] - lambda
  Wt <- split(Wt, rep(1:length(Wt), each = 1))
  # Solution of value function
  get_Vt <- function(Wt, Dt, Gt, Ft){
    Wt^2*Dt + 2*Wt*Gt + Ft
  }
  # Solution of control
  get_Zt <- function(Wt, Ht, Dt1, Gt1, Rf, lambda, mean_J){
    lapply(Wt,function(Wt) -(Dt1 * (Wt*Rf + (Rf-2)*lambda) +Gt1) * solve(Ht) %*% mean_J)
    # -solve(Ht) %*% (Dt1 * (Wt*Rf + (Rf-2)*lambda)) %*% mean_J
  }
  # Evolution of Wt
  evo_Wt <- function(Wt, Rf, lambda, Jt, Zt){
    # mapply(function(Wt, Jt, Zt) c(Wt*Rf +(Rf-2)*lambda + t(Jt) %*% Zt),Wt = Wt, Jt = Jt, Zt = Zt, SIMPLIFY = FALSE)
    c(Wt*Rf +(Rf-2)*lambda + t(Jt) %*% Zt)
  }

  V0 <- get_Vt(Wt = Wt[[1]], Dt = Dt[[1]], Gt = Gt[[1]], Ft = Ft[[1]])
  adi <- discount^(1:Tn)
  cumu_adi <- rev(cumsum(adi))
  adii <- cumu_adi +rev(adi)
  # ft <- sapply(Vt, mean)-lambda^2*adii
  V0 <- V0-adii[[1]]
  V0
}




em <- function(mean_Rr, sd_Rr){
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
  Rr <- replicate(Tn,  matrix(rnorm(M, mean_Rr, sd_Rr), ncol = 1), simplify = F)

  pb <- txtProgressBar(min = 0, max = 20, style = 3)
  pm[[1]] <- dytim(Rr, Rf, discount = discount)
  setTxtProgressBar(pb, 1)
  for(it in 2:20){
    Rr <- replicate(Tn,  matrix(rnorm(M, mean_Rr, sd_Rr), ncol = 1), simplify = F)
    # Rr <- step2(Rr)
    pm[[it]] <- dytim(Rr, Rf, para = pm[[it-1]][[1]], discount = discount)
    setTxtProgressBar(pb, it)
  }
  v <- sapply(pm, function(x) x[[2]])
  return(v)
}













a <- (3:7)/100 +1
b <- (1:5)/100
ablist <- expand.grid(a,b) %>% as.matrix() %>% split(rep(1:25, 2))
mean_Rr <- 1.05
sd_Rr <- 0.03
# parallel::detectCores()
# future::plan(future::multiprocess)
# future::plan(future::sequential)
# on.exit(future::plan(old_plan))
# profvis::profvis(em(mean_Rr = mean_Rr, sd_Rr = sd_Rr))
loop <- function(ab){
  mean_Rr <- ab[[1]]
  sd_Rr <- ab[[2]]
  V <- em(mean_Rr = mean_Rr, sd_Rr = sd_Rr)
  om <- our_method(mean_Rr = mean_Rr, sd_Rr = sd_Rr)
  list(V,om)
}

# pbapply <- future_map(a, function(x) future_map(b, loop, sd_Rr = x, .progress = T), .progress = T)

# end <- future_map(ablist, loop, .progress = T)
end <- NULL
for(i in 10:length(ablist)){
  e <- try(loop(ablist[[i]]))
  while(class(e) == "try-error"){
    e <- try(loop(ablist[[i]]))
  }
  end[[i]] <- e

  print(i)
}
saveRDS(end, "rainbow.rds")
#
# loop(ablist[[6]])
#
# our_method(1.03, 0.02)
# em(1.03, 0.02)
# loop(c(1.03, 0.02))
loop(ablist[[1]])
