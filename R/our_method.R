library(tidyverse)

discount = 1/1.01
lambda = 1/2

Tn <- 5
N <- 5
M <- 1e4

simulate_riskys <- function(){
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
    Rt[[t]] <- Rt[[t-1]]  +  Z %*% A #+t(replicate(M,mu))
  }

  Rt <- Rt %>%
    lapply(function(x) x+1) # simulated return for each resky assets
  return(Rt)
}

get_Rr <- function(Rt){
  lm_weight <- function(Rt){
    reg_data <- cbind(Rt[,1], Rt[,1] - Rt[,2:NCOL(Rt)])
    model <- lm(reg_data[,1]~reg_data[,-1])
    weights <- c(1-sum(coef(model)[-1]), coef(model)[-1])
    names(weights) <- NULL
    return(weights)
  }
  weights <- sapply(Rt, lm_weight)


  Rr <- mapply(function(Rt, omega) Rt %*% omega, Rt = Rt, omega = as.data.frame(weights))
  return(Rr)
}

Rt <- simulate_riskys()
Rr <- get_Rr(Rt)
Rf <- 1.02
W <- matrix(nrow = M, ncol = Tn+1)
W[,1] <- 1000

# First element in J
J1 <- (Rr-Rf) %>% split(rep(1:ncol(.), each = nrow(.)))
Jt <- lapply(J1, function(J1) lapply(J1, function(J1) matrix(c(J1, -1), 2)))
# J1t <- split(J1, rep(1:ncol(J1), each = nrow(J1)))
# test_meanJ <- lapply(Jt, function(Jt) matrix(rowMeans(do.call(cbind, Jt)), 2))
mean_J <- matrix(colMeans(Rr-Rf), nrow = 2, ncol = Tn, byrow = TRUE)
mean_J[2,] <- -1
mean_J <-  split(mean_J, rep(1:ncol(mean_J), each = nrow(mean_J))) %>%
  lapply(matrix)

## JJ
JJ1 <- (Rr-Rf)^2
JJ2 <- -(Rr-Rf)

mean_JJ1 <- colMeans(JJ1)
mean_JJ2 <- colMeans(JJ2)
mean_JJ <- mapply(function(JJ1, JJ2) matrix(c(JJ1, JJ2, JJ2, 1), 2), JJ1 = mean_JJ1, JJ2  = mean_JJ2, SIMPLIFY = F)

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

Wt <- W - lambda
Wt <- split(Wt, rep(1:ncol(Wt), each = nrow(Wt)))
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
  mapply(function(Wt, Jt, Zt) c(Wt*Rf +(Rf-2)*lambda + t(Jt) %*% Zt),Wt = Wt, Jt = Jt, Zt = Zt, SIMPLIFY = FALSE)
}
# 0 : Tn-1
Vt <- vector("list", Tn)
Zt <- vector("list", Tn)
for(i in 1:Tn){
  Vt[[i]] <- get_Vt(Wt = Wt[[i]], Dt = Dt[[i]], Gt = Gt[[i]], Ft = Ft[[i]])
  Zt[[i]] <- get_Zt(Wt = Wt[[i]], Ht = Ht[[i]], Dt1 = Dt[[i+1]], Gt1 = Gt[[i+1]], Rf = Rf, lambda = lambda, mean_J = mean_J[[i]])
  Wt[[i+1]] <- do.call(c,evo_Wt(Wt = Wt[[i]], Rf = Rf,lambda = lambda, Jt =Jt[[i]], Zt = Zt[[i]]))
}
# V <- mapply(get_Vt, Wt = as.data.frame(Wt), Dt = Dt, Gt = Gt, Ft = Ft, SIMPLIFY = FALSE)


adi <- discount^(1:Tn)
cumu_adi <- rev(cumsum(adi))
adii <- cumu_adi +rev(adi)
ft <- sapply(Vt, mean)-lambda^2*adii








