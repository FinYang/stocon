library(tidyverse)

discount = 1/1.05
lambda = 1/2

Tn <- 5
N <- 5
M <- 1e4


Rf <- 1.01
W <- matrix(nrow = M, ncol = Tn+1)
W[,1] <- 1000


# single Rr
mean_Rr <- 1.05
sd_Rr <- 0.017


mean_J1 <- mean_Rr - Rf
mean_J <- replicate(Tn,list(matrix(c(mean_J1, -1), 2)))

mean_JJ <- matrix(c((mean_J1^2+sd_Rr^2), -(mean_Rr-Rf), -(mean_Rr-Rf), 1), 2)
mean_JJ <- replicate(Tn,list(mean_JJ))

HDGF <- get_HDGF(Tn=Tn, lambda=lambda, mean_JJ = mean_JJ, mean_J = mean_J, Rf = Rf)
Ht <- HDGF[[1]]
Dt <- HDGF[[2]]
Gt <- HDGF[[3]]
Ft <- HDGF[[4]]

# W tilde
Wt <- W[1,] - lambda
Wt <- split(Wt, rep(1:length(Wt), each = length(Wt)))
Wt <- Wt[[1]]
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
# 0 : Tn-1
Vt <- vector("list", Tn)
Zt <- vector("list", Tn)
for(i in 1:Tn){
  Vt[[i]] <- get_Vt(Wt = Wt[[i]], Dt = Dt[[i]], Gt = Gt[[i]], Ft = Ft[[i]])
  Zt[[i]] <- get_Zt(Wt = Wt[[i]], Ht = Ht[[i]], Dt1 = Dt[[i+1]], Gt1 = Gt[[i+1]], Rf = Rf, lambda = lambda, mean_J = mean_J[[i]])
  # Wt[[i+1]] <- do.call(c,evo_Wt(Wt = Wt[[i]], Rf = Rf,lambda = lambda, Jt =Jt[[i]], Zt = Zt[[i]]))

  Wt[[i+1]] <- evo_Wt(Wt = Wt[[i]], Rf = Rf,lambda = lambda, Jt =mean_J[[i]], Zt = Zt[[i]])
}
V <- mapply(get_Vt, Wt = as.data.frame(Wt), Dt = Dt, Gt = Gt, Ft = Ft, SIMPLIFY = FALSE)

#
adi <- discount^(1:Tn)
cumu_adi <- rev(cumsum(adi))
adii <- cumu_adi +rev(adi)
# ft <- sapply(Vt, mean)-lambda^2*adii
V0 <- V0-adii[[1]]
V0





