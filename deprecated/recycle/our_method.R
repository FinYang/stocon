library(tidyverse)

discount = 1/1.05
lambda = 1/2


set.seed(2222)
# Simulate the return of risky assets
mu <- c(0.01, 0.003, 0.012, 0.004, 0.015)
vol <- c(0.01, 0.06, 0.03, 0.07, 0.005)
Rt <- sim_simple(mu, vol, Tn=4, N=5, independent = FALSE)
Rt <- Rt %>%
  lapply(function(x) x+1)


get_Rr <- function(Rt){
  weights <- sapply(Rt, weights_lm)
  Rr <- mapply(function(Rt, omega) Rt %*% omega, Rt = Rt, omega = as.data.frame(weights))
  return(Rr)
}

Rr <- get_Rr(Rt)


Rf <- 1.01
W <- matrix(nrow = M, ncol = Tn+1)
W[,1] <- 1000


# about J -----------------------------------------------------------------
# all using mean
# need to change to nonparametric methods ----------


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

