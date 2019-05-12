# Simulation function



#' Simple simulation
#'
#' Simulate returns using Choleski Decomposition where correlation depends on the order of assets
#'
#'
#' @param R1 Returns at time 0. Should be a M*N matrix where M and N are the number of entry and number of assets.
#' @param mu Vector with length N contains mean of the assets returns
#' @param vol Vector with length N contains volatilities of the assets returns
#' @param Tn Number of periods, excluding time 0.
#' @param N number of assets
#' @param M number of realization
#' @return List of returns
#' @author Yangzhuoran Yang
#' @importFrom magrittr %>%
#' @export
#'
sim_simple <- function(R1=NULL, mu=NULL, vol=NULL, Tn=10, N=5, M=1e4, for_weights = TRUE){
  if(length(R1 & mu & vol)==0){
    N <- 5
    R1 <- mu <- vol <- NULL
  } else if(!(NCOL(R1)==N & length(mu)==N & length(vol)==N)){
    warning("Conflicts in User given information. Droped")
    R1 <- mu <- vol <- NULL
    N <- 5
  }

  t_incre <- 1
  time <- seq(from=0, to=Tn, by=t_incre) # time span

  # Rt <- vector("list", T)
  rt <- matrix(nrow=M, ncol=N)
  Rt <- lapply(seq_len(Tn+1), function(X) rt)
  # return at time 0
  if(length(R1 & mu & vol) == 0){
    mu <- c(0.001, 0.0003, 0.0012, 0.0004, 0.0015)
    vol <- c(0.01, 0.06, 0.03, 0.07, 0.005)
    if(for_weights){
      for(i in 1:N)
        Rt[[1]][,i] <- rnorm(M, mean = mu[i], sd = vol[i])

    } else {
      for(i in 1:N)
        Rt[[1]][,i] <- mu[i]
    }


  } else {
    Rt[[1]] <- R1
  }
  # covariance matrix
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

  for(t in 2:(Tn+1)){
    Z <- matrix(rnorm(M*N), ncol=N)
    Rt[[t]] <- Rt[[t-1]] + t(replicate(M,mu))*t_incre + sqrt(t_incre) * Z %*% A
  }
  return(Rt)
}


# simulation CIR ------------------------------------------------------
#
#
# a <- c(rep(0.09, 5))
# sigma <- c(0.024, 0.026, 0.027, 0.06, 0.08)
# b <- c(0.07, 0.08, 0.09, 0.01, 0.02)
#
# d <- 4*a*b/sigma^2
# c <- sigma^2*(1-exp(-a*t_incre))/(4*a)
#
# Loop over time and assets -----------------------------------------------
#
# for(n in 1:N) {
#   if(d[n]>1){
#     for(t in 2:length(time)){
#       # Z<-randn(M,1)
#       Z <- rnorm(M)
#       # X<-chi2rnd(d-1,M,1)
#       X <- rchisq(M, d-1)
#       lambda <- Rt[[t-1]][,n]*exp(-a[n]*t_incre)/c[n]
#       Rt[[t]][,n]<-c[n]*((Z+sqrt(lambda))^2+X)
#     }
#   }else{
#     for(t in 2:length(time)){
#       lambda <- Rt[[t-1]][,n]*exp(-a[n]*t_incre)/c[n]
#       # P<-poissrnd(lambda/2,M,1)
#       P <- rpois(M, lambda/2)
#       # X=chi2rnd(d+2*N)
#       X <- rchisq(M, d[n]+2*P)
#       Rt[[t]][,n]<-c*X
#     }
#   }
# }






