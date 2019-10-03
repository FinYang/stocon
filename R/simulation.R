# Simulation function



#' Return Simulation from multivariate normal distribution
#'
#' Simulate self-independent returns using Choleski Decomposition where correlation depends on the order of assets
#' from multvariate normal distribution
#'
#' The correlation function determines the correlation between asset by the difference between the index of the assets.
#' The default function is \code{rho = exp(-par * |distance|)}.
#' User defined correlation function need to have two arguments \code{i} and  \code{j} to indicate the position of assets
#'
#'
#' @param mu Either a scalar or a vector with length N contains mean of the assets returns
#' @param vol Either a scalar or a vector with length N contains volatilities of the assets returns
#' @param Tn Number of periods, excluding time 0.
#' @param N Number of assets
#' @param M Number of realization
#' @param varcor The variance covariance matrix of the assets. Diagonal elements must equal to vol squared.
#' If supplied, arguments \code{par} and \code{rho_do} will be ignored. If NULL, will be calculated using the correlation function and vol.
#' @param par The parameter in the default correlation function. See Details.
#' @param rho_do The function of correlation between assets. See below for the default function.
#' @param dependent If the following simulated series are added on the basis of existing ones.
#' @return List of returns
#' @author Yangzhuoran Yang
#' @importFrom magrittr %>%
#' @export
#'
sim_simple <- function(dis_par = list(mu = 0.05, vol = 0.02, df = 2), distribution = c("norm", "t"), Tn = 10, N = 5, M = 1e4, varcov = NULL,
                       par = 0.2, rho_do = NULL, dependent = c(NA, "AR1"),
                       den_par = list(AR1 = list(beta_0 = 0.05, beta_1 = 0.5)),
                       return_varcov = FALSE, plus_one = TRUE){


  mu <- dis_par$mu
  vol <- dis_par$vol
  if(distribution[[1]] == "t"){
    df <- dis_par$df

  }
  if ((length(mu) == 1) && (length(vol) == 1)) {
    mu <- rep(mu, N)
    vol <- rep(vol, N)
  } else {
    if ((length(mu) != N) | (length(vol) != N)) {
      stop("Mean vector or Vol vector length does not match the number of assets")
    }
  }
  time <- seq(from = 0, to = Tn, by = 1) # time span

  # covariance matrix
  if (!is.null(varcov)) {
    if (!all(diag(varcov) == vol)) {
      stop("Diagonal elements of covariance matrix must equal to vol squared")
    }
  } else {
    if (is.null(rho_do)) {
      rho_do <- function(i, j, parr = par) {
        exp(-parr * abs(i - j))
      }
    }
    rho_m <- sapply(1:N, function(i) {
      sapply(1:N, function(j) rho_do(i = i, j = j))
    })
    # rho_m[lower.tri(rho_m, TRUE)] <- 0
    # rho_m <- rho_m + t(rho_m) + diag(rep(1, N))

    varcov <- vol %*% t(vol) * rho_m
  }
  A <- chol(varcov)
  # -
  if(is.na(dependent[[1]])){
    Rt <- NULL
    if (N == 1) {
      for (t in 1:(Tn + 1)) {
        Z <- matrix(rnorm(M * N), ncol = N)
        if(distribution[[1]] == "t"){
          W <- matrix(df/rchisq(M, df))

          Rt[[t]] <- replicate(M, mu) + Z %*% A * W
        } else {
          Rt[[t]] <- replicate(M, mu) + Z %*% A

        }
      }
    } else {
      for (t in 1:(Tn + 1)) {
        Z <- matrix(rnorm(M * N), ncol = N)
        Rt[[t]] <- t(replicate(M, mu)) + Z %*% A
      }
    }
  } else {
    if("AR1" %in% dependent[[1]]){
      beta_0 <- den_par$AR1$beta_0
      beta_1 <- den_par$AR1$beta_1
      rt <- matrix(nrow=M, ncol=N)
      Rt <- lapply(seq_len(Tn), function(X) rt)
      for(i in 1:N)
        Rt[[1]][,i] <- rnorm(M, mean = mu[i], sd = vol[i])
      for(t in 2:(Tn+1)){
        Z <- matrix(rnorm(M*N), ncol=N)
        if(distribution[[1]] == "t"){
          W <- matrix(df/rchisq(M, df))
          Rt[[t]] <- beta_0 + beta_1 * Rt[[t-1]]  +  Z %*% A * W
        } else {
          Rt[[t]] <- beta_0 + beta_1 * Rt[[t-1]]  +  Z %*% A
        }
      }
    }
    #   Rt[[1]] <- sapply(1:N, function(i) rnorm(M, mean = mu[[i]], sd = vol[[i]]))
    #   for (t in 2:(Tn + 1)) {
    #     Z <- matrix(rnorm(M * N), ncol = N)
    #     Rt[[t]] <- Rt[[t - 1]] + t(replicate(M, mu)) + Z %*% A
    #   }
  }

  if(plus_one){
    Rt <- lapply(Rt, `+`, 1)
  }
  if(return_varcov){
    return(list(Rt = Rt, varcov =varcov))
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
