library(tidyverse)
library(glmnet)
library(boot)
library(quadprog)


# big loop
# for(bl in 1:1){
v <- numeric(3)
set.seed(2222)
rt <- sim_simple(Tn=10, N=5, M=1e4)

qp_weights <- qp_weights(rt)
qp_result <- optim_dynam(Rt = rt, weights = qp_weights)
v[1] <- qp_result$v

qp_l_weights <- lasso_weights(rt, qp_lasso = TRUE)
qp_l_result <- optim_dynam(Rt = rt, weights = qp_l_weights)
v[2] <- qp_l_result$v
v




library(LowRankQP)
lasso.lr <- function(X,y,t) {
  n <- dim(X)[1]
  p <- dim(X)[2]
  Vn <- X %*% cbind (cbind (diag(p), -diag(p)), 0)
  yX = apply(sweep (X, MARGIN=1, -y, '*'), 2, sum)
  Zn = c (2*yX, -2*yX, 0)
  bOls = lm.fit(X, y)$coefficients
  u = c(abs(bOls), abs(bOls), sum(abs(bOls)))
  A = matrix (c(rep (1, 2*p), 1), nrow=1)
  b = c(min(t, sum(abs(bOls))))
  soln = LowRankQP(sqrt(2)*t(Vn), Zn, A, b, u, method="LU",verbose =F)
  return(round(soln$alpha[1:p] - soln$alpha[(p+1):(2*p)], digits=5))
}











neg_value_func(qp_result$theta[1,], 1, qp_result$theta,
           rt, M = NROW(rt[[1]]), Tn = length(rt) - 1,
           qp_weights, W = NULL, ValueOnly = T, utility = power_u,
           beta = 1/1.05, consu_function = boot::inv.logit)

value_func(qp_l_result$theta[1,], 1, qp_l_result$theta,
           rt, M = NROW(rt[[1]]), Tn = length(rt) - 1,
           qp_l_weights, W = NULL, ValueOnly = T, utility = power_u,
           beta = 1/1.05, consu_function = boot::inv.logit)

value_func(1, 1, matrix(rep(1,20), nrow = 10),
           rt, M = NROW(rt[[1]]), Tn = length(rt) - 1,
           qp_weights, W = NULL, ValueOnly = T, utility = power_u,
           beta = 1/1.05, consu_function = boot::inv.logit)
# }









##quadprog example"
Dmat       <- matrix(0,3,3)
diag(Dmat) <- 1
dvec       <- c(0,5,0)
Amat       <- matrix(c(-4,-3,0,2,1,0,0,-2,1),3,3)
bvec       <- c(-8,2,0)
qp <- quadprog::solve.QP(Dmat,dvec,Amat,bvec=bvec)
qpXT <- solveQPXT(Dmat,dvec,Amat,bvec=bvec)
range(qp$solution - qpXT$solution)

N <- 10
set.seed(2)
cr <- matrix(runif(N * N, 0, .05), N, N)
diag(cr) <- 1
cr <- (cr + t(cr)) / 2
set.seed(3)
sigs <- runif(N, min = .02, max = .25)
set.seed(5)
dvec <- runif(N, -.1, .1)
Dmat <- sigs %o% sigs * cr
Amat <- cbind(diag(N), diag(N) * -1)
bvec <- c(rep(-1, N), rep(-1, N))
resBase <- solveQPXT(Dmat, dvec, Amat, bvec)
##absolute value constraint on decision variable:
res <- solveQPXT(Dmat, dvec, Amat, bvec,
                 AmatPosNeg = matrix(rep(-1, 2 * N)), bvecPosNeg = -1)
sum(abs(res$solution[1:N]))

## penalty of L1 norm
resL1Penalty <- solveQPXT(Dmat, dvec, Amat, bvec, dvecPosNeg = -.005 * rep(1, 2 * N))
sum(abs(resL1Penalty$solution[1:N]))





##number of decision variables
N <- 5
##M is equal length(|b|) or 0; L = length(|b - b0|) or 0
M <- L <- 0
  M <- N

K <- ncol(Amat)
K1 <- ncol(AmatPosNeg)
K2 <- ncol(AmatPosNegDelta)

##expand original constraints to problem size
Amat <- rbind(
  Amat,
  matrix(0, M + L, K)
)

##create slack constraints: decision vector is: [b, |b|, |b - b0|]
AmatSlack <- cbind(
  rbind(
    diag(x = 1, N, M),
    diag(x = 1, M, M),
    matrix(0, L, M)
  ),
  rbind(
    diag(x = -1, N, M),
    diag(x = 1, M, M),
    matrix(0, L, M)
  ),
  rbind(
    diag(x = 1, N, L),
    matrix(0, M, L),
    diag(x = 1, L, L)
  ),
  rbind(
    diag(x = -1, N, L),
    matrix(0, M, L),
    diag(x = 1, L, L)
  )
)
b0 <- NULL
bvecSlack <- c(rep(0, M * 2), c(b0, b0 * -1))

##expand abs value constraint matrices
AmatAbs <- cbind(
  rbind(
    AmatPosNeg,
    matrix(0, 2 * L, K1)
  ),
  rbind(
    matrix(0, 2 * M, K2),
    AmatPosNegDelta
  )
)

##Map the problem as stated in docs to [b, |b|, |b - b0|] to reduce dimensionality
MAP <- cbind(
  rbind(
    diag(x = .5, N, M),
    diag(x = .5, M, M),
    matrix(0, L, M)
  ),
  rbind(
    diag(x = -.5, N, M),
    diag(x = .5, M, M),
    matrix(0, L, M)
  ),
  rbind(
    diag(x = .5, N, L),
    matrix(0, M, L),
    diag(x = .5, L, L)
  ),
  rbind(
    diag(x = -.5, N, L),
    matrix(0, M, L),
    diag(x = .5, L, L)
  )
)

AmatAbs <- MAP %*% AmatAbs

AMAT <- cbind(
  Amat,
  AmatSlack,
  AmatAbs
)

BVEC <- c(
  bvec,
  bvecSlack,
  bvecPosNeg,
  NULL
)














