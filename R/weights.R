# Weights function

#' portfolio quadratic programming for multiperiod data set
#' @param rt List of return
#' @param constr If numeric, take as the imposed constraint.
#' If NULL, nu constraint. Otherwise compute proper constraint using 10-fold cross-validation
#' @return Weights for each periods
#' @author Yangzhuoran Yang
#' @export
qp_weights <- function(rt, constr = 1){
  if(!is.list(rt)) rt <- list(rt)
  if(length(constr) == 0){
    qp_weights <- mapply(qp_weights_noc, rt)
  } else if(is.numeric(constr)){
    qp_weights <- mapply(qp_weights_do, rt, MoreArgs = list(constr = constr))
  } else {
    qp_weights <- mapply(cv.qp_weights_do, rt)
  }
  return(qp_weights)
}

# portfolio quadratic programming with no constraint
qp_weights_noc <- function(rt){
  n <- NCOL(rt)
  Dmat <- cov(rt)
  dvec <- matrix(numeric(n))
  Amat <- matrix(rep(1,n))
  bvec <- matrix(1)
  qp <- quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq=1)
  qp$solution
}

# protfolio quadratic programming with constraint selected using cross validation
cv.qp_weights_do <- function(rt){
  n <- NCOL(rt)
  Dmat <- cov(rt)
  dvec <- matrix(numeric(n))
  Amat <- matrix(rep(1,n))
  bvec <- matrix(1)
  qp <- quadprog::solve.QP(Dmat, dvec, Amat, bvec, meq=1)
  c_max <- sum(abs(qp$solution))
  if(c_max<1) c_max <- 1
  c_grid <- seq(1, c_max, 0.01)
  constr <- c_grid[which.min(sapply(c_grid, cv_risk, rt=rt))]
  qp_weights_do(rt, constr = constr)
}

# calculate cv risk given different constr
cv_risk <- function(constr, rt){
  folds <- caret::createFolds(seq(NROW(rt)))
  ev <- function(i, rt, folds, constr){
    wei <- suppressWarnings(try(qp_weights_do(rt[-folds[[i]],], constr = constr)))
    as.numeric(t(wei) %*% cov(rt[folds[[i]],]) %*% wei)
  }
  mean(sapply(1:10, ev, rt = rt, folds = folds, constr = constr))
}



#'  portfolio quadratic programming
#' @importFrom quadprogXT solveQPXT
#' @importFrom quadprog solve.QP
qp_weights_do <- function(rt, constr = 1, method = c("quadprogXT")){
  n <- NCOL(rt)

  if(method %in% "quadprogXT"){

    Dmat <- cov(rt)
    dvec <- numeric(n)
    Amat <- matrix(rep(1,n))
    AmatPosNeg <- matrix(rep(-1, 2 * n))
    bvecPosNeg <- -constr
    bvec <- 1
    qp <- try(quadprogXT::solveQPXT(Dmat, dvec, Amat, bvec, meq=1,
                                    AmatPosNeg = AmatPosNeg, bvecPosNeg = bvecPosNeg), silent = T)
    if("try-error" %in% class(qp)) {
      warning("Possible overflow error. Attempted to scale Dmat.")

      sc <- norm(Dmat,"2")
      qp <- try(quadprogXT::solveQPXT(Dmat/sc, dvec/sc, Amat, bvec, meq=1,
                                      AmatPosNeg = AmatPosNeg, bvecPosNeg = bvecPosNeg), silent = T)
      if("try-error" %in% class(qp)){
        warning("Attempt failed. Switch to `qp_orig`")
        return(qp_orig(rt, constr))
      }
    }

    qp_weights <- qp$solution
    return(qp_weights[1:n])
  } else {
    # c <- numeric(n*2)
    # Q <- cov(rt)
    # Q_nQ <- cbind(Q, -Q)
    # H <- rbind(Q_nQ, -Q_nQ)
    # A <- matrix( c(rep(1, n), rep(-1,n), rep(1, 2*n)), byrow = T, nrow = 2)
    # b <- matrix(c(1,0))
    # r <- matrix(c(0,constr))
    # l <- matrix(numeric(2*n))
    # u <- matrix(rep(1, 2*n))
    # kernlab::ipop(t(c), H, A, b, l, u, r)
    #
    # quadprog::solve.QP(Dmat = H, dvec = c,
    #          Amat = cbind(matrix(c(rep(-1, n), rep(1, n), rep(1, n), rep(-1, n), rep(-1, 2*n)), nrow = 2*n), diag(2*n)),
    #          bvec = c(-1,1, -constr, numeric(2*n)))
  }

}

# original qp with the number of parametar of 2^n
# not recomended in high dimensions
qp_orig <- function(rt, constr = 1){
  n <- NCOL(rt)
  Dmat <- cov(rt)
  dvec <- numeric(n)
  l1norm_A <- as.matrix(expand.grid(rep(list(c(-1,1)),n)))
  Amat <- rbind(rep(1,n), l1norm_A)
  bvec <- c(1, rep(-constr,NROW(l1norm_A)))
  qp <- quadprog::solve.QP(Dmat, dvec, t(Amat), bvec, meq=1)
  qp$solution
}




#' Lasso portfolio selection
#'
#' Default using lasso after quadratic programming. As a approximation to constrained risk
#' minimization
#'
#' @param Rt Lists of return rt
#' @param qp_lasso if TRUE, taking no-short-sale portfolio from quadratic programming as y in lasso
#' @param qp_weights can supply weights of o-short-sale portfolio from quadratic programming
#'
#' @return List of weights for each period
#' @author Yangzhuoran Yang
#' @importFrom magrittr %>%
#' @export
lasso_weights <- function(Rt, N = NCOL(Rt[[1]]), qp_lasso = TRUE, qp_weights = NULL){
  if(qp_lasso){
    if(length(qp_weights) == 0) qp_weights <- qp_weights(Rt)
    Rt_noshortsale <- mapply(function(Rt,w) Rt %*% w,
                                         Rt = Rt,
                                         w = as.data.frame(qp_weights))
  input_Rt <- mapply(function(ns, Rt) cbind(ns, Rt), ns = as.data.frame(Rt_noshortsale), Rt = Rt, SIMPLIFY = F)
  } else input_Rt <- Rt
  lasso_data <- lapply(input_Rt, function(Rt) cbind(Rt[,1], Rt[,2:NCOL(Rt)]-Rt[,1]))


  # lasso_weights ----
  lasso.weights.m <- mapply(function(data)
    glmnet::cv.glmnet(x=as.matrix(data[,2:NCOL(data)]),
                      y=as.matrix(data[,1]),
                      alpha = 1,
                      intercept = TRUE,
                      lambda = exp(seq(log(0.00001), log(3), length.out=200))) %>%
      list(),
    data = lasso_data)
  lasso.weights <- mapply(coef, lasso.weights.m)
  lasso.weights <- do.call(cbind, lasso.weights) %>% as.matrix() %>% .[-1,]

  if(qp_lasso){
    ns.weights <- apply(lasso.weights, 2, function(x) 1-sum(x))
    weights <- matrix(rep(ns.weights,N), byrow = T, nrow = N)*qp_weights + lasso.weights
  } else weights <- lasso.weights
  # lasso.weights <- rbind(1-colSums(lasso.weights), lasso.weights)
  colnames(weights) <- 0:(length(Rt)-1)
  rownames(weights) <- paste("Asset",1:N, sep = "_")

  return(weights)
}



