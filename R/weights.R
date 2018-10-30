# Weights function

#' no-short-sale portfolio quadratic programming
#' @author Yangzhuoran Yang
#' @export
qp_weights <- function(Rt){
  qp_weights <- mapply(qp_weights_do, Rt)
  return(qp_weights)
}


# no-short-sale portfolio quadratic programming

qp_weights_do <- function(Rt){
  Dmat <- cov(Rt)
  dvec <- rep(0,5)
  l1norm_A <- as.matrix(expand.grid(rep(list(c(-1,1)),5)))
  Amat <- rbind(rep(-1,5), l1norm_A)
  bvec <- rep(-1,NROW(l1norm_A)+1)
  qp <- quadprog::solve.QP(Dmat, dvec, t(Amat), bvec, meq=1)
  qp_weights <- qp$solution
  return(qp_weights)
}




#' Lasso portfolio selection
#'
#' Default using lasso after quadratic programming
#'
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


# Compute Rtw ----



# Rtw <- mapply(function(r,w) r %*% as.matrix(w), r=Rt, w=weights_list)


# Rtw <-t(weights) %*% t(Rt_bar)

