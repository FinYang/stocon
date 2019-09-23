

# para_beta <- numeric(3)
#' The beta function used in he EM algorithm
#'
#' It is bascally a third polynomial of the wealth.
#'
#' @param para_beta Vector of length 3 containing the parameters in beta function
#' @param w Scalar. The wealth in the beta function
#' @author Yangzhuoran Yang
#' @export
beta_function <- function(para_beta, w){
  para_beta[[1]] + para_beta[[2]]*w + para_beta[[3]] * w^2
}

# para_c <- numeric(2)
#' The comsumption function used in the EM algorithm
#'
#' It is the wealth times the logit function of a first order polynomial of wealth
#'
#' @param para_c Vector of length 2 containing the parameters in the comsumption function
#' @param w Scalar. The Wealth in the comsumption function.
#' @author Yangzhuoran Yang
#' @export
c_function <- function(para_c, w){
  w*exp(para_c[[1]]+para_c[[2]]*w)/(1+exp(para_c[[1]]+para_c[[2]]*w))
}



#' The Variance-Mean value function used in the EM algorithm
#'
#' @param update_par Matrix. Each row contains the parameters at a time point.
#' The number of column are the sum of the number of parameter in beta function,
#' the number of parameters in the consumption function, and the number of weights of the assets minus 1
#' (the sum of the weights is one so only need to estimate N-1 weights).
#' The positions in the vector should also be in that order.
#' @param i Integer. The position of the current time point
#' @param Rr List of matrix. The rertuns of the assets.
#' Each list is for one time point and each column in the matrix is one assets.
#' @param Rf Scalar. The risk free rate. e.g. 1.01
#' @param M Integer. The number of instances for one assets at a time
#' @param Tn Integer. The length of time horizon. Not including 0.
#' @param W Matrix. Each column is one time point and each row is one observation.
#' @param discount Scalar. The rate of discount. e.g. 1/1.01
#' @param lambda The scalar that make the calculation easiler.
#' @param returnW If TRUE, return the wealth. If False, return the value from the value function.
#'
#' @author Yangzhuoran Yang
#'
#' @export
value_varmean <- function(update_par = NULL, i = 1, para, Rr, Rf,
                          M = NROW(Rr[[1]]), Tn = length(Rr), W, discount=1/1.01,
                          lambda = 1/2, returnW = FALSE, detail = FALSE){
  if(!is.null(update_par)) para[i, ] <- update_par
  para_beta <- para[ ,1:3]
  para_c <- para[ ,4:5]
  if(NCOL(Rr[[1]]) > 1){
    para_w <- para[,6:(NCOL(para))]
    para_w <- cbind(para_w, apply(para_w, 1, function(x) 1- sum(x)))
    Rt <- Rr
    Rr <- mapply(function(r, w) r %*% as.matrix(w), r=Rt, w=split(para_w, 1:NROW(para_w)), SIMPLIFY = F)
  }
  Value <- numeric(M)
  # t-1 : Tn-1
  BETA <- list()
  C <- list()
  for(t in (i):(Tn)){
    BETA[[t-i+1]] <- beta_function(para_beta = para_beta[t, ], w = W[,t])
    C[[t-i+1]] <- c_function(para_c = para_c[t, ], w = W[,t])
    W[,t+1] <- (W[,t]-BETA[[t-i+1]])*Rf + BETA[[t-i+1]]*Rr[[t]] - C[[t-i+1]]
    Value <- Value + discount^(t-i+1)*(C[[t-i+1]]^2-2*lambda*C[[t-i+1]])
  }
  Value <- Value + discount^(Tn-i+1)*(W[,Tn+1]^2-2*lambda*W[,Tn+1])
  Value <- mean(Value)
  if(returnW) return(W)
  if(detail) return(list(W=W, value = Value))
  return(Value)
}

#' The iteration fucntion used in the EM algorithm
#'
#' This is only for one round of E step and M step
#'
#' @param Rr List of matrix. The rertuns of the assets.
#' Each list is for one time point and each column in the matrix is one assets.
#' @param Rf Scalar. The risk free rate. e.g. 1.01
#' @param valuefunction Function. The function that returns the value being optimised.
#' The first argument needs to be the parameter to optimise.
#' @param M Integer. The number of instances for one assets at a time
#' @param Tn Integer. The length of time horizon. Not including 0.
#' @param ini_W Scalar. The initial wealth.
#' @param discount Scalar. The rate of discount. e.g. 1/1.01
#' @param lambda The scalar that make the calculation easiler.
#' @param para Matrix. Initial value of all the paramters. Each column is a parameter and each row is a time point.
#' @author Yangzhuoran Yang
#' @seealso \code{EM}
#' @export
round_EM <- function(Rr, Rf, valuefunction = value_varmean, M = NROW(Rr[[1]]),
                  Tn = length(Rr)-1,
                  ini_W = 1000, discount=1/1.01,
                  lambda = 1/2, para = NULL){
  N <- NCOL(Rr[[1]])
  if(is.null(para)){
    para <- matrix(0,nrow = Tn+1, ncol = 5)
    para <- cbind(para,matrix(1/N , nrow = Tn+1, ncol =  N-1))
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
  colnames(para) <- c(paste0("beta_", 1:3), paste0("c_", 1:2), paste0("weight", 1:(N-1)))
  # value <- list()
  value <- numeric(Tn)
  for(i in 1:(Tn)){
    value[[i]] <- valuefunction(para = para, i=i, Rr = Rr, Rf = Rf, W = W,
                           M = M, Tn = Tn, discount = discount, lambda = lambda)
    # value[[i]]$t <- i-1
  }
  wealth <- valuefunction(para = para, i=i, Rr = Rr, Rf = Rf, W = W,
                          M = M, Tn = Tn, discount = discount, lambda = lambda, returnW = TRUE)
  ex_para <- para
  para <- para[-NROW(para),]
  para_beta <- para[,1:3]
  para_c <- para[,4:5]
  weights <- para[,6:NCOL(para)]
  weights <- cbind(weights, 1-rowSums(weights))
  colnames(weights)[[NCOL(weights)]] <- paste0("weight", NCOL(weights))

  BETA <- mapply(beta_function, para_beta = split(para_beta, seq_len(NROW(para_beta))), w = as.data.frame(wealth)[,-NCOL(wealth)])
  C <- mapply(c_function, para_c = split(para_c, seq_len(NROW(para_c))), w = as.data.frame(wealth)[,-NCOL(wealth)])



  return(list(para = ex_para, value = value, wealth = wealth, BETA = BETA, C=C))
}


#' Wrapper to iterate E and M steps.
#'
#'
#' @param Rr List of matrix. The rertuns of the assets.
#' Each list is for one time point and each column in the matrix is one assets.
#' @param Rf Scalar. The risk free rate. e.g. 1.01
#' @param valuefunction Function. The function that returns the value being optimised.
#' The first argument needs to be the parameter to optimise.
#' @param M Integer. The number of instances for one assets at a time
#' @param Tn Integer. The length of time horizon. Not including 0.
#' @param ini_W Scalar. The initial wealth.
#' @param discount Scalar. The rate of discount. e.g. 1/1.01
#' @param lambda The scalar that make the calculation easiler.
#' @param para Matrix. Initial value of all the paramters. Each column is a parameter and each row is a time point.
#' @author Yangzhuoran Yang
#' @seealse \code{round_EM}
#'
#'
#'  @export
EM <- function(Rr, Rf, max_iteration = 100, valuefunction = value_varmean, M = NROW(Rr[[1]]),
               Tn = length(Rr)-1,
               ini_W = 1000, discount=1/1.01,
               lambda = 1/2, para = NULL, early_stop = ceiling(max_iteration/5)){
  pm <- NULL

  pm[[1]] <- round_EM(Rr, Rf, valuefunction, M, Tn, ini_W, discount, lambda, para)
  pb <- txtProgressBar(min = 2, max = 100, style = 3)
  n_unchange <- 0
  for(it in 2:max_iteration){
    pm[[it]] <-  round_EM(Rr, Rf, valuefunction, M, Tn, ini_W, discount, lambda, para= pm[[it-1]][[1]])
    setTxtProgressBar(pb, it)

    if(abs(pm[[it]]$value[[1]] - pm[[it-1]]$value[[1]])/ abs(pm[[it-1]]$value[[1]]) <0.01){
      n_unchange <- n_unchange +1
      if(n_unchange == early_stop) break
    }

  }
  close(pb)
  return(pm)
}





