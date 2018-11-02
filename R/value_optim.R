#' Value Function
#'
#'
#'
#'
#' @param updatetheta updated theta
#' @param i starting period
#' @param theta matrix contain all theta
#' @param Rt returns
#' @param M number of realization
#' @param Tn Number of periods, excluding time 0.
#' @param weights weight matrix where each column contain weights for one period
#' @param W wealth
#' @param ValueOnly logical. If FALSE, returns both value and wealth
#' @param utility utility function
#' @param beta discount factor. Currently set to be the same 1/1.05
#' @param consu_function consumption has been set to inverse logit
#'
#' @return list contain value and wealth (optional)
#' @importFrom boot inv.logit
#' @author Yangzhuoran Yang
#' @export
#'
value_func <- function(updatetheta, i, theta,
                       Rt, M = NROW(Rt[[1]]), Tn = length(Rt) - 1,
                       weights, W = NULL, ValueOnly = FALSE, utility = power_u,
                       beta = 1/1.05, consu_function = boot::inv.logit){

  if(i==1 & length(W)==0){
    W <- matrix(nrow = Tn+1, ncol = M)
    W[1,] <- 1000
  }



  weights_list <- weights %>% data.frame() %>% as.list()

  theta[i,] <- updatetheta
  Value <- numeric(M)
  consu <- matrix(nrow = Tn, ncol = M)

  for(t in i:Tn){
    consu[t,] <- W[t,]*consu_function(theta[t,1]+theta[t,2]*W[t,])
    W[t+1,] <- (W[t,]-consu[t,])*(1+ t(Rt[[t]] %*% weights_list[[t]]) )
    Value <- Value + utility(consu[t,])*beta^(t-1)
  }
  Value <- sum(Value + utility(W[NROW(W),])*beta^Tn)
  if(ValueOnly) return(Value)
  return(list(W=W,Value=Value))
}

neg_value_func <- function(updatetheta, i, theta,
                           Rt, M = NROW(Rt[[1]]), Tn = length(Rt) - 1,
                           weights, W = NULL, ValueOnly = T, utility = power_u,
                           beta = 1/1.05, consu_function = inv.logit){
  -value_func(updatetheta, i, theta,
              Rt, M , Tn ,
              weights, W , ValueOnly, utility,
              beta, consu_function)
}


#' Dynamic optimization
#'
#'
#' @param Rt returns
#' @param weights weight matrix where each column contain weights for one period
#' @param W_ini initial wealth currently set to be the same 1000
#' @param utility utility function
#' @param beta discount factor. Currently set to be the same 1/1.05
#' @param consu_function consumption has been set to inverse logit
#'
#' @return list contain theta and value
#' @importFrom boot inv.logit
#'
#' @author Yangzhuoran Yang
#' @export
#'
optim_dynam <- function(Rt, M = NROW(Rt[[1]]), Tn = length(Rt) - 1,
                        weights,  W_ini = 1000, utility = power_u,
                        beta = 1/1.05, consu_function = inv.logit){
  # wealth ----

  W <- matrix(nrow = Tn+1, ncol = M)
  W[1,] <- W_ini # initial wealth currently set to be the same


  theta  <- numeric(Tn*2) %>% matrix(ncol = 2)

  #####Valuefunction

  W_V <- value_func(0, 1, theta, Rt, M=M, Tn=Tn,  weights, W, utility = utility,
                    beta = beta, consu_function = consu_function)
  W <- W_V$W

  for(t in Tn:1){
    # theta(t,:)=fminunc(@(updatetheta)value_func(updatetheta,t, theta, R,weight, M, 10, W),theta(t,:));
    theta[t,] <- optim(theta[t,], neg_value_func, i=t, theta = theta,
                       Rt = Rt, weights = weights, W=W,
                       ValueOnly = TRUE, utility = utility, beta = beta,
                       consu_function = consu_function)$par
    # update?
    # W <- value_func(theta[t,], i=t, theta = theta,
    #                 Rt = Rt, weights = weights, W=W, ValueOnly = F, utility = utility)$W
  }
  w_value <- value_func(theta[1,], 1, theta, Rt, M=M, Tn=Tn,  weights, W,
             utility = utility, ValueOnly = F, beta = beta, consu_function = consu_function)
  return(list(theta = theta, v = w_value$Value, W = w_value$W))
}




# old V----
# Vt <- matrix(nrow=M, ncol=Tn+1)
#
# Vt[,Tn+1] <- beta^Tn * u(W[,Tn+1])
# v_cmu <- Vt[,Tn+1]
# for(i in Tn:1){
#
#   V <- function(thet, v_cmu){
#     consum <- W[,i]*inv.logit(thet[1]+thet[2]*Rtw[,i])
#     (mean(beta^i*u(consum)) + mean(v_cmu)) %>% return()
#   }
#   theta[i,] <- optim(c(0,0), V, v_cmu=v_cmu)$par
#   for(l in i:Tn){
#     consu[,l] <- W[,l]*inv.logit(theta[l,1]+theta[l,2]*Rtw[,l])
#     W[,l+1] <- (1+Rtw[,l])*(W[,l]-consu[,l])
#     Vt[,l] <- beta^(l-1)*u(consu[,l])
#   }
#   Vt[,Tn+1] <- beta^Tn * u(W[,Tn+1])
#
#   v_cum <- rowSums(Vt[,i:(Tn+1)])
# }


