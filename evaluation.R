

new_stoconOCPA <- function(ft, BETA, C, weights, Wt, data, specification){
  structure(list(result=list(ft = ft, BETA = BETA, C=C, weights = weights, Wt = Wt),data = data, specification = specification), class = c("stoconOCPA"))
}
new_stoconEM <- function(x){
  structure(x, class = c("stoconEM"))
}

new_stoconMODEL <- function(x){
  structure(x, class = c(class(x), "stoconMODEL"))
}


validate_stoconMODEL <- function(model){
  if(!"stoconMODEL" %in% class(model))
    stop("Not a stoconMODEL")

  if(!identical(names(model), c("ft", "Zt", "weights", "Wt", "specification")))
    stop("Missing model results")

  if(!identical(names(model$specification), c("Tn", "M", "ini_W", "discount", "lambda_c", "lambda_w")))
    stop("Missing specifications")
  invisible(model)
}


fitted.stoconEM <- function(model, test_set, ...){
  pre <- model[[length(model)]]
  structure(list(result=list(ft = pre$value, Zt = cbind(pre$BETA, pre$C), weights = weights, Wt = Wt),data = data, specification = specification), class = c("stoconOCPA"))
}


fitted.stoconMODEL <- function(model, test_set, ...){
  weights <- model$result$weights
  Rf <- model$specification$Rf

  Rr <- mapply(function(Rt, omega) Rt %*% omega, Rt = test_set, omega = as.data.frame(weights))
  # J1 <- model$data$J1
  # mean_J <- model$data$mean_J
  # mean_JJ <- model$data$mean_JJ
  lambda_c <- model$specification$lambda_c
  lambda_w <- model$specification$lambda_w
  discount <- model$specification$discount
  Tn <- model$specification$Tn


  # HDGF <- get_HDGF(Tn=Tn, lambda_c=lambda_c, lambda_w = lambda_w, mean_JJ = mean_JJ, mean_J = mean_J, Rf = Rf, discount = discount)
  # Ht <- HDGF[[1]]
  # Dt <- HDGF[[2]]
  # Gt <- HDGF[[3]]
  # Ft <- HDGF[[4]]

  # Zt <- lapply(split(model$result$Zt, seq_len(NROW(model$result$Zt))), matrix)
  BETA <- model$result$BETA
  C <- model$result$C

  Wt <- array(dim=dim(Rr))
  Wt[,1] <- model$specification$ini_W-lambda_w
  for(t in 1:(Tn)){
    Wt[,t+1] <- (Wt[,t]-BETA[,t])*Rf + BETA[,t]*Rr[,t] - C[,t]
  }
  get_Value <- function(i, Wt){
    Value <- numeric(NROW(Rr))
    for(t in i:(Tn)){
      # BETA[[t-1+1]] <- beta_funct1on(para_beta = para_beta[t, ], w = W[,t])
      # C[[t-1+1]] <- c_funct1on(para_c = para_c[t, ], w = W[,t])
      Wt[,t+1] <- (Wt[,t]-BETA[,t])*Rf + BETA[,t]*Rr[,t] - C[,t]
      Value <- Value + discount^(t-i+1)*(C[,t]^2-2*lambda_c*C[,t])
    }
    Value <- Value + discount^(Tn-i+1)*(Wt[,Tn+1]^2-2*lambda_w*Wt[,Tn+1])
    Value <- mean(Value)

    return(Value)
  }
  value <- numeric(Tn)
  for(i in 1:(Tn)){
    value[[i]] <- get_Value(i, Wt = Wt)
    # value[[i]]$t <- i-1
  }

  # Vt <- numeric(Tn)
  # for(i in seq_len(Tn)){
  #   Vt[[i]] <- get_Vt(Wt = Wt[[i]], Dt = Dt[[i]], Gt = Gt[[i]], Ft = Ft[[i]])
  #   # Zt[[i]] <- get_Zt(Wt = Wt[[i]], Ht = Ht[[i]], Dt1 = Dt[[i+1]], Gt1 = Gt[[i+1]], Rf = Rf, lambda_c = lambda_c, lambda_w = lambda_w, mean_J = mean_J[[i]])
  #
  #   Wt[,i+1] <- evo_Wt(Wt = Wt[,i], Rf = Rf,lambda_c = lambda_c, lambda_w = lambda_w, J1 =J1[,i], Zt = Zt[[i]])
  # }
  return(value)
}
# fitted.stoconOCPA <- function(model, test_set, ...){
#   print(model$Zt)
# }


# length.stocon
print.stoconOCPA <- function(model){
  cat("stocon Model: OCPA\n")
  cat("Value function\n")
  cat(model$result$ft)

}


a <- EM(Rt, 1.01, lambda = 50)
b <- OCPA(train_set, 1.01, lambda_c = 50, lambda_w = 50)
fitted(model, test_set = test_set)









x <- data.frame(a=1:10)

deviden <- function(x, y){

  # abc <- list(...)
  x
}
deviden(x, w=3)

?mutate


mutate(x, bbc=11:20)

lapply(1:10, function(x, n, a) x+n+a, n=3, a=19)


