

OCPA <- function(Rt, Rf, Tn = length(Rt)-1, M = NROW(Rt[[1]]),
                 ini_W = 1000, discount = 1/1.01,
                 lambda = 1/2){

  weights <- sapply(Rt, weights_lm)
  Rr <- mapply(function(Rt, omega) Rt %*% omega, Rt = Rt, omega = as.data.frame(weights))
  # bw <- apply(Rr, 2, function(x) diff(range(x))/(length(x)/10))
  bw <- apply(Rr-Rf, 2, bw.nrd0)


  J1 <- Rr-Rf


  mean_J <- lapply(seq_len(NCOL(Rr)),function(i) get_J(Rr[,i], Rf=Rf))
  mean_JJ <- lapply(seq_len(NCOL(Rr)),function(i) get_JJ(Rr[,i], bw[[i]], Rf=Rf))
  # mean_JJ <- lapply(seq_len(NCOL(Rr)),function(i) get_JJ(Rr[,i], 0.015))



  HDGF <- get_HDGF(Tn=Tn, lambda=lambda, mean_JJ = mean_JJ, mean_J = mean_J, Rf = Rf, discount = discount)
  Ht <- HDGF[[1]]
  Dt <- HDGF[[2]]
  Gt <- HDGF[[3]]
  Ft <- HDGF[[4]]

  Wt <- array(dim=dim(Rr))
  Wt[,1] <- ini_W-lambda

  Vt <- numeric(Tn)
  Zt <- list()
  for(i in seq_len(Tn)){
    Vt[[i]] <- get_Vt(Wt = Wt[[i]], Dt = Dt[[i]], Gt = Gt[[i]], Ft = Ft[[i]])
    Zt[[i]] <- get_Zt(Wt = Wt[[i]], Ht = Ht[[i]], Dt1 = Dt[[i+1]], Gt1 = Gt[[i+1]], Rf = Rf, lambda = lambda, mean_J = mean_J[[i]])
    # Wt[[i+1]] <- do.call(c,evo_Wt(Wt = Wt[[i]], Rf = Rf,lambda = lambda, Jt =Jt[[i]], Zt = Zt[[i]]))

    Wt[,i+1] <- evo_Wt(Wt = Wt[,i], Rf = Rf,lambda = lambda, J1 =J1[,i], Zt = Zt[[i]])
  }
  adi <- discount^(0:(Tn-1))
  cumu_adi <- rev(cumsum(adi))
  adii <- cumu_adi +rev(adi)
  # ft <- sapply(Vt, mean)-lambda^2*adii
  tt <- 0:(Tn-1)
  pow <- sapply(tt, function(x) seq(1, Tn - x, 1))
  adi <- lapply(pow, function(x) discount^x)
  cumu_adi <- sapply(adi, sum)
  adii <- cumu_adi + sapply(adi, function(x) tail(x, 1))
  adiil <- lambda^2 + adii
  # ft <- sapply(Vt, mean)-lambda^2*adii
  ft <- Vt-adiil

  return(list(ft = ft, Zt = Zt))

}


# get_Rr <- function(Rt){
#   weights <- sapply(Rt, weights_lm)
#   Rr <- mapply(function(Rt, omega) Rt %*% omega, Rt = Rt, omega = as.data.frame(weights))
#   return(Rr)
# }

get_J <- function(Rr, Rf){
  matrix(c(mean(Rr-Rf), -1))
}

get_JJ <- function(Rr, bw, Rf){
  y <- Rr-Rf
  two_one <- one_two <- mean(-y)
  # one_one <- mean(bw^2+(mean(y))^2)
  one_one <- bw^2+mean(y^2)
  matrix(c(one_one, one_two, two_one, 1 ), 2)
}

inverse_2by2 <- function(mat){
  (1/(mat[1,1]*mat[2,2]-mat[1,2]*mat[2,1])) * matrix(c(mat[2,2], -mat[1,2], -mat[2,1], mat[1,1]), byrow = T, nrow = 2)
}

# Calculate HDGF
get_HDGF <- function(Tn, lambda, mean_JJ, mean_J, Rf, discount){
  # H  0 : Tn-1 = Tn
  # D  0 : Tn   = Tn+1
  # G  0 : Tn   = Tn+1
  # F  0 : Tn   = Tn+1
  Ht <- vector("list", Tn)
  Dt <- vector("list", Tn+1)
  Gt <- vector("list", Tn+1)
  Ft <- vector("list", Tn+1)
  Dt[[Tn+1]] <- 1
  Gt[[Tn+1]] <- 0
  Ft[[Tn+1]] <- 0
  for(i in Tn:1){
    Ht[[i]] <- matrix(c(0,0,0,1), 2) + Dt[[i+1]]*mean_JJ[[i]]
    mul <- (1-Dt[[i+1]]*t(mean_J[[i]]) %*% inverse_2by2(Ht[[i]]) %*% mean_J[[i]])
    Dt[[i]] <- c(discount*Rf^2*Dt[[i+1]] * mul)
    Gt[[i]] <- c((discount*Rf*Gt[[i+1]] + discount*Rf*(Rf -2)*lambda*Dt[[i+1]]) * mul )
    Ft[[i]] <- c((discount*(Rf-2)^2*lambda^2*Dt[[i+1]] + 2*discount*(Rf-2)*lambda*Gt[[i+1]] + discount*Ft[[i+1]]) * mul)
  }
  return(list(Ht, Dt, Gt, Ft))

}

get_Vt <- function(Wt, Dt, Gt, Ft){
  Wt^2*Dt + 2*Wt*Gt + Ft
}
# Solution of control
get_Zt <- function(Wt, Ht, Dt1, Gt1, Rf, lambda, mean_J){
  # lapply(Wt,function(Wt) -(Dt1 * (Wt*Rf + (Rf-2)*lambda) +Gt1) * solve(Ht) %*% mean_J)
  -(Dt1 * (Wt*Rf + (Rf-2)*lambda) +Gt1) * inverse_2by2(Ht) %*% mean_J
  # -solve(Ht) %*% (Dt1 * (Wt*Rf + (Rf-2)*lambda)) %*% mean_J
}

evo_Wt <- function(Wt, Rf, lambda, J1, Zt){
  mapply(function(W, J){Jt <- matrix(c(J, -1)); c(W*Rf +(Rf-2)*lambda + t(Jt) %*% Zt)},
         W = Wt, J = J1)

  # mapply(function(Wt, Jt, Zt) c(Wt*Rf +(Rf-2)*lambda + t(Jt) %*% Zt),Wt = Wt, Jt = Jt, Zt = Zt, SIMPLIFY = FALSE)

}



