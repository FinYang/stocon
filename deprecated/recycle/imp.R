# Rr <- Rt <- sim_simple(N=100)
pm <- EM(Rt, 1.01, lambda = 50)
pa <- OCPA(Rt,1.01, lambda_c = 50, lambda_w = 50)
pa$ft
v <- sapply(pm, function(x) x[[2]][[1]])
ggplot2::qplot(y= v, x=seq_along(v), geom = "line")
plotly::ggplotly()

v


pa$Zt
pm[[length(pm)]]$BETA
pm[[length(pm)]]$C

# -------------------------------------------------------------------------



OCPA_analytic <- function(Rt, Rf, Tn = length(Rt)-1, M = NROW(Rt[[1]]),
                 ini_W = 1000, discount = 1/1.01,
                 lambda_c = 1/2, lambda_w = 1/2){

  weights <- sapply(Rt, weights_lm)
  Rr <- mapply(function(Rt, omega) Rt %*% omega, Rt = Rt, omega = as.data.frame(weights))
  # bw <- apply(Rr, 2, function(x) diff(range(x))/(length(x)/10))
  bw <- apply(Rr-Rf, 2, bw.nrd0)


  J1 <- Rr-Rf


  mean_J <- lapply(seq_len(NCOL(Rr)),function(i) get_J(Rr[,i], Rf=Rf))
  mean_JJ <- lapply(seq_len(NCOL(Rr)),function(i) get_JJ(Rr[,i], bw[[i]], Rf=Rf))
  # mean_JJ <- lapply(seq_len(NCOL(Rr)),function(i) get_JJ(Rr[,i], 0.015))



  HDGF <- get_HDGF(Tn=Tn, lambda_c=lambda_c, lambda_w = lambda_w, mean_JJ = mean_JJ, mean_J = mean_J, Rf = Rf, discount = discount)
  Ht <- HDGF[[1]]
  Dt <- HDGF[[2]]
  Gt <- HDGF[[3]]
  Ft <- HDGF[[4]]

  Wt <- array(dim=dim(Rr))
  Wt[,1] <- ini_W-lambda_w

  Vt <- numeric(Tn)
  Zt <- list()
  for(i in seq_len(Tn)){
    Vt[[i]] <- get_Vt(Wt = Wt[[i]], Dt = Dt[[i]], Gt = Gt[[i]], Ft = Ft[[i]])
    Zt[[i]] <- get_Zt(Wt = Wt[[i]], Ht = Ht[[i]], Dt1 = Dt[[i+1]], Gt1 = Gt[[i+1]], Rf = Rf, lambda_c = lambda_c, lambda_w = lambda_w, mean_J = mean_J[[i]])

    Wt[,i+1] <- evo_Wt(Wt = Wt[,i], Rf = Rf,lambda_c = lambda_c, lambda_w = lambda_w, J1 =J1[,i], Zt = Zt[[i]])
  }
  # adi <- discount^(0:(Tn-1))
  # cumu_adi <- rev(cumsum(adi))
  # adii <- cumu_adi +rev(adi)
  tt <- 0:(Tn-1)
  pow <- sapply(tt, function(x) seq(1, Tn - x, 1))
  adi <- lapply(pow, function(x) discount^x)
  cumu_adi <- sapply(adi, sum)
  adii <- cumu_adi*lambda_c^2 + lambda_w^2*sapply(adi, function(x) tail(x, 1))
  # adiil <- lambda^2 + adii
  # ft <- sapply(Vt, mean)-lambda^2*adii
  ft <- Vt-adii
  Zt_c <- lapply(Zt, function(x) c(x[[1]],x[[2]]+lambda_c))
  return(list(ft = ft, Zt = Zt_c, weights = weights, Wt = Wt+lambda_w))

}
