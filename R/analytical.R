

#' @export
OCPA_analytical <- function(N, Rf, Rt = NULL, varcov = NULL, distribution = c("norm", "t"),
                            dis_par = list(mu = 0.05, vol = 0.02, df=2),par = 0.2,  Tn = 9, M = 10000,
                            ini_W = 1000, discount = 1/1.01,
                            lambda_c = 1/2, lambda_w = 1/2, lambda = NULL, ...){
  if(!is.null(lambda)){
    lambda_c <- lambda_w <- lambda
  }
  if((!is.null(Rt))&&(!is.null(varcov))){
    N <- NCOL(Rt[[1]])
    Tn <- length(Rt)-1
    M <- NROW(Rt[[1]])

  } else {
    tem <- sim_simple(dis_par = dis_par, distribution = distribution, N=N, Tn=Tn, M=M, return_varcov = TRUE, ...)
    Rt <- tem[[1]]
    varcov <- tem[[2]]

  }

  weights <- sapply(Rt, weights_lm)

  Rr <- mapply(function(Rt, omega) Rt %*% omega, Rt = Rt, omega = as.data.frame(weights))
  Rr <- Rr[,-NCOL(Rr)]
  J1 <- Rr-Rf
  mu <- dis_par$mu
  vol <- dis_par$vol
  if ((length(mu) == 1) && (length(vol) == 1)) {
    mu <- rep(mu, N)
    vol <- rep(vol, N)
  }
  Rmu <- sapply(seq_len(NCOL(weights)), function(i) t(weights[,i]) %*% matrix(mu))
  Rvar <- sapply(seq_len(NCOL(weights)), function(i) t(weights[,i]) %*% varcov %*% weights[,i])



  get_J <- function(Rmu, Rf){
    matrix(c(Rmu-Rf, -1))
  }

  get_JJ <- function(Rmu, Rvar, Rf){
    two_one <- one_two <- Rf-Rmu
    # one_one <- mean(bw^2+(mean(y))^2)
    one_one <- Rvar +Rmu^2 -2*Rmu*Rf +Rf^2
    matrix(c(one_one, one_two, two_one, 1 ), 2)
  }
  mean_J <- lapply(seq_along(Rmu),function(i) get_J(Rmu[[i]], Rf=Rf))
  mean_JJ <- lapply(seq_along(Rmu),function(i) get_JJ(Rmu[[i]], Rvar[[i]], Rf=Rf))



  HDGF <- get_HDGF(Tn=Tn, lambda_c=lambda_c, lambda_w = lambda_w, mean_JJ = mean_JJ, mean_J = mean_J, Rf = Rf, discount = discount)
  Ht <- HDGF[[1]]
  Dt <- HDGF[[2]]
  Gt <- HDGF[[3]]
  Ft <- HDGF[[4]]

  Wt <- matrix(ncol = Tn+1, nrow = M)
  Wt[,1] <- ini_W-lambda_w

  Vt <- list()
  Zt <- list()
  for(i in seq_len(Tn)){
    Vt[[i]] <- get_Vt(Wt = Wt[,i], Dt = Dt[[i]], Gt = Gt[[i]], Ft = Ft[[i]])
    Zt[[i]] <- get_Zt(Wt = Wt[,i], Ht = Ht[[i]], Dt1 = Dt[[i+1]], Gt1 = Gt[[i+1]], Rf = Rf, lambda_c = lambda_c, lambda_w = lambda_w, mean_J = mean_J[[i]])
    # Vt[[i]] <- get_Vt(Wt = Wt[[i]], Dt = Dt[[i]], Gt = Gt[[i]], Ft = Ft[[i]])
    # Zt[[i]] <- get_Zt(Wt = Wt[[i]], Ht = Ht[[i]], Dt1 = Dt[[i+1]], Gt1 = Gt[[i+1]], Rf = Rf, lambda_c = lambda_c, lambda_w = lambda_w, mean_J = mean_J[[i]])

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

  # ft <- Vt-adii
  ft <- mapply(function(Vt, adii)Vt-adii, Vt=Vt, adii=adii, SIMPLIFY = F  ) %>% sapply(mean)

  BETA <- sapply(Zt,function(z) sapply(z, function(x) x[[1]]))
  C <- sapply(Zt,function(z) sapply(z, function(x) x[[2]] +lambda_c))

  # Zt_c <- lapply(Zt, function(x) c(x[[1]],x[[2]]+lambda_c))
  # Zt_c <- do.call(base::rbind,Zt_c)
  # colnames(Zt_c) <- c("BETA", "C")
  specification <- list(Rf=Rf, Tn = Tn, M = M,
                        ini_W = ini_W, discount = discount,
                        lambda_c = lambda_c, lambda_w = lambda_w)
  data <- list(Rt=Rt, J1=J1, mean_J = mean_J, mean_JJ = mean_JJ)
  # new_stoconMODEL(ft, Zt_c, weights, Wt+lambda_w)
  return(new_stoconMODEL(new_stoconOCPA(ft, BETA, C, weights[,-NCOL(weights)], Wt+lambda_w, data, specification)))

}

