






# Calculate HDGF
HDGF <- function(Tn, lambda, mean_JJ, mean_J, Rf){
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
    mul <- (1-Dt[[i+1]]*t(mean_J[[i]]) %*% solve(Ht[[i]]) %*% mean_J[[i]])
    Dt[[i]] <- c(discount*Rf^2*Dt[[i+1]] * mul)
    Gt[[i]] <- c((discount*Rf*Gt[[i+1]] + discount*Rf*(Rf -2)*lambda*Dt[[i+1]]) * mul )
    Ft[[i]] <- c((discount*(Rf-2)^2*lambda^2*Dt[[i+1]] + 2*discount*(Rf-2)*lambda*Gt[[i+1]] + discount*Ft[[i+1]]) * mul)
  }
  return(list(Ht, Dt, Gt, Ft))

}


