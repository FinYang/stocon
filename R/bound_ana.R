#
# Rt <- sim_simple(dis_par = list(mu=1:5, vol=5:1, df=2), distribution = "t", dependent = TRUE) %>%
#   list2array() %>% aperm(c(2,3,1)) %>% array2list()
# Rt <- Rt[[1]]
# lbound_for_each <- function(Rt, Rf, Tn = NCOL(Rt)-1, N = NROW(Rt),
#                             ini_W = 1000, discount = 1/1.01,
#                             lambda_c = 1/2, lambda_w = 1/2, lambda = NULL){
#   J <- Rt-Rf
#
#   get_J_ <- function(J){
#     matrix(c(J, -1))
#   }
#
#   get_JJ_ <- function(J){
#     J %*% t(J)
#   }
#
#   J <- lapply(seq_len(NCOL(J)),function(i) get_J_(J[,i]))
#   JJ <- lapply(J, get_JJ_)
#
#   INN <- matrix(0, N+1, N+1)
#   INN[N+1, N+1] <- 1
#      det(INN + JJ[[Tn]])
#
#
#   get_HDGF_ <- function(Tn, lambda_c, lambda_w, JJ, J, Rf, discount, N){
#     # H  0 : Tn-1 = Tn
#     # D  0 : Tn   = Tn+1
#     # G  0 : Tn   = Tn+1
#     # F  0 : Tn   = Tn+1
#     Ht <- vector("list", Tn)
#     Dt <- vector("list", Tn+1)
#     Gt <- vector("list", Tn+1)
#     Ft <- vector("list", Tn+1)
#     Dt[[Tn+1]] <- 1
#     Gt[[Tn+1]] <- 0
#     Ft[[Tn+1]] <- 0
#     INN <- matrix(0, N+1, N+1)
#     INN[N+1, N+1] <- 1
#     for(i in Tn:1){
#       Ht[[i]] <- INN + Dt[[i+1]]*JJ[[i]]
#       mul <- (1-Dt[[i+1]]*t(J[[i]]) %*% solve(Ht[[i]]) %*% J[[i]])
#       Dt[[i]] <- c(discount*Rf^2*Dt[[i+1]] * mul)
#       Gt[[i]] <- c((discount*Rf*Gt[[i+1]] + discount*Rf*((Rf -1)*lambda_c-lambda_w)*Dt[[i+1]]) * mul )
#       Ft[[i]] <- c((discount*((Rf -1)*lambda_c-lambda_w)^2*Dt[[i+1]] + 2*discount*((Rf -1)*lambda_c-lambda_w)*Gt[[i+1]] + discount*Ft[[i+1]]) * mul)
#     }
#     return(list(Ht, Dt, Gt, Ft))
#
#   }
#
#
#   HDGF <- get_HDGF_(Tn=Tn, lambda_c=lambda_c, lambda_w = lambda_w, JJ = JJ, J = J, Rf = Rf, discount = discount, N=N)
#   Ht <- HDGF[[1]]
#   Dt <- HDGF[[2]]
#   Gt <- HDGF[[3]]
#   Ft <- HDGF[[4]]
#
#
#
#
#
#
#
#
#
# }
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
