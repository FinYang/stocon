
#' @export
list2array <- function(xlist){
  d1 <- sapply(xlist, NROW) %>%
    unique()
  if(length(d1) != 1) stop("Different row number")
  d2 <- sapply(xlist, NCOL) %>%
    unique()
  if(length(d2) != 1) stop("Different col number")
  xlist %>%
    unlist() %>%
    array(dim = (c(d1, d2, length(xlist))))
}
#' @export
array2list <- function(xarray){
  plyr:::splitter_a(xarray,3, .id = NULL)
}
#' @export
L_bound <- function(Rt, Rf, M = NROW(Rt[[1]]),
                    Tn = length(Rt)-1,
                    ini_W = 1000, discount=1/1.01,
                    lambda = 1/2, polynomial = TRUE, parallel = TRUE){
  N <- NCOL(Rt[[1]])

  Rt_array <- list2array(Rt)
  Rt_list <- array2list(aperm(Rt_array, c(3,2,1)))

  if(polynomial){


    init_para <- function(){
      para <- matrix(0,nrow = Tn+1, ncol = 5)
      para <- cbind(para,matrix(1/N , nrow = Tn+1, ncol =  N-1))
      return(c(para))
    }

    value_for_each <- function(para,x){
      # W <- matrix(nrow = M, ncol = Tn+1)
      W <- numeric(Tn+1)
      W[[1]] <- ini_W
      para <- matrix(c(para), nrow= Tn+1)
      para_beta <- para[ ,1:3]
      para_c <- para[ ,4:5]
      if(NCOL(x) > 1){
        para_w <- para[,6:(NCOL(para))]
        para_w <- cbind(para_w, apply(para_w, 1, function(x) 1- sum(x)))
        Rr <- mapply(function(r, w) r %*% as.matrix(w), r=split(x, 1:NROW(x)), w=split(para_w, 1:NROW(para_w)))
      }
      Value <- numeric(Tn+1)
      # t-1 : Tn-1
      BETA <- numeric(Tn)
      C <- numeric(Tn)
      for(j in 1:Tn){
        BETA[[j]] <- beta_function(para_beta = para_beta[j, ], w = W[[j]])
        C[[j]] <- c_function(para_c = para_c[j, ], w = W[[j]])
        W[[j+1]] <- (W[[j]]-BETA[[j]])*Rf + BETA[[j]]*Rr[[j]] - C[[j]]
        Value[[j]] <- discount^(j)*(C[[j]]^2-2*lambda*C[[j]])

      }
      Value[[Tn+1]] <- discount^(Tn)*(W[[Tn+1]]^2-2*lambda*W[[Tn+1]])
      out <- sum(Value)
      return(out)

    }
  } else {
    init_para <- function(){
      para <- matrix(0,nrow = Tn+1, ncol = 2)
      para <- cbind(para,matrix(1/N , nrow = Tn+1, ncol =  N-1))
      return(c(para))
    }

    value_for_each <- function(para,x){
      # W <- matrix(nrow = M, ncol = Tn+1)
      W <- numeric(Tn+1)
      W[[1]] <- ini_W
      para <- matrix(c(para), nrow= Tn+1)
      # para_beta <- para[ ,1:3]
      # para_c <- para[ ,4:5]
      if(NCOL(x) > 1){
        para_w <- para[,3:(NCOL(para))]
        para_w <- cbind(para_w, apply(para_w, 1, function(x) 1- sum(x)))
        Rr <- mapply(function(r, w) r %*% as.matrix(w), r=split(x, 1:NROW(x)), w=split(para_w, 1:NROW(para_w)))
      }
      Value <- numeric(Tn+1)
      # t-1 : Tn-1
      BETA <- para[,1, drop = TRUE]
      C <- para[,2, drop = TRUE]
      for(j in 1:Tn){
        # BETA[[j]] <- beta_function(para_beta = para_beta[j, ], w = W[[j]])
        # C[[j]] <- c_function(para_c = para_c[j, ], w = W[[j]])
        W[[j+1]] <- (W[[j]]-BETA[[j]])*Rf + BETA[[j]]*Rr[[j]] - C[[j]]
        Value[[j]] <- discount^(j)*(C[[j]]^2-2*lambda*C[[j]])

      }
      Value[[Tn+1]] <- discount^(Tn)*(W[[Tn+1]]^2-2*lambda*W[[Tn+1]])
      out <- sum(Value)
      return(out)

    }
  }
  optim_for_each <- function(x, max_iteration=100){
    para <- list()
    para[[1]] <- init_para()
    var_path <- numeric(max_iteration)
    path <- list()
    n_rep <- 0
    for(i in seq_len(max_iteration)){
      path[[i]] <- optim(para[[i]], value_for_each, x=x)
      para[[i+1]] <- path[[i]]$par
      var_path[[i]] <- path[[i]]$value
      if(i>1 && (((var_path[[i]]-var_path[[i-1]])/var_path[[i]])<0.01) ){
        n_rep <- n_rep+1
      }
      if(n_rep>3) break
    }
    optim(para[[length(para)]], value_for_each, x=x)$value
    # return(var_path)
  }
  if(parallel){

    old_plan <- future::plan(future::multiprocess)
    on.exit(future::plan(old_plan))
    out <- furrr::future_map(Rt_list, optim_for_each, .progress = TRUE)
  } else {

    out <- pbapply::pbsapply(Rt_list, optim_for_each)
  }
  mean(do.call(base::c, out))
}


















