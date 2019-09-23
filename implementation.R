## ---- para ----

library(tidyverse)
Tn <- 5
N <- 100
M <- 1e4

Rf <- 1.01
discount <- 1/1.05

## ---- em ----


# Rr <- simulate_riskys()
# Rr <- step2(Rr)
# mean_Rr <- 1.05
# sd_Rr <- 0.03
# Rr <- replicate(Tn,  matrix(rnorm(M, mean_Rr, sd_Rr), ncol = 1), simplify = F)
# set.seed(22222)
# Rt <- sim_simple(Tn = Tn, N=100, M=M) %>%
#   lapply(function(x) x+1)
# saveRDS(Rt,"Rt.rds")
Rt <- readRDS("Rt.rds")
pm <- NULL
time <- numeric(102)
time[[1]] <- Sys.time()
pm[[1]] <- round_EM(Rt, Rf, valuefunction = value_varmean, discount = discount, lambda = 50)
time[[2]] <- Sys.time()
pb <- txtProgressBar(min = 1, max = 100, style = 3)
for(it in 2:100){
  # Rr <- replicate(Tn,  matrix(rnorm(M, mean_Rr, sd_Rr), ncol = 1), simplify = F)
  # Rt <- sim_simple(Tn = Tn, N=5, M=M)
  # Rr <- step2(Rr)
  pm[[it]] <- round_EM(Rt, Rf, para = pm[[it-1]][[1]], discount = discount, lambda = 50)
  time[[it+1]] <- Sys.time()
  setTxtProgressBar(pb, it)
}
time[[102]] <- Sys.time()
v <- sapply(pm, function(x) x[[2]][[1]])
saveRDS(pm, "em_assets_100_lambda50.rds")
ggplot2::qplot(y= v, x=seq_along(v), geom = "line")
plotly::ggplotly()

## ---- OCPA ----

Rt <- readRDS("Rt.rds")

test <- OCPA(Rt, Rf)
test2 <- OCPA(Rt, Rf, lambda_c=100, lambda_w=50)
View(OCPA(Rt, Rf, lambda_c=100, lambda_w=50))


pml2 <- round_EM(Rt, Rf, para = pml[[1]], discount = discount)
