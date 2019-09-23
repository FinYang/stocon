library(tidyverse)
library(lubridate)
library(naniar)
russell <- read_csv("data-raw/russell-200.csv", col_types = paste0("c",paste0(rep("d", 195), collapse = ""))) %>%
  mutate(Date = dmy(Date))
 # vis_miss(russell, sort_miss=TRUE, warn_large_data = F)

# russell %>% apply(2, function(x) diff(x, 250))

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
array2list <- function(xarray){
  plyr:::splitter_a(xarray,3, .id = NULL)
}

simple_return <- function(x, p = 250){
  diff(x, 250)/x[1:(length(x)-250)]
}

rolling_window <- function(x, n=10, p=250){
  sapply(1:n, function(nc) x[seq(nc,(nc+p-1),1)])
}
  # sapply(1:(1000-250+1), function(nc) seq(nc,(nc+250-1),1)) %>% View

russellR <- russell[,-1] %>%
  apply(2, simple_return)

p_na <- russellR %>%
  apply(2, function(x) sum(is.na(x))/length(x))
russellR <- russellR[,which(p_na<0.5)] %>%
  as.data.frame() %>%
  drop_na()




Rt <- lapply(seq_len(NCOL(russellR)), function(i) rolling_window(russellR[,i]) ) %>%
  list2array() %>%
  aperm(c(1,3,2)) %>%
  array2list()





# vis_miss(as.data.frame(russellR), sort_miss=TRUE, warn_large_data = F)




# add_yearly <- function(xcol, w = 250){
#   sapply(seq_along(xcol), function(nc) sum(xcol[nc:(nc+w-1)]))
# }

#
#
#
# NA_ind <- russell %>% apply(2, anyNA) %>% which()
# russell <- select(russell, -NA_ind)


# russellR <- russell[,-1] %>%
#   transmute_all(function(x) c(diff(log(x)), NA)) %>%
#   mutate(Date = russell$Date) %>%
  # drop_na() %>%
  # filter(year(Date) %in% 2012:2017 ) %>%
  # split(year(.$Date)) %>%
  # lapply(function(x) select(x, -Date) %>%
  #          (function(y){y+1})() %>%
  #          (function(z){z[1:250,]})) %>%
  # lapply(as.matrix)
# test <- OCPA(russellR, Rf = 1.00001, ini_W = 10000)
