library(tidyverse)
library(lubridate)
library(naniar)
russell <- read_csv("data-raw/russell-200.csv", col_types = paste0("c",paste0(rep("d", 195), collapse = ""))) %>%
  mutate(Date = dmy(Date))
 # vis_miss(russell, sort_miss=TRUE, warn_large_data = F)
NA_ind <- russell %>% apply(2, anyNA) %>% which()
russell <- select(russell, -NA_ind)


russellR <- russell[,-1] %>%
  transmute_all(function(x) c(diff(log(x)), NA)) %>%
  mutate(Date = russell$Date) %>%
  drop_na() %>%
  filter(year(Date) %in% 2012:2017 ) %>%
  split(year(.$Date)) %>%
  lapply(function(x) select(x, -Date) %>%
           (function(y){y+1})() %>%
           (function(z){z[1:250,]})) %>%
  lapply(as.matrix)
# test <- OCPA(russellR, Rf = 1.00001, ini_W = 10000)
