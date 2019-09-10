library(tidyverse)
library(lubridate)
library(naniar)
russell <- read_csv("data-raw/russell-200.csv", col_types = paste0("c",paste0(rep("d", 195), collapse = ""))) %>%
  mutate(Date = dmy(Date))
 # vis_miss(russell, sort_miss=TRUE, warn_large_data = F)
NA_ind <- russell %>% apply(2, anyNA) %>% which()
russell <- select(russell, -NA_ind)


russellR <- russell[,-1] %>%
  transmute_all(function(x) c(diff(x), NA)/x) %>%
  drop_na()
