
library(ggforce)
library(tidyverse)


pm <- readRDS("em_single_asset.rds")
upb <- 70
v <- sapply(pm, function(x) x[[2]])
v0 <- 3350.08
v <- v[1:upb]
data <- data.frame(x=seq_along(v), y=v/1000, labels = NA)
data2 <- data.frame(x= rep(length(v)+1, 2), y=c(3.55, 2.9), labels = c("PA", "EM"), vzoom = T)

ggplot(data, aes(y=y, x=x))+
  geom_line() +
  labs(x = "Number of Iterations",
       y = expression(paste("Value at time 0, ", V[0], "('000)"))) +
  facet_zoom(x = x >3, ylim = c(2,5),horizontal = F,  zoom.data = vzoom) +
  geom_label(data = data2,
            mapping = aes(label = labels))+
  geom_hline(mapping = aes(yintercept = v0/1000), color = "blue")
