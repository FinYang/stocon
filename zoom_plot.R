
library(ggforce)



pm <- readRDS("em_single_asset.rds")
v <- sapply(pm, function(x) x[[2]])
v0 <- 3350.08

library(tidyverse)
data <- data.frame(x=seq_along(v), y=v/1000, vzoom = T)
ggplot(data, aes(y=y, x=x))+
  geom_line() +
  labs(x = "Number of Iterations",
       y = expression(paste("Value at time 0, ", V[0], "('000)"))) +
  geom_hline(yintercept = v0/1000, color = "blue") +
  geom_label(data = data.frame(x= c(100, 100), y=c(3.7, 2.7), label = c("PA", "EM"), vzoom = T),
            mapping = aes(x=x, y=y, label = label)) +
  facet_zoom(x = x >75, ylim = c(2,5), horizontal = F)

ggplot(mapping = aes(y=v/1000, x=seq_along(v)))+
  geom_line(show.legend = T) +
  labs(x = "Number of Iterations",
       y = expression(paste("Value at time 0, ", V[0], "('000)"))) +
  geom_hline(yintercept = v0/1000, color = "blue") +
  # geom_label(data = data.frame(x= c(100, 100), y=c(3.7, 2.7), label = c("PA", "EM"), vzoom = T),
            # mapping = aes(x=x, y=y, label = label)) +
  facet_zoom(x = x >75, ylim = c(2,5), horizontal = F, zoom.data = vzoom) +
  geom_label_repel(data = data.frame(x= c(100, 100), y=c(3.7, 2.7), label = c("PA", "EM"), vzoom = T), aes(label = labels))



