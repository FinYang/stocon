
library(ggforce)
library(tidyverse)


pm <- readRDS("em_single_asset.rds")
upb <- 70
v <- sapply(pm, function(x) x[[2]])
v0 <- 3350.08
v <- v[1:upb]
data <- data.frame(x=seq_along(v), y=v/1000, labels = NA, oury = v0/1000)
data2 <- data.frame(x= rep(length(v)+1, 2), y=c(3.55, 2.9), labels = c("PA", "EM"), vzoom = T)

ggplot(data, aes(y=y, x=x))+
  geom_line(linetype = 5) +
  geom_line(aes(y = oury)) +
  labs(x = "Number of Iterations",
       y = expression(paste("Value at time 0, ", V[0], "('000)"))) +
  facet_zoom(x = x >3, ylim = c(2,5),horizontal = F,  zoom.data = vzoom) +
  geom_label(data = data2,
             mapping = aes(label = labels))#+
# geom_hline(mapping = aes(yintercept = v0/1000, linetype = "dash"))








data <- mapply(function(end, ab) data.frame(em = end[[1]]/1000, pa = end[[2]]/1000, mean_Rr = ab[[1]], sd_Rr = ab[[2]], x=1:20),
               end = end, ab = ablist, SIMPLIFY = F)
data <- do.call(rbind, data) %>%
  gather(key = "algorithm", value = "value", em:pa)%>%
  # mutate(label = paste0(expression(mu), "=", mean_Rr-1, expression(sigma), "=", sd_Rr))
  mutate(label_mu = paste0("=",round(mean_Rr-1, 2)), label_sd = paste0("=", sd_Rr), labels = paste0("mu", label_mu,", sigma",label_sd))
text_data <- mapply(function(end, ab) data.frame(em = end[[1]][[20]]/1000, pa = end[[2]]/1000, mean_Rr = ab[[1]], sd_Rr = ab[[2]], x=20),
                    end = end, ab = ablist, SIMPLIFY = F)
text_data <- do.call(rbind, text_data) %>%
  mutate(mean_Rr_1 = round(mean_Rr-1, 2)) %>%
  # mutate(label = paste0(expression(mu), "=", mean_Rr-1, expression(sigma), "=", sd_Rr))
  mutate(label_mu = paste0("=",round(mean_Rr-1, 2)), label_sd = paste0("=", sd_Rr), labels = paste0("* mu", label_mu,", *sigma",label_sd))
# ggplot(filter(data, mean_Rr == ablist[[1]][[1]], sd_Rr == ablist[[1]][[2]] )) +
# ggplot(filter(data,  sd_Rr == ablist[[1]][[2]])) +
ggplot(data) +
  geom_line(aes(x=x, y=value, linetype = algorithm, color = paste("mu ==", round(mean_Rr-1, 2),  "~sigma==", sd_Rr))) +
  scale_color_hue(guide = FALSE)+
  scale_linetype_manual(name = "Algorithm", breaks = c("pa", "em"), labels = c("PA", "EM"), values = c(5,1)) +
  theme(legend.position = c(0.8,0.7)) +
  geom_text(data = text_data, parse = TRUE, mapping = aes(x=x, y=pa, label = paste("mu ==", round(mean_Rr-1, 2),  "~sigma==", sd_Rr),
                                                          color =paste("mu ==", round(mean_Rr-1, 2),  "~sigma==", sd_Rr) ))+
  labs(x = "Number of Iterations",
       y = expression(paste("Value at time 0, ", V[0], "('000)")))
