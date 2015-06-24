## plot various skewness figures

# look at the skew-t distribution:
library(skewt)
library(ggplot2)
params <- expand.grid(x = seq(-7, 7, length.out = 100), df = c(2, 5, 1e6),
  gamma = round(exp(c(-0.6, -0.3, 0, 0.3, 0.6)), 1))
o <- plyr::mdply(params, dskt)

o$nu_label <- paste("nu =", o$df)
o$skew_label <- paste("skew =", o$gamma)

o$nu_label <- reorder(o$nu_label, o$df)

plain_theme <- theme_bw() + theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  strip.background = element_blank())

p <- ggplot(o, aes(x, V1, colour = nu_label)) + geom_line() +
  facet_wrap(~skew_label, nrow = 1) + ylab("Probability density") +
  xlab("x") + plain_theme + labs(colour=expression(nu))
  print(p)
ggsave("skew-t-illustration.pdf", width = 9, height = 2)

source("5-shape-data.R")

gomp_hat_skew <- gomp_hat_skew %>% group_by(taxonomic_class) %>%
  arrange(log_skew_50) %>%
  mutate(plot_id = 1:n()) %>% as.data.frame()

p <- ggplot(gomp_hat_skew, aes(1/nu_50, log_skew_50)) +
  geom_segment(aes(y = log_skew_25, yend = log_skew_75, x = 1/nu_50, xend = 1/nu_50),
    alpha = 0.08, lwd = 0.3) +
  #geom_segment(aes(y = log_skew_5, yend = log_skew_95, x = nu_50, xend = nu_50),
    #alpha = 0.06, lwd = 0.2) +
  geom_segment(aes(y = log_skew_50, yend = log_skew_50, x = 1/nu_25, xend = 1/nu_75),
    alpha = 0.09, lwd = 0.3) +
  geom_point(alpha = 0.3, pch = 21) +
  theme_bw() + geom_hline(yintercept = 0, lty = 2) +
  scale_x_continuous(breaks = 1/c(100, 10, 5, 3, 2), labels = c(100, 10, 5, 3, 2)) +
  geom_hline(yintercept = median(gomp_hat_skew$log_skew_50), col = "red") +
  xlab(expression(nu)) + ylab("log(Skewness parameter)")
#print(p)
ggsave("skewness-vs-nu.pdf", width = 6, height = 5)

