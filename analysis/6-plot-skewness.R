# plot various skewness figures

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
  xlab("x") + plain_theme + labs(colour='')
  print(p)
ggsave("skew-t-illustration.pdf", width = 9, height = 2)


