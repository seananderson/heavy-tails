# Plot various skewness figures

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
  geom_segment(aes(y = log_skew_50, yend = log_skew_50, x = 1/nu_25, xend = 1/nu_75),
    alpha = 0.09, lwd = 0.3) +
  geom_point(alpha = 0.3, pch = 21) +
  theme_bw() + geom_hline(yintercept = 0, lty = 2) +
  scale_x_continuous(breaks = 1/c(100, 10, 5, 3, 2), labels = c(100, 10, 5, 3, 2)) +
  geom_hline(yintercept = median(gomp_hat_skew$log_skew_50), col = "red") +
  xlab(expression(nu)) + ylab("log(Skewness parameter)")
ggsave("skewness-vs-nu.pdf", width = 6, height = 5)

## extract samples from all posteriors of skewness parameter:

sk <- readRDS("skew_samples.rds")
sk <- sk %>% group_by(main_id) %>%
  mutate(
    p10 = sum(nu <= 10) / n(),
    p50 = sum(nu <= 10) / n(),
    median_nu = median(nu),
    p10_gt50 = ifelse(p10 > 0.5, TRUE, FALSE),
    p50_gt50 = ifelse(p50 > 0.2, TRUE, FALSE),
    heavy = ifelse(p10_gt50, "1 heavy", ifelse(p50_gt50, "2 slighly heavy", "3 normal")),
    heavy3 = ifelse(p10_gt50, "1 slightly heavy/heavy",
      ifelse(p50_gt50, "1 slightly heavy/heavy", "3 normal")),
    heavy2 = ifelse(p10_gt50, "1 heavy",
      ifelse(p50_gt50, "3 slightly heavy/normal", "3 slightly heavy/normal")),
    h1 = ifelse(median_nu <= 10, "1 heavy",
      ifelse(median_nu <= 70, "2 moderate", "3 normal"))) %>%
  as.data.frame()

pal <- c(RColorBrewer::brewer.pal(5, "YlOrRd")[c(4,2)], "#000000")
library(ggplot2)
gg_skew <- function(heavy_type) {
  ggplot(sk, aes_string("log_skew", fill = heavy_type, col = heavy_type)) +
    geom_density(alpha = 0.2, lwd = 0.5) +
    scale_x_continuous(breaks = log(c(0.2, 0.5, 1, 2, 5)), labels = c(0.2, 0.5, 1, 2, 5)) +
    coord_cartesian(xlim = c(-2.3, 2.3)) + xlab("Skew parameter") + plain_theme +
    geom_vline(xintercept = 0, col = "#00000050") +
    scale_fill_manual(values = pal) +
    scale_colour_manual(values = pal)
}
p1 <- gg_skew("heavy")
p2 <- gg_skew("heavy2")
p3 <- gg_skew("heavy3")
p4 <- gg_skew("h1")
pdf("skewness-densities.pdf", width = 10, height = 10)
gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2)
dev.off()

pdf("skewness-densities1.pdf", width = 5, height = 4)
p4 %>% print()
dev.off()

normal <- filter(sk, h1 == "3 normal")$log_skew %>% sort
moderate <- filter(sk, h1 == "2 moderate")$log_skew %>% sort
heavy <- filter(sk, h1 == "1 heavy")$log_skew %>% sort

ticks <- c(0.2, 0.5, 1, 2, 5)
par(yaxs = "i")
pal <- c(RColorBrewer::brewer.pal(6, "YlOrRd")[c(5,4)], "#4D4D4D")
dn <- density(normal, from = min(sk$log_skew), to = max(sk$log_skew))
dm <- density(moderate, from = min(sk$log_skew), to = max(sk$log_skew))
dh <- density(heavy, from = min(sk$log_skew), to = max(sk$log_skew))
plot(1, 1, xlim = c(-2, 2), ylim = c(0, 1.0), type = "n", axes = FALSE, ylab = "", xlab = "")
pfunc <- function(x, y, col, alpha = 50) {
  polygon(c(x, rev(x)), c(y, rep(0, length(y))), border = NA,
    col = paste0(col, alpha))
  lines(x, y, col  = col, lwd = 1.8)
}
pfunc(dn$x, dn$y, col = pal[3], alpha = "15")
pfunc(dm$x, dm$y, col = pal[2], alpha = "15")
pfunc(dh$x, dh$y, col = pal[1], alpha = "15")
axis(1, at = log(ticks), ticks)

group_by(sk, heavy) %>%
  summarise(p_lt0 = sum((log_skew < 0))/n(),
    q5 = exp(quantile(log_skew, prob = 0.05)),
    q50 = exp(quantile(log_skew, prob = 0.5)),
    q95 = exp(quantile(log_skew, prob = 0.95)))

group_by(sk, heavy2) %>%
  summarise(p_lt0 = sum((log_skew < 0))/n())

group_by(sk, heavy3) %>%
  summarise(p_lt0 = sum((log_skew < 0))/n())

group_by(sk, heavy) %>%
  summarise(quantile(log_skew))

d1 <- density(heavy,
  from = min(c(heavy, normal)),
  to   = max(c(heavy, normal)))
d2 <- density(normal,
  from = min(c(heavy, normal)),
  to   = max(c(heavy, normal)))
z <- d2$y/sum(d2$y) - d1$y/sum(d1$y)
plot(d1$x, z, type = "l")
abline(h=0)
sum(z[z>0])

