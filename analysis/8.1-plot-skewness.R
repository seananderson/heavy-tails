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
  # print(p)
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

## extract samples from all posteriors of skewness parameter:

sk <- readRDS("skew_samples.rds")
sk <- sk %>% group_by(main_id) %>%
  mutate(median_nu = median(nu)) %>%
  mutate(h1 = ifelse(median_nu <= 10, "1 heavy",
      ifelse(median_nu <= 70, "2 moderate", "3 normal"))) %>%
  as.data.frame()

normal <- filter(sk, h1 == "3 normal")$log_skew %>% sort
moderate <- filter(sk, h1 == "2 moderate")$log_skew %>% sort
heavy <- filter(sk, h1 == "1 heavy")$log_skew %>% sort

ticks <- c(0.1, 0.2, 0.5, 1, 2, 5)
# pdf("skewness-densities-base.pdf", width = 5, height = 4)
# par(yaxs = "i")
pal <- c(RColorBrewer::brewer.pal(6, "YlOrRd")[c(5,4)], "#4D4D4D")
dn <- density(normal, from = min(sk$log_skew), to = max(sk$log_skew))
dm <- density(moderate, from = min(sk$log_skew), to = max(sk$log_skew))
dh <- density(heavy, from = min(sk$log_skew), to = max(sk$log_skew))
# plot(1, 1, xlim = c(-2, 2), ylim = c(0, 1.0), type = "n", axes = FALSE, ylab = "", xlab = "")
pfunc <- function(x, y, col, alpha = 50, fill = TRUE) {
  if (fill) {
    polygon(c(x, rev(x)), c(y, rep(0, length(y))), border = NA,
      col = paste0(col, alpha))
  }
  lines(x, y, col = col, lwd = 1.8)
}
# pfunc(dn$x, dn$y, col = pal[3], alpha = "15")
# pfunc(dm$x, dm$y, col = pal[2], alpha = "15")
# pfunc(dh$x, dh$y, col = pal[1], alpha = "15")
# axis(1, at = log(ticks), ticks)
# abline(v = 0, lty = 2)
# dev.off()

group_by(sk, h1) %>%
  summarise(p_lt0 = sum((log_skew < 0))/n(),
    q5 = exp(quantile(log_skew, prob = 0.05)),
    q50 = exp(quantile(log_skew, prob = 0.5)),
    q95 = exp(quantile(log_skew, prob = 0.95)))

params <- expand.grid(x = seq(-7, 7, length.out = 200), df = c(5, 20, 1e6),
  gamma = c(0.65, 0.75, 0.85))
o <- plyr::mdply(params, dskt) %>%
  rename(dens = V1)
on <- data.frame(x = unique(params$x))
on$dens <- dnorm(on$x)

pfunc2 <- function(nu, skew, col, label1 = "", label2 = "") {
  p <- filter(o, df == nu & gamma == skew) %>%
    arrange(x)
  plot(1, 1, xlim = c(-7, 7), ylim = c(0, 0.4), type = "n",
    axes = FALSE, ylab = "", xlab = "")
  lines(on$x, on$dens, col = paste0(pal[3], 40), lwd = 1.8, lty = 1)
  lines(p$x, p$dens, col = col, lwd = 1.8)
  mtext(label1, side = 1, col = "grey40", line = 0.1, cex = 0.85)
  mtext(label2, side = 1, col = "grey40", line = 1.6, cex = 0.85)
}

add_label <- function(xfrac = -0.02, yfrac = -0.1, label = "", pos = 4, ...) {
  par(xpd = NA)
  u <- par("usr")
  x <- u[1] + xfrac * (u[2] - u[1])
  y <- u[4] - yfrac * (u[4] - u[3])
  text(x, y, label, pos = pos, ...)
  par(xpd = FALSE)
}

#dev.off()

# 1. I would make all the axes numbers smaller
# 2. Add some more white space between a and b, between the left and the
# right, and below each of the sparkline plots.
# 3. Need at least three numbers on each y-axis for the sparklines since
# in log-scale, but they could be much smaller in size.
# 4. Red grouse looks weird, why does the red lower dashed line end up
# *above* the black dashed lower line in the last year? It might be that
# you need to run the MCMC for a much longer run and retain a much
# higher number of samples (100,000) to get an accurate estimate of the
# lower 99.5-ile. Or might just be the true behavior given the higher
# red median projected.
# 5. The text needs to mention the numbers that are the ratio of red
# dashed to black dashed after five years for at least these five cases.
# It's hard to infer from the log-scale of the plot but looks
# potentially substantial.
# 6. in a the nu is listed first, in c the globulin (my default name
# when I don't know the greek letter name) is listed first and nu
# second.

l <- rbind(
  c(1, 2, 3, 5, 5, 5),
  c(4, 4, 4, 6, 6, 6),
  c(4, 4, 4, 7, 7, 7),
  c(4, 4, 4, 8, 8, 8)
  )
pdf("skew-fig.pdf", width = 5, height = 4)
par(mar = c(1,0,2.5,1), oma = c(3.7, 1.5, .5, 0), cex = 0.6)
par(tck = -0.02)
par(mgp = c(2, 0.8, 0))
layout(l)
#layout.show(8)

# illustrations of skew-nu vs. normal:
pfunc2(nu = 5, skew = 0.65, col = pal[1], label1 = expression(nu == 5), label2 = expression(gamma == 0.65))
add_label(label = "A", cex = 1.8, font = 2, yfrac = -0.64)
add_label(label = "Process deviation shape", cex = 1.4, font = 1,
  xfrac = 0.2, col = "grey20", yfrac = -0.64)
pfunc2(nu = 20, skew = 0.75, col = pal[2], label1 = expression(nu == 20), label2 = expression(gamma == 0.75))
pfunc2(nu = 1e6, skew = 0.85, col = pal[3], label1 = expression(nu == infinity), label2 = expression(gamma == 0.85))

# posteriors of skewness:
par(mar = c(0, 0, 5, .5))
plot(1, 1, xlim = c(-2.5, 2), ylim = c(0, 1.0), type = "n",
  axes = FALSE, ylab = "", xlab = "", yaxs = "i")
#abline(v = 0, lty = 1, col = "#00000070")
segments(0, 0, 0, 0.75, lty = 1, col = "#00000070")
pfunc(dn$x, dn$y, col = pal[3], alpha = "15")
pfunc(dm$x, dm$y, col = pal[2], alpha = "15")
pfunc(dh$x, dh$y, col = pal[1], alpha = "15")
par(mgp = c(2, 0.5, 0))
par(tck = -0.02)
axis(1, at = log(ticks), ticks, cex.axis = 1.1, col.axis = "grey40", col = "grey40")
mtext(expression(Skewness~parameter~(gamma)), side = 1, col = "grey20", line = 2.5, cex = 0.9)
#mtext("Skewed", side = 1, col = "grey40", line = -6.7, adj = 0.05)
#mtext("down", side = 1, col = "grey40", line = -5.3, adj = 0.05)

#mtext("Skewed up", side = 1, col = "grey40", line = -6.7, adj = 0.87)
# mtext("up", side = 1, col = "grey40", line = -5.3, adj = 0.87)

add_label(label = "B", cex = 1.8, font = 2, yfrac = -0.05)
add_label(label = "Skewness parameter posteriors", cex = 1.4, font = 1,
  yfrac = -0.05, xfrac = 0.05, col = "grey20")

par(xpd = NA)
text(log(1.2), 0.63, "Normal tails", col = pal[3], pos = 4, cex = 1.1)
text(log(1), 0.73, "Slightly heavy tails", col = pal[2], pos = 4, cex = 1.1)
text(log(0.79), 0.83, "Heavy-tailed populations", col = pal[1], pos = 4, cex = 1.1)
par(xpd = FALSE)

# spark lines with projections:
# from 8-sim-extinction-risk.R:
# 20579 - grey heron - severe winters
# 10113 willow grouse - parasites and predators
# 7099 hare - but re-run with - nu 27, skew 0.5
# red grouse - 10039 parasites and predators

ids <- c(10113, 7099, 1235, 20579)
n_project <- 6L - 1L

get_gomp_res <- function(pop, sigma_proc, lambda, b) {
  res <- rep(NA, length(pop))
  for(i in 2:length(pop)) {
    res[i] <- log(pop)[i] - (lambda + b * log(pop[i-1]))
  }
   l <- qnorm(0.0001, 0, sd = sigma_proc)
   u <- qnorm(0.9999, 0, sd = sigma_proc)
  #l <- sigma_proc * -4
  #u <- sigma_proc * 4
  return(list(res = res, l = l, u = u))
}

heavy_res <- plyr::ldply(ids, function(i) {
  x <- filter(gpdd, main_id == i)$population_untransformed
  sig <- filter(gomp_hat_base, main_id == i)$sigma_proc_50
  lambda <- filter(gomp_hat_base, main_id == i)$lambda_50
  b <- filter(gomp_hat_base, main_id == i)$b_50
  qq <- get_gomp_res(x, sigma_proc = sig, lambda = lambda, b = b)
  with(qq, data.frame(main_id = i, res = c(res, rep(NA, n_project)), l, u,
      pop = c(x, rep(NA, n_project)), year = c(-(length(x) - 2):0, 1:(n_project + 1))))
})

heavy_res <- mutate(heavy_res, bs_down = ifelse(res < l, TRUE, FALSE))
heavy_res <- mutate(heavy_res, bs_up = ifelse(res > u, TRUE, FALSE))

# merge in the projections:
before_and_projections <- readRDS("before_and_projections.rds")
plot_dat <- before_and_projections %>% arrange(main_id, year) %>%
  filter(main_id %in% ids)

# ratios of projections:
# qqq <- before_and_projections %>% group_by(model, main_id) %>% summarise(lower = q_low3[year == 5]) %>% as.data.frame %>% group_by(main_id) %>% summarise(ratio = lower[model == "heavy-skew"] / lower[model == "normal"]);ggplot(qqq, aes(log(ratio))) + geom_density() + xlim(-2.5, 2.5)
# round(100 - exp(quantile(qqq$ratio, probs = c(0.25, 0.5, 0.75))) * 100, 0)

plot_dat2 <- inner_join(plot_dat, heavy_res)

par(mar = c(0, 4.5, 2.5, 0))
par(mgp = c(2, 0.35, 0))

spark <- function(dat) {
  ylim <- log(range(na.omit(c(dat$q_low3, dat$med))))
  plot(1, 1, xlim = c(-54, 6), ylim = ylim,
    type = "n", xlab = "", ylab = "", axes = FALSE)

  rect(1, min(ylim), 6, max(ylim), col = "grey90", border = NA)

  dat_d <- filter(dat, year <= 1, model == "normal")
  dat_pr_n <- filter(dat, year >= 1, model == "normal")
  dat_pr_h <- filter(dat, year >= 1, model == "heavy-skew")

  lines(dat_d$year, log(dat_d$med), col = "grey20")
  points(dat_d$year[dat_d$bs_down], log(dat_d$med[dat_d$bs_down]), col = pal[1], pch = 20)
  points(dat_d$year[dat_d$bs_up], log(dat_d$med[dat_d$bs_up]), col = "blue", pch = 20)
  if (min(ylim) < log(2000)) {
    tick_ <- c(1, 10, 100, 1000, 10000)
  } else {
    tick_ <- c(1000, 3000, 5000, 9000)
  }
  axis(2, las = 1, col = "grey40", col.axis = "grey40", at = log(tick_), labels = tick_, cex.axis = 1.1)

  # projection:
  lines(dat_pr_n$year, log(dat_pr_n$med), col = "grey20")
  lines(dat_pr_n$year, log(dat_pr_n$q_low3), col = "grey20", lty = "22")

  lines(dat_pr_h$year, log(dat_pr_h$med), col = pal[1])
  lines(dat_pr_h$year, log(dat_pr_h$q_low3), lty = "22", col = pal[1])
}

ii <<- 0
for (i in ids) {
  ii <<- ii + 1
  this_dat <- filter(plot_dat2, main_id == i)
  spark(this_dat)
  if (ii == 1) {
    add_label(label = "C", cex = 1.8, font = 2, yfrac = -0.5, xfrac = -0.16)
    add_label(label = "Example time series", cex = 1.4, font = 1,
      xfrac = -0.08, col = "grey20", yfrac = -0.5)
  }
  meta <- filter(gomp_hat_skew, main_id == unique(this_dat$main_id))
  sk_param <- sprintf("%.1f", round(exp(meta$log_skew_50), 1))
  nu_param <- sprintf("%.0f", round(meta$nu_50, 0))
  #add_label(label = paste0(meta$common_name, ", sk ", sk_param, ", nu ", nu_param),
    #cex = 1.1, col = "grey20")
  add_label(label = bquote(.(meta$common_name)~gamma == .(sk_param)~nu == .(nu_param)),
    cex = 1.1, col = "grey20")
}

par(tcl = -0.22)
par(mgp = c(2, 0.35, 0))
  if (ii == 4)
    axis(1, col = "grey40", col.axis = "grey40", cex.axis = 1.1)

mtext("Years", side = 1, line = 2.2, col = "grey20", cex = 0.9)

dev.off()

## #########################
## # test other ways of getting residual black swans:
## id <- 20579
## id <- 10113
## id <- 7099
## id <- 12 # normal
## # pdf("bs-threshold-methods.pdf", width = 7, height = 12)
## par(mfrow = c(1, 2))
## # cn <- filter(gomp_hat_normal, main_id == 20579)
## # resids <- get_gomp_res(filter(gpdd, main_id == 20579)$population_untransformed, sigma_proc = cn$sigma_proc_50, lambda = cn$lambda_50, b = cn$b_50)$res
## #
## # hist(resids[-1], breaks = 20)
## # x <- seq(-0.5, 0.5, length.out = 200)
## # lines(x, dnorm(x, mean = 0, sd = cn$sigma_proc_50)*5.5, col = "red")
## # bs <- qnorm(0.001, 0, sd = cn$sigma_proc_50)
## # bs <- cn$sigma_proc_50 * -2
## # abline(v = bs, col = "blue")
## # plot(resids, type = "o")
## # abline(h = bs, col = "blue")
##
## # with tails:
## cn <- filter(gomp_hat_skew, main_id == id)
## cn <- filter(gomp_hat_skew, main_id == id)
## resids <- get_gomp_res(filter(gpdd, main_id == id)$population_untransformed, sigma_proc = cn$sigma_proc_50, lambda = cn$lambda_50, b = cn$b_50)$res
## bs1 <- qnorm(0.0001, 0, sd = cn$sigma_proc_50)
## bs2 <- cn$sigma_proc_50 * -5
## bs2 <- c(bs2, -bs2)
## bs1 <- c(bs1, -bs1)
## # par(mfrow = c(1, 2))
## hist(resids[-1], breaks = 20, xlim = range(c(bs1, bs2, resids[-1])))
## # x <- seq(-0.5, 0.5, length.out = 200)
## # lines(x, dnorm(x, mean = 0, sd = cn$sigma_proc_50)*3.5, col = "red")
## abline(v = bs1, col = "blue")
## abline(v = bs2, col = "purple")
## plot(resids, type = "o", ylim = range(c(bs1, bs2, resids[-1])))
## abline(h = bs1, col = "blue")
## abline(h = bs2, col = "purple")
##
##
## # with tails, but re-fit normal to residuals:
## cn <- filter(gomp_hat_base, main_id == 20579)
## resids <- get_gomp_res(filter(gpdd, main_id == 20579)$population_untransformed, sigma_proc = cn$sigma_proc_50, lambda = cn$lambda_50, b = cn$b_50)$res
## # par(mfrow = c(1, 2))
## hist(resids[-1], breaks = 20)
## x <- seq(-0.5, 0.5, length.out = 200)
## norm_sd <- sd(resids[-1])
## lines(x, dnorm(x, mean = 0, sd = norm_sd)*5.5, col = "red")
## bs <- qnorm(0.001, 0, sd = norm_sd)
## abline(v = bs, col = "blue")
## plot(resids, type = "o")
## abline(h = bs, col = "blue")
##
## # or by taking a value that's at the 1% of the t distribution:
## cn <- filter(gomp_hat_base, main_id == 20579)
## resids <- get_gomp_res(filter(gpdd, main_id == 20579)$population_untransformed, sigma_proc = cn$sigma_proc_50, lambda = cn$lambda_50, b = cn$b_50)$res
## # par(mfrow = c(1, 2))
## hist(resids[-1], breaks = 20)
## bs <- qt(0.01, df = cn$nu_50) * cn$sigma_proc_50
## abline(v = bs, col = "blue")
## plot(resids, type = "o")
## abline(h = bs, col = "blue")
