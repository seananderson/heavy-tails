# cut from 5.5-modelling.R
# this version is a streamlined version with the final analysis
# see 5.6-... for the intermediate work

library("rstan")
options(mc.cores = 2L)
source("5-shape-data.R")

# scale by 2 standard deviations and subtract mean:
gelm_scale <- function(x) {
  (x - mean(x, na.rm = TRUE)) / (2 * sd(x, na.rm = TRUE))
}

ghb <- dplyr::select(gomp_hat_base, Lifesp, sigma_proc_50, b_50, lambda_50,
  dataset_length, p10, main_id)
ghb <- na.omit(ghb)
ghb$log_Lifesp_scaled <- gelm_scale(log(ghb$Lifesp))
ghb$log_sigma_proc_50_scaled <- gelm_scale(log(ghb$sigma_proc_50))
ghb$b_50_scaled <- gelm_scale(ghb$b_50)
ghb$lambda_50_scaled <- gelm_scale(ghb$lambda_50)
ghb$log_dataset_length_scaled <- gelm_scale(log(ghb$dataset_length))
ghb$heavy <- ifelse(ghb$p10 > 0.5, 1, 0)

dat <- readRDS("gomp-base-mean-sd.rds") # means and sds instead of medians
look <- gpdd[,c("main_id", "taxonomic_class", "taxonomic_order", "exact_name")]
look <- look[!duplicated(look), ]
dat <- inner_join(dat, look)
dat <- inner_join(dat, ghb[,c("main_id", "dataset_length",
    "log_dataset_length_scaled",
    "log_Lifesp_scaled", "p10")])

dat$log_sigma_scaled_mean <- gelm_scale(dat$mean_log_sigma_proc)
dat$log_sigma_scaled_sd <- dat$sd_log_sigma_proc / (2 * sd(dat$mean_log_sigma_proc))
dat$lambda_scaled_mean <- gelm_scale(dat$mean_lambda)
dat$lambda_scaled_sd <- dat$sd_lambda / (2 * sd(dat$mean_lambda))
dat$b_scaled_mean <- gelm_scale(dat$mean_b)
dat$b_scaled_sd <- dat$sd_b / (2 * sd(dat$mean_b))
dat <- droplevels(dat)
dat$taxonomic_order <- as.factor(dat$taxonomic_order)
dat$taxonomic_class <- as.factor(dat$taxonomic_class)
dat$exact_name <- as.factor(dat$exact_name)
dat$order_id <- as.numeric(dat$taxonomic_order)
dat$sp_id <- as.numeric(dat$exact_name)
dat$class_id <- as.numeric(dat$taxonomic_class)

d <- dat
saveRDS(d, file = "beta-modelling-dat.rds")

if(!file.exists("betareg4.rds")) {
  stan_beta4 <- stan_model("betareg4.stan")
  saveRDS(stan_beta4, file = "betareg4.rds")
} else {
  stan_beta4 <- readRDS("betareg4.rds")
}

if(!file.exists("beta-stan-samples.rds")) {
  m.stan.beta4 <- sampling(stan_beta4,
    data = list(
      N = nrow(d),
      n_order = max(d$order_id),
      n_class = max(d$class_id),
      n_sp = max(d$sp_id),
      x1 = as.numeric(d$lambda_scaled_mean),
      x2 = as.numeric(d$b_scaled_mean),
      x3 = as.numeric(d$log_sigma_scaled_mean),
      x4 = as.numeric(d$log_dataset_length_scaled),
      x5 = as.numeric(d$log_Lifesp_scaled),
      x1_sigma = as.numeric(d$lambda_scaled_sd),
      x2_sigma = as.numeric(d$b_scaled_sd),
      x3_sigma = as.numeric(d$log_sigma_scaled_sd),
      order_id = d$order_id,
      class_id = d$class_id,
      sp_id = as.numeric(d$sp_id),
      y = d$p10),
    pars = c("b1", "b2", "b3", "b4", "b5", "mu_a",
      "sigma_a_class", "sigma_a_order", "sigma_a_sp", "phi",
      "a_class", "a_order"),
    iter = 2000, chains = 4, thin = 1, control = list(adapt_delta = 0.999))
  png("trace-beta.png", width = 700, height = 1000)
  rstan::traceplot(m.stan.beta4, pars = c("b1", "b2", "b3", "b4", "b5",
      "mu_a", "sigma_a_class", "sigma_a_order", "sigma_a_sp", "phi"),
    nrow = 5, ncol = 2,
    inc_warmup = FALSE)
  dev.off()
  saveRDS(m.stan.beta4, file = "beta-stan-samples.rds")
}

m <- readRDS("beta-stan-samples.rds") # or reload this
sink("beta-stan-samples-2.txt")
print(m)
sink()

means <- plyr::laply(extract(m), mean)[1:5]
ord <- order(means)

coefs <- c(expression(lambda), expression(b), expression(log(sigma)),
  expression(log(Timeseries~length)), expression(log(Lifespan)))

pdf("stan-beta-correlates.pdf", width = 4, height = 6)
par(mfrow = c(5, 1), mar = c(0,0,0,0), cex = 0.9, xpd = FALSE,
  oma = c(4, 4, 1, 1))
par(tck = -0.02, mgp = c(2, 0.5, 0), col.axis = "grey25", col = "grey25")
for(i in ord) {
  x <- density(extract(m)[[i]])
  plot(1, 1, xlim = c(-1, 1), ylim = c(0, 4),
    main = "", axes = FALSE, yaxs = "i", type = "n")
  polygon(c(x$x, rev(x$x)), c(x$y, rep(0, length(x$y))), border = "grey50",
    col = "grey90", lwd = 1.5)
  box()
  abline(v = 0, lty = 1, col = "grey60")
  legend(-1.2, 3.8, legend = coefs[i], bty = "n", adj = c(0, 0))
  axis(2, at = c(0, 2), las = 1)
}
axis(1)
mtext("Coefficient value\n(per 2 SDs of predictor)", side = 1, line = 2.4,
  cex = 0.9, outer = TRUE)
mtext("Probability density", side = 2, line = 1.8, cex = 0.9, outer = TRUE)
dev.off()

mu_a <- extract(m, pars = "mu_a")[[1]]
b4 <- extract(m, pars = "b4")[[1]]

# and get predictions for log sigma
b3 <- extract(m, pars = "b3")[[1]]

x <- seq(min(d$log_sigma_scaled_mean), max(d$log_sigma_scaled_mean),
  length.out = 100)
zsd <- 2 * sd(dat$mean_log_sigma_proc)
x_raw <- exp(x * zsd + mean(dat$mean_log_sigma_proc))
pred_sigma <- data.frame(sigma_proc_mean = x_raw, log_sigma_scaled_mean = x)

p <- matrix(ncol = length(x), nrow = length(mu_a))
for(i in 1:length(mu_a)) {
  p[i,] <- plogis(mu_a[i] + b3[i] * x)
}
pred_sigma <- data.frame(pred_sigma, plyr::adply(p, 2, quantile,
  probs = c(0.05, 0.5, 0.95)))

# and get % increase in Pr(nu < 10) with an increase in 10 data points:
x <- seq(min(d$log_dataset_length_scaled), max(d$log_dataset_length_scaled),
  length.out = 4000)
zsd <- 2 * sd(log(ghb$dataset_length))
x_raw <- exp(x * zsd + mean(log(ghb$dataset_length)))
ticks <- c(20, 40, 60, 80, 100, 120)
pred_dataset_length <- data.frame(dataset_length = x_raw, log_dataset_length_scaled = x)
p <- matrix(ncol = length(x), nrow = length(mu_a))
for(i in 1:length(mu_a)) {
  p[i,] <- plogis(mu_a[i] + b4[i] * x)
}
pred_dataset_length <- data.frame(pred_dataset_length,
  plyr::adply(p, 2, quantile, probs = c(0.05, 0.5, 0.95)))

p_inc <- plyr::ldply(seq(20, 110, 10), function(i) {
 p <- dplyr::filter(pred_dataset_length, dataset_length > i)[1,"X50."]
 data.frame(N = i, p = p)
 })
p_perc_inc <- vector(length = (nrow(p_inc) - 1))
for(i in 1:(nrow(p_inc)-1)) {
  p_perc_inc[i] <- p_inc$p[i+1] / p_inc$p[i]
}

# pdf("perc-prob-increase-with-n.pdf", width = 4.5, height = 3.7)
# par(mar = c(4, 6, 1.5, 1), oma = c(0, 0, 0, 0), cex = 0.8)
# par(tck = -0.015, mgp = c(2.8, 0.4, 0), col.axis = "grey25", col = "grey25", las = 1)
# plot(p_inc$N[-length(p_inc$N)], 100 * (p_perc_inc - 1), xlab = "",
#   ylab = paste("Expected percent increase in probability of detecting\nheavy",
#     "tails with 10 additional time steps of data"),
#   pch = 21, col = "grey30", bg = "grey80", las = 1,
#   ylim = c(0, max(100 * (p_perc_inc-1))))
# mtext("Initial time-series length", side = 1, line = 2, cex = 0.8)
# dev.off()
#
saveRDS(p_inc, file = "prob_inc_heavy_with_n.rds")
