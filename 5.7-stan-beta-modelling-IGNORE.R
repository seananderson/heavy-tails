# cut from 5.5-modelling.R

source("5-shape-data.R")

gelm_scale <- function(x) {
  (x - mean(x, na.rm = TRUE)) / (2 * sd(x, na.rm = TRUE))
}
gomp_hat_base$log_Lifesp_scaled <- gelm_scale(log(gomp_hat_base$Lifesp))
gomp_hat_base$log_sigma_proc_50_scaled <- gelm_scale(log(gomp_hat_base$sigma_proc_50))
gomp_hat_base$b_50_scaled <- gelm_scale(gomp_hat_base$b_50)
gomp_hat_base$lambda_50_scaled <- gelm_scale(gomp_hat_base$lambda_50)
gomp_hat_base$log_dataset_length_scaled <- gelm_scale(log(gomp_hat_base$dataset_length))
gomp_hat_base$log_Len_scaled <- gelm_scale(log(gomp_hat_base$Len))
gomp_hat_base$heavy <- ifelse(gomp_hat_base$p10 > 0.5, 1, 0)

dat <- readRDS("gomp-base-mean-sd.rds")
look <- gpdd[,c("main_id", "taxonomic_class", "taxonomic_order", "exact_name")]
look <- look[!duplicated(look), ]
dat <- inner_join(dat, look)
dat <- inner_join(dat, gomp_hat_base[,c("main_id", "dataset_length",
    "log_dataset_length_scaled",
    "log_Lifesp_scaled", "p10")])

dat$log_sigma_scaled_mean <- gelm_scale(dat$mean_log_sigma_proc)
dat$log_sigma_scaled_sd <- dat$sd_log_sigma_proc / (2 * sd(dat$mean_log_sigma_proc))
dat$lambda_scaled_mean <- gelm_scale(dat$mean_lambda)
dat$lambda_scaled_sd <- dat$sd_lambda / (2 * sd(dat$mean_lambda))
dat$b_scaled_mean <- gelm_scale(dat$mean_b)
dat$b_scaled_sd <- dat$sd_b / (2 * sd(dat$mean_b))
dat <- na.omit(dat)
dat <- droplevels(dat)
dat$taxonomic_order <- as.factor(dat$taxonomic_order)
dat$taxonomic_class <- as.factor(dat$taxonomic_class)
dat$exact_name <- as.factor(dat$exact_name)
dat$order_id <- as.numeric(dat$taxonomic_order)
dat$sp_id <- as.numeric(dat$exact_name)
dat$class_id <- as.numeric(dat$taxonomic_class)

temp.dat <- gomp_hat_base[,c("log_Lifesp_scaled", "log_sigma_proc_50_scaled",
  "log_dataset_length_scaled", "b_50_scaled", "lambda_50_scaled",
  "taxonomic_class", "taxonomic_order", "heavy", "p10")]
temp.dat <- na.omit(temp.dat)
temp.dat$taxonomic_class <- as.factor(temp.dat$taxonomic_class)
temp.dat$taxonomic_order <- as.factor(temp.dat$taxonomic_order)

d <- temp.dat
d$order_id <- as.numeric(d$taxonomic_order)
d$class_id <- as.numeric(d$taxonomic_class)

if(!file.exists("stan-beta.rds")) {
  stan_beta <- stan_model("betareg.stan")
  saveRDS(stan_beta, file = "stan-beta.rds")
} else {
  stan_beta <- readRDS("stan-beta.rds")
}

m.stan.beta <- sampling(stan_beta,
  data = list(
    N = nrow(d),
    n_order = max(d$order_id),
    n_class = max(d$class_id),
    x = as.numeric(d$log_dataset_length_scaled),
    order_id = d$order_id,
    class_id = d$class_id,
    y = d$p10),
  pars = c("b", "mu_a", "sigma_a_class", "sigma_a_order", "phi"),
  iter = 2000, chains = 3)

sink("m.stan.beta.txt")
print(m.stan.beta)
sink()
saveRDS(m.stan.beta, file = "m.stan.beta.rds")

# test:
library(glmmADMB)
m.admb.2 <- glmmadmb(p10 ~ log_dataset_length_scaled + (1 | taxonomic_class / taxonomic_order), family = "beta", data = temp.dat)

# whole model with measurement error:
if(!file.exists("stan-beta2.rds")) {
  stan_beta2 <- stan_model("betareg2.stan")
  saveRDS(stan_beta2, file = "stan-beta2.rds")
} else {
  stan_beta2 <- readRDS("stan-beta2.rds")
}

m.stan.beta2 <- sampling(stan_beta2,
  data = list(
    N = nrow(d),
    n_order = max(d$order_id),
    n_class = max(d$class_id),
    x1 = as.numeric(d$lambda_50_scaled),
    x2 = as.numeric(d$b_50_scaled),
    x3 = as.numeric(d$log_sigma_proc_50_scaled),
    x4 = as.numeric(d$log_dataset_length_scaled),
    x5 = as.numeric(d$log_Lifesp_scaled),
    x1_sigma = rep(1.0, nrow(d)),
    x2_sigma = rep(1.0, nrow(d)),
    x3_sigma = rep(1.0, nrow(d)),
    order_id = d$order_id,
    class_id = d$class_id,
    y = d$p10),
  pars = c("b1", "b2", "b3", "b4", "b5", "mu_a",
    "sigma_a_class", "sigma_a_order", "phi"),
  iter = 400, chains = 3)

d <- dat

m.stan.beta3 <- sampling(stan_beta2,
  data = list(
    N = nrow(d),
    n_order = max(d$order_id),
    n_class = max(d$class_id),
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
    y = d$p10),
  pars = c("b1", "b2", "b3", "b4", "b5", "mu_a",
    "sigma_a_class", "sigma_a_order", "phi"),
  iter = 500, chains = 3)

#sink("m.stan.beta.txt")
#print(m.stan.beta)
#sink()
#saveRDS(m.stan.beta, file = "m.stan.beta.rds")

# compare to load("m.admb.2.rda") # bang on!

# now fit with true uncertainty around predictors:

stan_beta3 <- stan_model("betareg3.stan")
stan_beta4 <- stan_model("betareg4.stan")

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
    "a_class", "a_order", "mu"),
  iter = 2000, chains = 4)
saveRDS(m.stan.beta4, file = "beta-stan-samples.rds")
m <- readRDS("beta-stan-samples.rds")


means <- plyr::laply(extract(m), mean)[1:5]
ord <- order(means)

coefs <- c(expression(lambda), expression(b), expression(log(sigma)), expression(log(Timeseries~length)), expression(log(Lifespan)))

pdf("stan-beta-correlates.pdf", width = 4, height = 6)
par(mfrow = c(5, 1), mar = c(0,0,0,0), cex = 0.8, xpd = FALSE, oma = c(4, 4, 1, 1))
par(tck = -0.02, mgp = c(2, 0.5, 0), col.axis = "grey25", col = "grey25")
for(i in ord) {
  x <- density(extract(m)[[i]])
  plot(1, 1, xlim = c(-1, 1), ylim = c(0, 4),
    main = "", axes = FALSE, yaxs = "i", type = "n")
  polygon(c(x$x, rev(x$x)), c(x$y, rep(0, length(x$y))), border = "grey50", col = "grey90", lwd = 1)
  box()
  abline(v = 0, lty = 1, col = "grey60")
  legend(-1.2, 3.8, legend = coefs[i], bty = "n", adj = c(0, 0))
  axis(2, at = c(0, 2), las = 1)
}
axis(1)
mtext("Scaled coefficient value", side = 1, line = 2.2, cex = 0.8, outer = TRUE)
mtext("Probability density", side = 2, line = 2, cex = 0.8, outer = TRUE)
dev.off()

mu_a <- extract(m, pars = "mu_a")[[1]]
b4 <- extract(m, pars = "b4")[[1]]

x <- seq(-1, 4, length.out = 100)
plot(x, x, ylim = c(0, 1), type = "n")
for(i in sample(1:length(mu_a), 300)) {
  y_hat <- mu_a[i] + b4[i] * x
  lines(x, plogis(y_hat), col = "#00000020")
}
