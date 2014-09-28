library(rstan)
stan_gomp2_ar1 <- readRDS("stan-gomp2-ar1.rds")
stan_gomp2_ar1_obs <- readRDS("stan-gomp2-ar1-obs.rds")
stan_gomp2_ar1_obs_equal <- readRDS("stan-gomp2-ar1-obs-equal.rds.rds")
stan_gomp2_ar1_obs_est <- readRDS("stan-gomp2-ar1-obs-est.rds")

# Comparisons:
# Base model used in main figures:
# Gompertz AR1, sigma_proc estimated, sigma_obs assumed = 0
#
# Sensitivity to measurement error assumptions:
# Gompertz AR1, sigma_proc estimated, sigma_obs assumed = 0.1
# Gompertz AR1, sigma_proc estimated, sigma_obs assumed = 0.3
# Gompertz AR1, sigma_proc and sigma_obs estimated and assumed equal
#
# Alternative functional form:
# Ricker-logistic AR1 [need to code], sigma_proc estimated,
# sigma_obs assumed = 0
#
# Very simple model (same as b = 1, phi = 0):
# Model growth rates as stationary, sigma_proc estimated, sigma_obs assumed = 0
#
# Plan: small Bash script that sources the ID scripts and then submits 2 or 3
# jobs for each running script. Could also specify the assumed observation
# error this way too. Crank up the samples a bit.
# Make a cross plot figure showing each of these against the base model or
# each of these against each of these.

# gpdd <- readRDS("gpdd-clean.rds")
# x <- subset(gpdd, main_id == 20527)
# pop <- x$population_untransformed
# y <- diff(log(x$population_untransformed))
# par(mfrow = c(2, 1))
# plot(log(pop), type = "o")
# plot(y, type = "o")
# stan_gomp_ar1 <- readRDS("stan-gomp-ar1.rds")
# stan_gomp <- readRDS("stan-gomp.rds")
  warmup <- 400
  iterations <- 800
  # z <- seq(-1, 1, length.out = 100);plot(z, dcauchy(z, 0, 1), type = "l", ylim = c(0, 0.4))

  ###
  # OR simulate something
  # set.seed(1)
  m <- list()

  load("nu_effective_seeds.rda")

  for(j in 1:1) {
    set.seed(nu_3_seeds_N50$seeds[j])
    phi_sim <- 0.2
    nu_sim <- 3
    sigma_sim <- 0.65
    eps_sim <- rep(NA, 50)
    eps_sim[1] <- metRology::rt.scaled(1, df = nu_sim, mean = 0, sd = sigma_sim)
    for(i in 2:length(eps_sim)) {
      eps_sim[i] <- eps_sim[i-1] * phi_sim +
        metRology::rt.scaled(1, df = nu_sim, mean = 0, sd = sigma_sim)
    }
    # optional obs. error:
    eps_sim <- eps_sim + rnorm(length(eps_sim), 0, 0.3)
    #plot(eps_sim)
    lamb_sim <- 1.5
    b_sim <- 0.5
    y1 <- 3
    pop <- rep(NA, length(eps_sim))
    pop[1] <- y1
    for(i in 2:length(pop)) {
      pop[i] <- lamb_sim + b_sim * pop[i-1] + eps_sim[i]
    }
    pop <- exp(pop)
    par(mfrow = c(1, 1))
    # plot(log(pop), type = "l")
    x <- data.frame(population_untransformed = pop)
    ###

    m[[j]] <- sampling(stan_gomp2_ar1_obs_equal,
      data = list(N = nrow(x), y = log(x$population_untransformed),
        nu_rate = 0.01,
        b_lower = -1, b_upper = 2),
      pars = c("lambda", "b", "phi", "nu", "sigma_obs_proc"), iter = iterations,
      chains = 4, warmup = warmup)


    #   phi_out[j] <- median(extract(sm2_ar1, pars = "phi")[[1]])
    #   l[j] <- quantile(extract(sm2_ar1, pars = "phi")[[1]], probs = 0.05)
    #   u[j] <- quantile(extract(sm2_ar1, pars = "phi")[[1]], probs = 0.95)
  }
##############
out <- plyr::ldply(m, function(x) {
  m <- as.data.frame(lapply(extract(x), quantile, probs = 0.5))
  m$lp__ <- NULL
  #   l <- as.data.frame(lapply(extract(x), quantile, probs = 0.05))
  #   u <- as.data.frame(lapply(extract(x), quantile, probs = 0.95))
  m
})
out <- reshape2::melt(out)
true_df <- reshape2::melt(data.frame(lambda = lamb_sim, sigma_proc = sigma_sim,
  b = b_sim, phi = phi_sim, nu = nu_sim))

# ggplot(out, aes(variable, value)) + geom_violin(scale = "width") + facet_wrap(~variable, scales = "free") + geom_point(alpha = 0.5, col = "red") + geom_hline(data = true_df, aes(yintercept = value), lty = 2)

ggplot(out, aes(variable, value)) + geom_boxplot() + facet_wrap(~variable, scales = "free") + geom_point(alpha = 0.5, col = "red", position = position_jitter(width = 0.03)) + geom_hline(data = true_df, aes(yintercept = value), lty = 2, lwd = 1)


# par(mfrow = c(1, 2))
# plot(phi_out, ylim = c(-0.5, 1));abline(h = phi_sim, col = "red")
# segments(1:length(phi_out), l, 1:length(phi_out), u)
# boxplot(phi_out);abline(h = phi_sim, col = "red")

######

# sm <- sampling(stan_gomp,
#   data = list(N = nrow(x), y = log(x$population_untransformed),
#     nu_rate = 0.01,
#     b_lower = -1, b_upper = 2),
#   pars = c("lambda", "sigma_proc", "nu", "b"), iter = iterations,
#   chains = 3, warmup = warmup)
#
# sm_ar1 <- sampling(stan_gomp_ar1,
#   data = list(N = nrow(x), y = log(x$population_untransformed),
#     nu_rate = 0.01,
#     b_lower = -1, b_upper = 2),
#   pars = c("lambda", "sigma_proc", "nu", "b", "phi"), iter = iterations,
#   chains = 3, warmup = warmup)
#
# sm2_ar1 <- sampling(stan_gomp2_ar1,
#   data = list(N = nrow(x), y = log(x$population_untransformed),
#     nu_rate = 0.01,
#     b_lower = -1, b_upper = 2),
#   pars = c("lambda", "sigma_proc","b", "phi"), iter = iterations,
#   chains = 3, warmup = warmup)




# e <- extract(sm)
# p <- lapply(e, median)
#
# lambda <- p$lambda
# sigma_proc <- p$sigma_proc
# nu <- p$nu
# b <- p$b
# phi <- 0
# res <- rep(NA, length(pop))
#
# for(i in 3:length(pop)) {
#   res[i] <- log(pop)[i] - (
#     lambda + b * log(pop[i-1]) +
#       phi * (log(pop[i-1]) - (lambda + b * log(pop[i-2])))
#   )
# }
# (arima(res, order = c(1L, 0L, 0L)))
# par(mfrow = c(2, 1))
# plot(res, type = "l")
# abline(h = 0, lty = 2)
#
# e <- extract(sm_ar1)
# p <- lapply(e, median)
#
# lambda <- p$lambda
# sigma_proc <- p$sigma_proc
# nu <- p$nu
# b <- p$b
# phi <- p$phi
# res <- rep(NA, length(pop))
#
# for(i in 3:length(pop)) {
#   res[i] <- log(pop)[i] - (
#     lambda + b * log(pop[i-1]) +
#       phi * (log(pop[i-1]) - (lambda + b * log(pop[i-2])))
#   )
# }
# (arima(res, order = c(1L, 0L, 0L)))
# plot(res, type = "l")
# abline(h = 0, lty = 2)
