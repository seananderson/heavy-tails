stan_logistic_ar1 <- readRDS("stan-logistic-ar1.rds")
# saveRDS(stan_gomp_ar1, file = "stan-gomp-ar1.rds")

# gpdd <- readRDS("gpdd-clean.rds")
library(rstan)
# x <- subset(gpdd, main_id == 20527)
# pop <- x$population_untransformed
# y <- diff(log(x$population_untransformed))
# par(mfrow = c(2, 1))
# plot(log(pop), type = "o")
# plot(y, type = "o")
#stan_gomp_ar1 <- readRDS("stan-gomp-ar1.rds")
#stan_gomp <- readRDS("stan-gomp.rds")
warmup <- 100
iterations <- 200
# z <- seq(-1, 1, length.out = 100);plot(z, dcauchy(z, 0, 1), type = "l", ylim = c(0, 0.4))

###
# OR simulate something
#set.seed(1)
m <- list()

for(j in 1:1) {
set.seed(1)
phi_sim <- 0.4
nu_sim <- 3
sigma_sim <- 0.10
eps_sim <- rep(NA, 50)
eps_sim[1] <- metRology::rt.scaled(1, df = nu_sim, mean = 0, sd = sigma_sim)
  for(i in 2:length(eps_sim)) {
    eps_sim[i] <- eps_sim[i-1] * phi_sim +
      metRology::rt.scaled(1, df = nu_sim, mean = 0, sd = sigma_sim)
  }
#plot(eps_sim)
r_sim <- 0.3
K_sim <- exp(2.5)
y1 <- 2
pop <- rep(NA, length(eps_sim))
pop[1] <- y1
for(i in 2:length(pop)) {
  pop[i] <- pop[i-1] + r_sim * (1 - exp(pop[i-1]) / K_sim) + eps_sim[i]
}
pop <- exp(pop)
par(mfrow = c(1, 1))
plot(log(pop), type = "l")
x <- data.frame(population_untransformed = pop)
#y <- log(x$population_untransformed) / max(log(x$population_untransformed))
y <- log(x$population_untransformed)
K_upper <- max(exp(y)) * 2
#plot(y)
###
###

x$r_obs <- c(NA, diff(log(x$population_untransformed)))
m[[j]] <- sampling(stan_logistic,
  data = list(N = nrow(x) - 1, r_obs = x$r_obs[-1],
    Nt = exp(y[-nrow(x)]) / max(exp(y[-nrow(x)])),
    #Nt = exp(y[-nrow(x)]),
    nu_rate = 0.01, r_upper = 20, K_upper = 2),
  pars = c("r", "sigma_proc","K", "nu"), iter = iterations,
  chains = 4, warmup = warmup)

#m[[j]] <- sampling(stan_logistic_ar1,
  #data = list(N = nrow(x), y = y, nu_rate = 0.01, r_upper = 20, K_upper = K_upper),
  #pars = c("r", "sigma_proc","K", "nu", "phi"), iter = iterations,
  #chains = 2, warmup = warmup)


#   phi_out[j] <- median(extract(sm2_ar1, pars = "phi")[[1]])
#   l[j] <- quantile(extract(sm2_ar1, pars = "phi")[[1]], probs = 0.05)
#   u[j] <- quantile(extract(sm2_ar1, pars = "phi")[[1]], probs = 0.95)
}

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
