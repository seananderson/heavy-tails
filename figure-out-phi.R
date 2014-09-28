stan_model <-
  'data {
  int<lower=0> N; // rows of data
  vector[N] y; // vector to hold observations
  real<lower=0> nu_rate; // rate parameter for nu exponential prior
  real b_lower;
  real b_upper;
}
parameters {
  real lambda;
  real<lower=b_lower, upper=b_upper> b;
  real<lower=0> sigma_proc;
  real<lower=2> nu;
  real<lower=-1, upper=1> phi;
}
transformed parameters {
  vector[N] epsilon;    // error terms
  epsilon[1] <- 0;
  for (i in 2:N) {
       epsilon[i] <- y[i] - (lambda + b * y[i - 1])
                          - (phi * epsilon[i - 1]);
  }
}
model {
  nu ~ exponential(nu_rate);
  lambda ~ normal(0, 10);
  sigma_proc ~ cauchy(0, 2.5);
  phi ~ normal(0, 1);
  for (i in 2:N) {
    y[i] ~ student_t(nu,
                     lambda + b * y[i - 1]
                     + phi * epsilon[i - 1],
                     sigma_proc);
  }
}
'
stan_gomp2_ar1 <- stan_model(model_code = stan_model)
# saveRDS(stan_gomp_ar1, file = "stan-gomp-ar1.rds")

# gpdd <- readRDS("gpdd-clean.rds")
library(rstan)
# x <- subset(gpdd, main_id == 20527)
# pop <- x$population_untransformed
# y <- diff(log(x$population_untransformed))
# par(mfrow = c(2, 1))
# plot(log(pop), type = "o")
# plot(y, type = "o")
stan_gomp_ar1 <- readRDS("stan-gomp-ar1.rds")
stan_gomp <- readRDS("stan-gomp.rds")
warmup <- 200
iterations <- 400
# z <- seq(-1, 1, length.out = 100);plot(z, dcauchy(z, 0, 1), type = "l", ylim = c(0, 0.4))

###
# OR simulate something
#set.seed(1)
m <- list()

load("nu_effective_seeds.rda")

for(j in 1:20) {
  set.seed(nu_5_seeds_N50$seeds[j])
phi_sim <- 0.2
nu_sim <- 100
sigma_sim <- 0.65
eps_sim <- rep(NA, 50)
eps_sim[1] <- metRology::rt.scaled(1, df = nu_sim, mean = 0, sd = sigma_sim)
  for(i in 2:length(eps_sim)) {
    eps_sim[i] <- eps_sim[i-1] * phi_sim +
      metRology::rt.scaled(1, df = nu_sim, mean = 0, sd = sigma_sim)
  }
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

m[[j]] <- sampling(stan_gomp2_ar1,
  data = list(N = nrow(x), y = log(x$population_untransformed),
    nu_rate = 0.01,
    b_lower = -1, b_upper = 2),
  pars = c("lambda", "sigma_proc","b", "phi", "nu"), iter = iterations,
  chains = 3, warmup = warmup)


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
