gpdd <- readRDS("gpdd-clean.rds")
library(rstan)
x <- subset(gpdd, main_id == 20527)
pop <- x$population_untransformed
y <- diff(log(x$population_untransformed))
par(mfrow = c(2, 1))
plot(log(pop), type = "o")
plot(y, type = "o")
stan_gomp_ar1 <- readRDS("stan-gomp-ar1.rds")
stan_gomp <- readRDS("stan-gomp.rds")
warmup <- 400
iterations <- 800


###
# OR simulate something
phi_sim <- 0.9
eps_sim <- arima.sim(n = 100, list(ar = phi_sim, 1, 1))
lamb_sim <- 10
b_sim <- 0
y1 <- 1
pop <- rep(NA, 100)
pop[1] <- y1
for(i in 2:100) {
  pop[i] <- lamb_sim + b_sim * pop[i-1] + eps_sim[i]
}
pop <- exp(pop)
par(mfrow = c(1, 1))
plot(log(pop), type = "l")
x <- data.frame(population_untransformed = pop)
###

sm <- sampling(stan_gomp,
  data = list(N = nrow(x), y = log(x$population_untransformed),
    nu_rate = 0.01,
    b_lower = -1, b_upper = 2),
  pars = c("lambda", "sigma_proc", "nu", "b"), iter = iterations,
  chains = 3, warmup = warmup)

sm_ar1 <- sampling(stan_gomp_ar1,
  data = list(N = nrow(x), y = log(x$population_untransformed),
    nu_rate = 0.01,
    b_lower = -1, b_upper = 2),
  pars = c("lambda", "sigma_proc", "nu", "b", "phi"), iter = iterations,
  chains = 3, warmup = warmup)

e <- extract(sm)
p <- lapply(e, median)

lambda <- p$lambda
sigma_proc <- p$sigma_proc
nu <- p$nu
b <- p$b
phi <- 0
res <- rep(NA, length(pop))

for(i in 3:length(pop)) {
  res[i] <- log(pop)[i] - (
    lambda + b * log(pop[i-1]) +
    phi * (log(pop[i-1]) - (lambda + b * log(pop[i-2])))
    )
}
(arima(res, order = c(1L, 0L, 0L)))
par(mfrow = c(2, 1))
plot(res, type = "l")
abline(h = 0, lty = 2)

e <- extract(sm_ar1)
p <- lapply(e, median)

lambda <- p$lambda
sigma_proc <- p$sigma_proc
nu <- p$nu
b <- p$b
phi <- p$phi
res <- rep(NA, length(pop))

for(i in 3:length(pop)) {
  res[i] <- log(pop)[i] - (
    lambda + b * log(pop[i-1]) +
    phi * (log(pop[i-1]) - (lambda + b * log(pop[i-2])))
    )
}
(arima(res, order = c(1L, 0L, 0L)))
plot(res, type = "l")
abline(h = 0, lty = 2)
