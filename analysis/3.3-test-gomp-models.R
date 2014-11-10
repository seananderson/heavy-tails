# this file runs some simulation tests of the Stan model
# I have set the parameter values to ballpark median values from the GPDD
# this file tests the ability for the Gompertz to perform well
# when we feed it process deviations that have been selected to have the true
# population value of nu (within a certain CV)
#
# Note that phi is being ignored now

library("rstan")
stan_gomp_obs <- readRDS("stan-gomp-obs.rds")
stan_gomp <- readRDS("stan-gomp.rds")

sim_gomp <- function(id = 99, nu = 10, sigma_proc = 0.65,
  sigma_obs_true = 0.001, sigma_obs_assumed = 0.001, N = 50, b = 0.50,
  lambda = 1.5, y1 = 3, iterations = 2000, warmup = 1000, plot = FALSE,
  max_iter = 8000, chains = 4, seed) {

  file_base <- paste0("/global/scratch/anderson/heavy/sim4/sm-",
    id, "-sigma_obs_true-", sigma_obs_true, "-sigma_obs_assumed-",
    sigma_obs_assumed, "-N-", N, "-b-", b, "-sigma_proc-", sigma_proc,
    "-nu-", nu, "-y1-", y1, "-lambda-", lambda)

  # simulate data:
  y <- vector(length = N)
  set.seed(seed)
  proc_error <- metRology::rt.scaled(N, df = nu, mean = 0, sd = sigma_proc)
  for(i in 2:N) {
    y[i] <- lambda + b * y[i-1] + proc_error[i-1]
  }
  y <- rnorm(N, mean = y, sd = sigma_obs_true)

  if(plot == TRUE) plot(y)

  # fit model:
  min_neff <- 0
  max_rhat <- 999

  if(!file.exists(paste0(file_base, ".rds"))) {
    while((max_rhat > 1.05 | min_neff < 200) & iterations <= max_iter) {
      if(sigma_obs_assumed <= 0.01) { # use the faster Stan model without obs error
        sm <- sampling(stan_gomp,
          data = list(N = N, y = y, nu_rate = 0.01,
            b_lower = -1, b_upper = 2),
          pars = c("lambda", "sigma_proc", "nu", "b"), iter = iterations,
          chains = chains, warmup = warmup)
      } else {
        sm <- sampling(stan_gomp_obs,
          data = list(N = N, y = y, nu_rate = 0.01,
            b_lower = -1, b_upper = 2, sigma_obs = sigma_obs_assumed),
          pars = c("lambda", "sigma_proc", "nu", "b"), iter = iterations,
          chains = chains, warmup = warmup)
      }
      max_rhat <- max(summary(sm)$summary[, "Rhat"])
      min_neff <- min(summary(sm)$summary[, "n_eff"])
      warmup <- warmup * 2
      iterations <- iterations * 2
    }
    if(file.exists("/global/scratch/anderson/heavy/")) {
      saveRDS(sm, file = paste0(file_base, ".rds"))
      sink(file = paste0(file_base, ".txt"))
      print(sm)
      sink()
    }
  } else { # already run this one
    message("reloading past model run")
    sm <- readRDS(paste0(file_base, ".rds"))
  }
  nu_samples <- extract(sm, pars = "nu")[[1]]
  p_0.1 <- length(nu_samples[nu_samples < 10]) / length(nu_samples)
  p_0.2 <- length(nu_samples[nu_samples < 20]) / length(nu_samples)
  m <- as.data.frame(lapply(extract(sm), quantile, probs = 0.5))
  l <- as.data.frame(lapply(extract(sm), quantile, probs = 0.1))
  u <- as.data.frame(lapply(extract(sm), quantile, probs = 0.9))
  names(l) <- paste0(names(l), "_l")
  names(u) <- paste0(names(u), "_u")
  out <- data.frame(id = id, m, l, u, p_0.1 = p_0.1, p_0.2 = p_0.2,
    max_rhat = max_rhat, min_neff = min_neff, N = N,
    b_true = b, lambda_true = lambda,
    nu_true = nu, sigma_proc_true = sigma_proc,
    sigma_obs_true = sigma_obs_true, sigma_obs_assumed = sigma_obs_assumed)
  out$nu_true <- nu
  out
}

iters <- 1:20

load("nu_effective_seeds.rda")
nu_norm_seeds <- 1:200

check_nu_1 <- plyr::ldply(iters, function(i) sim_gomp(i, nu = 3,
    seed = nu_3_seeds_N50$seeds[i]))
check_nu_2 <- plyr::ldply(iters, function(i) sim_gomp(i, nu = 5,
    seed = nu_5_seeds_N50$seeds[i]))
check_nu_4 <- plyr::ldply(iters, function(i) sim_gomp(i, nu = 1e9,
    seed = nu_norm_seeds[i]))
check_nu_base <- rbind(check_nu_1, check_nu_2, check_nu_4)

# now with 100 data points:
check_nu_1 <- plyr::ldply(iters, function(i) sim_gomp(i, nu = 3,
    seed = nu_3_seeds_N100$seeds[i], N = 100))
check_nu_2 <- plyr::ldply(iters, function(i) sim_gomp(i, nu = 5,
    seed = nu_5_seeds_N100$seeds[i], N = 100))
check_nu_4 <- plyr::ldply(iters, function(i) sim_gomp(i, nu = 1e9,
    seed = nu_norm_seeds[i], N = 100))
check_nu_N100 <- rbind(check_nu_1, check_nu_2, check_nu_4)

# and with 50 data points and 0.3 observation error (but not estimated) to see
# how this obscures the signal
check_nu_1 <- plyr::ldply(iters, function(i) sim_gomp(i, nu = 3,
    seed = nu_3_seeds_N50$seeds[i], sigma_obs_true = 0.2,
  sigma_obs_assumed = 0.001))
check_nu_2 <- plyr::ldply(iters, function(i) sim_gomp(i, nu = 5,
    seed = nu_5_seeds_N50$seeds[i], sigma_obs_true = 0.2,
  sigma_obs_assumed = 0.001))
check_nu_4 <- plyr::ldply(iters, function(i) sim_gomp(i, nu = 1e9,
    seed = nu_norm_seeds[i], sigma_obs_true = 0.2,
  sigma_obs_assumed = 0.001))
check_nu_sigma_obs_ignored <- rbind(check_nu_1, check_nu_2, check_nu_4)

# and try with 0.3 sigma_obs that is assumed known
check_nu_1 <- plyr::ldply(iters, function(i) sim_gomp(i, nu = 3,
    seed = nu_3_seeds_N50$seeds[i], sigma_obs_true = 0.2, sigma_obs_assumed = 0.2, iterations = 4000, warmup = 2000))
check_nu_2 <- plyr::ldply(iters, function(i) sim_gomp(i, nu = 5,
    seed = nu_5_seeds_N50$seeds[i], sigma_obs_true = 0.2, sigma_obs_assumed = 0.2, iterations = 4000, warmup = 2000))
check_nu_4 <- plyr::ldply(iters, function(i) sim_gomp(i, nu = 1e9,
    seed = nu_norm_seeds[i], sigma_obs_true = 0.2, sigma_obs_assumed = 0.2, iterations = 4000, warmup = 2000))
check_nu_sigma_obs_assumed <- rbind(check_nu_1, check_nu_2, check_nu_4)

check_nu <- rbind(check_nu_base, check_nu_N100, check_nu_sigma_obs_ignored,
 check_nu_sigma_obs_assumed)

saveRDS(check_nu, file = "check_nu.rds")
