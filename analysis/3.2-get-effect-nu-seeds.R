library(rstan)
stan_t <- readRDS("stan-t.rds")

get_effective_nu_seeds <- function(nu_true = 5, cv = 0.2, N = 50,
  seed_N = 20, start_seed = 0, sigma_proc = 0.5) {

  sigma_proc <- sigma_proc
  nu <- nu_true
  i <- start_seed
  j <- 1
  seeds <- rep(0, length = seed_N)
  nu_ests <- rep(0, length = seed_N)

  while(seeds[length(seeds)] == 0) {
    nu_close <- FALSE
    while(nu_close == FALSE) {
      i <- i + 1
      set.seed(i)
      y <- metRology::rt.scaled(N, df = nu, mean = 0, sd = 0.45)
      sm <- sampling(stan_t,
        data = list(N = N, y = y, nu_rate = 0.01),
        pars = c("sigma_proc", "nu", "mu"), iter = 800,
        chains = 4, warmup = 400)
      med_nu_hat <- median(extract(sm, pars = "nu")[[1]])
      # 88 is a very weird draw with a giant deviation as the 2nd value
      # which creates lambda plots that mess up all the axes
      # so, skipping seed 88 for nu = 5, N = 50
      if(med_nu_hat > (nu - cv) & med_nu_hat < (nu + cv) & i != 88) {
        nu_close <- TRUE
        seeds[j] <- i # a valid seed
        nu_ests[j] <- med_nu_hat # capture the estimated nu
        j <- j + 1
      } else {
        message(paste0("nu hat was ", med_nu_hat))
      }
    }
  }
  return(list(seeds = seeds, nu_ests = nu_ests))
}

nu_3_seeds_N50 <- get_effective_nu_seeds(seed_N = 20, nu_true = 3)
nu_5_seeds_N50 <- get_effective_nu_seeds(seed_N = 20, nu_true = 5)

nu_3_seeds_N100 <- get_effective_nu_seeds(seed_N = 20, nu_true = 3, N = 100)
nu_5_seeds_N100 <- get_effective_nu_seeds(seed_N = 20, nu_true = 5, N = 100)

save(nu_3_seeds_N50, nu_5_seeds_N50, nu_3_seeds_N100, nu_5_seeds_N100,
  file = "nu_effective_seeds.rda")
