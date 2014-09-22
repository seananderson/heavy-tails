# this file runs the Gompertz model without autocorrelation
# I'm calling this the 'base' run throughout

stan_gomp <- readRDS("stan-gomp.rds")
gpdd <- readRDS("gpdd-clean.rds")
library(rstan)

main_id_vec <- unique(gpdd$main_id)

out_gomp <- plyr::d_ply(subset(gpdd, main_id %in% main_id_vec), "main_id",
  function(x) {
    max_rhat <- 999
    min_neff <- 0
    iterations <- 500
    warmup <- 250

    while((max_rhat > 1.1 | min_neff < 100) & iterations <= 4000) {

      sm <- sampling(stan_gomp,
        data = list(N = nrow(x), y = log(x$population_untransformed),
          nu_rate = 0.01,
          b_lower = -1, b_upper = 2),
        pars = c("lambda", "sigma_proc", "nu", "b"), iter = iterations,
        chains = 4, warmup = warmup)

      # check:
      max_rhat <- max(summary(sm)$summary[, "Rhat"])
      min_neff <- min(summary(sm)$summary[, "n_eff"])

      if((max_rhat > 1.1 | min_neff < 100)) {
        warmup <- warmup * 2
        iterations <- iterations * 2
      }
    }

    if(file.exists("/global/scratch/anderson/heavy/")) {
      saveRDS(sm, file = paste0("/global/scratch/anderson/heavy/base/sm-",
        unique(x$main_id), ".rds"))
      sink(paste0("/global/scratch/anderson/heavy/base/sm-",
        unique(x$main_id), ".txt"))
      print(sm)
      sink()
    }

  })
