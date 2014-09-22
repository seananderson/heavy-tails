# this file runs the Gompertz model with AR1 residual coefficient
# but also adds an assumed quantity of observation error
# it does not estimate the observation error

stan_gomp_ar1_obs <- readRDS("stan-gomp-ar1-obs.rds")
gpdd <- readRDS("gpdd-clean.rds")
library(rstan)

# Mammals only:
main_id_vec <- unique(subset(gpdd, taxonomic_class == "Mammalia")$main_id)
#main_id_vec <- sample(main_id_vec, 4)

out_gomp <- plyr::d_ply(subset(gpdd, main_id %in% main_id_vec), "main_id",
  function(x) {
    max_rhat <- 999
    min_neff <- 0
    iterations <- 2000
    warmup <- 1000
    if(!file.exists(paste0("/global/scratch/anderson/heavy/ar1-obs0.2/sm-",
      unique(x$main_id), ".rds"))) {
      while((max_rhat > 1.1 | min_neff < 100) & iterations <= 4000) {
        sm <- sampling(stan_gomp_ar1_obs,
          data = list(N = nrow(x), y = log(x$population_untransformed),
            nu_rate = 0.01,
            b_lower = -1, b_upper = 2, sigma_obs = 0.2),
          pars = c("lambda", "sigma_proc", "nu", "b", "phi"), iter = iterations,
          chains = 4, warmup = warmup)

        # check:
        max_rhat <- max(summary(sm)$summary[, "Rhat"])
        min_neff <- min(summary(sm)$summary[, "n_eff"])

        warmup <- warmup * 2
        iterations <- iterations * 2
      }

      if(file.exists("/global/scratch/anderson/heavy/ar1-obs0.2/")) {
        saveRDS(sm, file = paste0("/global/scratch/anderson/heavy/ar1-obs0.2/sm-",
          unique(x$main_id), ".rds"))
        sink(paste0("/global/scratch/anderson/heavy/ar1-obs0.2/sm-",
          unique(x$main_id), ".txt"))
        print(sm)
        sink()
      }
    }
  })
