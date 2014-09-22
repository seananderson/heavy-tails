# this file runs the models *without* AR1 coefficient on the residuals

stan_gomp <- readRDS("stan-gomp.rds")
gpdd <- readRDS("gpdd-clean.rds")
library(rstan)

main_id_vec <- unique(gpdd$main_id)

out_gomp <- plyr::d_ply(subset(gpdd, main_id %in% main_id_vec), "main_id",
  function(x) {

    max_rhat <- 999
    min_neff <- 0
    iterations<- 1000
    warmup <- 500
    file_base <- "/global/scratch/anderson/heavy/base-5.0/sm-"

    if(!file.exists(paste0(file_base, unique(x$main_id), ".rds"))) {
      print("already done")
      print(unique(x$main_id))
      while((max_rhat > 1.05 | min_neff < 200) & iterations <= 4000) {

        sm <- sampling(stan_gomp,
          data = list(N = nrow(x), y = log(x$population_untransformed),
            nu_rate = 0.01,
            b_lower = -1, b_upper = 2),
          pars = c("lambda", "sigma_proc", "nu", "b"), iter = iterations,
          chains = 4, warmup = warmup)

        # check:
        max_rhat <- max(summary(sm)$summary[, "Rhat"])
        min_neff <- min(summary(sm)$summary[, "n_eff"])

        warmup <- warmup * 2
        iterations <- iterations * 2
      }
      saveRDS(sm, file = paste0(file_base, unique(x$main_id), ".rds"))
      sink(paste0(file_base, unique(x$main_id), ".txt"))
      print(sm)
      sink()
    }
  })

# trash ones we no longer want:
#f <- list.files("/global/scratch/anderson/heavy/ar1-4.0", pattern = "*\\.rds")
#fn <- sub("sm-([0-9]+)\\.rds", "\\1", f)
#fn[!fn %in% gpdd$main_id]
#file.rename(

