fit_gpdd_model <- function(gpdd_dat, model, sub_folder,
  root_folder = "~/scratch/heavy", file_prefix = "sm",
  stan_dat = paste0("list(N = nrow(x), y = log(x$population_untransformed), ",
    "nu_rate = 0.01, b_lower = -1, b_upper = 2)"),
  pars = c("lambda", "sigma_proc", "nu", "b", "phi"), max_rhat_allowed = 1.05,
  min_neff_allowed = 200, iterations = 2000, max_iterations = 8000,
  iteration_increment = 2, warmup = 1000, chains = 4, overwrite = FALSE,
  .parallel = TRUE, refresh = -1) {

  library(rstan)

  if(.parallel) {
    library(doParallel)
    library(foreach)
    registerDoParallel(cores = 2)
    # getDoParWorkers() # check
  }

  if(!file.exists(paste0(root_folder, "/", sub_folder)))
    dir.create(paste0(root_folder, "/", sub_folder), recursive = TRUE)

  max_rhat <- 999
  min_neff <- 0
  file_base <- paste0(root_folder, "/", sub_folder, "/", file_prefix, "-")

  plyr::d_ply(gpdd_dat, "main_id", function(x) {
    this_file <- paste0(file_base, unique(x$main_id))

    # some models (Ricker-logistic) may model growth rates as the response
    # or use data scaled to the maximum for computational reasons:
    x$r_obs <- c(NA, diff(log(x$population_untransformed)))
    x$population_untrans_scaled <-
       x$population_untransformed / max(x$population_untransformed)
    N <- nrow(x)

    if(!file.exists(paste0(this_file, ".rds")) | overwrite) {
      while((max_rhat > max_rhat_allowed | min_neff < min_neff_allowed) &
          iterations <= max_iterations) {
        d <- eval(parse(text = stan_dat))
        # refresh = -1 supresses the running messages
        sm <- sampling(model, data = d, refresh = refresh,
          pars = pars, iter = iterations, chains = chains, warmup = warmup)
        max_rhat <- max(summary(sm)$summary[, "Rhat"])
        min_neff <- min(summary(sm)$summary[, "n_eff"])
        warmup <- warmup * iteration_increment
        iterations <- iterations * iteration_increment
      }
      saveRDS(sm, file = paste0(this_file, ".rds"))
      sink(paste0(this_file, ".txt"))
      print(sm)
      sink()
    }
  }, .parallel = .parallel)
}
