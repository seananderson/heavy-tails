# This file runs the models with an AR1 coefficient on the residuals and low assumed
# observation error

stan_gomp2_ar1_obs <- readRDS("stan-gomp2-ar1-obs.rds")
gpdd <- readRDS("gpdd-clean.rds")

stan_dat <- "list(N = nrow(x), y = log(x$population_untransformed), nu_rate = 0.01, b_lower = -1, b_upper = 2, sigma_obs = 0.1)"

fit_gpdd_model(gpdd_dat = gpdd, model = stan_gomp2_ar1_obs,
  sub_folder = "gomp-ar1-obs-0.1", stan_dat = stan_dat)

out <- plyr::ldply(unique(gpdd$main_id), extract_model,
  subfolder = "gomp-ar1-obs-0.1", get_phi = TRUE, type = "gompertz")

saveRDS(out, file = "gomp_hat_ar1_obs0.1")
