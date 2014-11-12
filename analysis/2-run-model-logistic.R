# This file runs the models with a logistic growth model

source("1.5-compile-fit-function.R")
source("1.6-extract-function.R")

model <- readRDS("stan-logistic.rds")
gpdd <- readRDS("gpdd-clean.rds")
id <- "logistic"

stan_dat <- "list(N = nrow(x) - 1, r_obs = x$r_obs[-1], nu_rate = 0.01, K_upper = 2, r_upper = 10, Nt = x$population_untrans_scaled[-N])"
stan_pars <- c("r", "K", "nu", "sigma_proc")

fit_gpdd_model(gpdd_dat = gpdd, model = model,
  sub_folder = id, stan_dat = stan_dat, pars = stan_pars)

out <- plyr::ldply(unique(gpdd$main_id), extract_model,
  sub_folder = id, get_phi = FALSE, type = "logistic")

saveRDS(out, file = paste0(id, "-hat.rds"))
