# This file runs the models with a logistic growth model

source("1.5-compile-fit-function.R")
source("1.6-extract-function.R")

model <- readRDS("stan-rate.rds")
gpdd <- readRDS("gpdd-clean.rds")
id <- "rate"

stan_dat <- "list(N = nrow(x) - 1, r_obs = x$r_obs[-1], nu_rate = 0.01)"
stan_pars <- c("lambda", "nu", "sigma_proc")

fit_gpdd_model(gpdd_dat = gpdd, model = model,
  sub_folder = id, stan_dat = stan_dat, pars = stan_pars)

out <- plyr::ldply(unique(gpdd$main_id), extract_model,
  sub_folder = id, get_phi = FALSE, type = "rate")

saveRDS(out, file = paste0(id, "-hat.rds"))
