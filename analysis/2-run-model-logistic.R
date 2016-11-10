# This file runs the models with a logistic growth model

source("1.5-compile-fit-function.R")
source("1.6-extract-function.R")

model <- readRDS("stan-logistic2.rds")
gpdd <- readRDS("gpdd-clean.rds")
id <- "logistic"

stan_dat <- paste("list(N = nrow(x), y = log(x$population_untransformed),",
  "nu_rate = 0.01, K_upper = 2 * max(x$population_untransformed), r_upper = 10)")
stan_pars <- c("r", "K", "nu", "sigma_proc")

library("rstan")
options(mc.cores = 1)
fit_gpdd_model(gpdd_dat = gpdd, model = model,
  sub_folder = id, stan_dat = stan_dat, pars = stan_pars)

out <- plyr::ldply(unique(gpdd$main_id), extract_model,
  sub_folder = id, get_phi = FALSE, type = "logistic", root_folder = ".")

saveRDS(out, file = paste0(id, "-hat.rds"))
