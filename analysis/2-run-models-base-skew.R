# This file runs the base Gompertz models (with skewed t process errors)

source("1.5-compile-fit-function.R")
source("1.6-extract-function.R")

model <- readRDS("stan-gomp-skew.rds")
gpdd <- readRDS("gpdd-clean.rds")
id <- "gomp-base-skew"

stan_pars <- c("lambda", "sigma_proc", "nu", "b", "log_skew")

library("rstan")
options(mc.cores = 4L)
fit_gpdd_model(gpdd_dat = gpdd, model = model, sub_folder = id,
  pars = stan_pars, iterations = 20000, max_iterations = 20000, warmup = 10000)

out <- plyr::ldply(unique(gpdd$main_id), extract_model,
  sub_folder = id, get_phi = FALSE, type = "gompertz",
  get_nu = TRUE, get_skew = TRUE)

saveRDS(out, file = paste0(id, "-hat.rds"))
