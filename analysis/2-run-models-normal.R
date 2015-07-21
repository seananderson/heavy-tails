# This file runs the base Gompertz models with normally distributed process noise

source("1.5-compile-fit-function.R")
source("1.6-extract-function.R")

model <- readRDS("stan-gomp-normal.rds")
gpdd <- readRDS("gpdd-clean.rds")
id <- "gomp-base-normal"

stan_pars <- c("lambda", "sigma_proc", "b")

fit_gpdd_model(gpdd_dat = gpdd, model = model, sub_folder = id,
  pars = stan_pars, iterations = 20000, max_iterations = 20000, warmup = 10000)

out <- plyr::ldply(unique(gpdd$main_id), extract_model,
  sub_folder = id, get_phi = FALSE, type = "gompertz", get_nu = FALSE)

saveRDS(out, file = paste0(id, "-hat.rds"))
