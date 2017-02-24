# This file runs the Gompertz models with a more informative prior on nu that
# favours lower values

source("1.5-compile-fit-function.R")
source("1.6-extract-function.R")

model <- readRDS("stan-gomp.rds")
gpdd <- readRDS("gpdd-clean.rds")
id <- "gomp-base-stronger"

stan_pars <- c("lambda", "sigma_proc", "nu", "b")

fit_gpdd_model(gpdd_dat = gpdd, model = model, sub_folder = id,
  pars = stan_pars,
  stan_dat = paste0("list(N = nrow(x), y = log(x$population_untransformed), ",
    "nu_rate = 0.02, b_lower = -1, b_upper = 2)"))

out <- plyr::ldply(unique(gpdd$main_id), extract_model,
  sub_folder = id, get_phi = FALSE, type = "gompertz")

saveRDS(out, file = paste0(id, "-hat.rds"))
