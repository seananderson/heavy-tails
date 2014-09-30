# This file runs the Gompertz models with an AR1 coefficient on the residuals
# and low assumed observation error

source("1.5-compile-fit-function.R")
source("1.6-extract-function.R")

model <- readRDS("stan-gomp-obs-equal.rds")
gpdd <- readRDS("gpdd-clean.rds")
id <- "gomp-obs-equal"

stan_dat <- "list(N = nrow(x), y = log(x$population_untransformed), nu_rate = 0.01, b_lower = -1, b_upper = 2)"

stan_pars <- c("lambda", "sigma_proc", "nu", "b")

fit_gpdd_model(gpdd_dat = gpdd, model = model,
  sub_folder = id, stan_dat = stan_dat, iterations = 2000, warmup = 1000,
  max_iterations = 4000, pars = stan_pars)

out <- plyr::ldply(unique(gpdd$main_id), extract_model,
  sub_folder = id, get_phi = FALSE, type = "gompertz")

saveRDS(out, file = paste0(id, "-hat.rds"))
