# This file runs the Gompertz models with an AR1 coefficient on the residuals

source("1.5-compile-fit-function.R")
source("1.6-extract-function.R")

options(mc.cores = min(c(4L, parallel::detectCores())))

model <- readRDS("stan-gomp2-ar1.rds")
gpdd <- readRDS("gpdd-clean.rds")
id <- "gomp-ar1"

stan_dat <- "list(N = nrow(x), y = log(x$population_untransformed), nu_rate = 0.01, b_lower = -1, b_upper = 2, phi_sd = 0.5)"

fit_gpdd_model(gpdd_dat = gpdd, model = model, sub_folder = id,
  iterations = 2000, max_iterations = 128000*2, warmup = 1000, chains = 4, root_folder = ".",
  stan_dat = stan_dat)

out <- plyr::ldply(unique(gpdd$main_id), extract_model,
  sub_folder = id, get_phi = TRUE, type = "gompertz", root_folder = ".")

library(dplyr)
ids_non_converged <- filter(out,  max_rhat > 1.05 | min_neff < 200)$main_id
ids_non_converged

saveRDS(out, file = paste0(id, "-hat.rds"))

