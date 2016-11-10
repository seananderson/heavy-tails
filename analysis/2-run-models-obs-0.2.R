# This file runs the Gompertz models with an AR1 coefficient on the residuals
# and low assumed observation error

source("1.5-compile-fit-function.R")
source("1.6-extract-function.R")

model <- readRDS("stan-gomp-obs.rds")
gpdd <- readRDS("gpdd-clean.rds")
id <- "gomp-obs-0.2"

stan_dat <- "list(N = nrow(x), y = log(x$population_untransformed), nu_rate = 0.01, b_lower = -1, b_upper = 2, sigma_obs = 0.2)"

stan_pars <- c("lambda", "sigma_proc", "nu", "b")

library("rstan")
options(mc.cores = 4L)

fit_gpdd_model(gpdd_dat = gpdd, model = model,
 sub_folder = id, stan_dat = stan_dat, iterations = 2000, warmup = 1000,
 max_iterations = 4000, pars = stan_pars)
out <- plyr::ldply(unique(gpdd$main_id), extract_model,
 sub_folder = id, get_phi = FALSE, type = "gompertz")
saveRDS(out, file = paste0(id, "-hat.rds"))

out <- readRDS(paste0(id, "-hat.rds"))
nc <- subset(out, max_rhat > 1.05 | min_neff < 200)$main_id

id <- "gomp-obs-0.2-2"
# didn't converge:
fit_gpdd_model(gpdd_dat = subset(gpdd, main_id %in% nc), model = model,
  sub_folder = id, stan_dat = stan_dat, iterations = 8000, warmup = 8000/2,
  max_iterations = 256000, pars = stan_pars, thin = 10)

out2 <- plyr::ldply(nc, extract_model,
  sub_folder = id, get_phi = FALSE, type = "gompertz", root_folder = ".")
saveRDS(out2, file = paste0(id, "-hat.rds"))

library(dplyr)
out <- readRDS("gomp-obs-0.2-hat.rds")
out2 <- readRDS("gomp-obs-0.2-2-hat.rds")
out <- out[!out$main_id %in% out2$main_id, ]
out <- bind_rows(out, out2) %>% as.data.frame()

saveRDS(out, file = "gomp-obs-0.2-extra-iterations-hat.rds")
