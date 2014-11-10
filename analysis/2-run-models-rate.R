# This file runs the models with a logistic growth model

library("rstan")
# and a basic model of just the growth rates
# same as setting b = 1 in Gompertz

stan_model <-
'data {
  int<lower=0> N;
  vector[N] r_obs;
  real<lower=0> nu_rate;
}
parameters {
  real lambda;
  real<lower=0> sigma_proc;
  real<lower=2> nu;
}
model {
  nu ~ exponential(nu_rate);
  lambda ~ normal(0, 10);
  sigma_proc ~ cauchy(0, 2.5);
  r_obs ~ student_t(nu, lambda, sigma_proc);
}
'
model <- stan_model(model_code = stan_model)
saveRDS(stan_rate, file = "stan-rate.rds")

source("1.5-compile-fit-function.R")
source("1.6-extract-function.R")

gpdd <- readRDS("gpdd-clean.rds")
id <- "rate"

stan_dat <- "list(N = nrow(x) - 1, r_obs = x$r_obs[-1], nu_rate = 0.01)"
stan_pars <- c("lambda", "nu", "sigma_proc")

fit_gpdd_model(gpdd_dat = gpdd, model = model,
  sub_folder = id, stan_dat = stan_dat, pars = stan_pars,
  iterations = 2000, warmup = 1000)

out <- plyr::ldply(unique(gpdd$main_id), extract_model,
  sub_folder = id, get_phi = FALSE, type = "rate")

saveRDS(out, file = paste0(id, "-hat.rds"))
