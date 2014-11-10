# This file runs the Gompertz models with an AR1 coefficient on the residuals
# and low assumed observation error

library("rstan")
# with fixed observation error (and no autocorrelation):
stan_model <-
  'data {
  int<lower=0> N; // rows of data
  vector[N] y; // vector to hold observations
  real<lower=0> nu_rate; // rate parameter for nu exponential prior
  real b_lower;
  real b_upper;
  real<lower=0> sigma_obs;
}
parameters {
  real lambda;
  real<lower=b_lower, upper=b_upper> b;
  real<lower=0> sigma_proc;
  real<lower=2> nu;
  vector[N] U; // states
}
model {
  nu ~ exponential(nu_rate);
  lambda ~ normal(0, 10);
  sigma_proc ~ cauchy(0, 2.5);
  for (i in 2:N) {
    U[i] ~ student_t(nu, lambda + b * U[i-1], sigma_proc);
  }
  y ~ normal(U, sigma_obs);
}
'
model <- stan_model(model_code = stan_model)

source("1.5-compile-fit-function.R")
source("1.6-extract-function.R")

gpdd <- readRDS("gpdd-clean.rds")
id <- "gomp-obs-0.2"

stan_dat <- "list(N = nrow(x), y = log(x$population_untransformed), nu_rate = 0.01, b_lower = -1, b_upper = 2, sigma_obs = 0.2)"

stan_pars <- c("lambda", "sigma_proc", "nu", "b")

fit_gpdd_model(gpdd_dat = gpdd, model = model,
  sub_folder = id, stan_dat = stan_dat, iterations = 2000, warmup = 1000,
  max_iterations = 4000, pars = stan_pars)

out <- plyr::ldply(unique(gpdd$main_id), extract_model,
  sub_folder = id, get_phi = FALSE, type = "gompertz")

saveRDS(out, file = paste0(id, "-hat.rds"))
