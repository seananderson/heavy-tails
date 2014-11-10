# This file runs the models with a logistic growth model

# logistic-growth model, parameterized as in de Valpine and
# Hastings 2002 Ecol. Mono. or Clark et al. 2010, Meth. Ecol. Evol.
stan_model <-
'data {
  int<lower=0> N; // rows of data
  vector[N] r_obs; // observed growth rates for time t
  vector[N] Nt; // numbers at time t
  real<lower=0> nu_rate; // rate parameter for nu exponential prior
  real<lower=0> K_upper;
  real<lower=0> r_upper;
}
parameters {
  real lambda;
  real<lower=0, upper=K_upper> K;
  real<lower=0, upper=r_upper> r;
  real<lower=0> sigma_proc;
  real<lower=2> nu;
}
model {
  nu ~ exponential(nu_rate);
  sigma_proc ~ cauchy(0, 2.5);
  r_obs ~ student_t(nu, r * (1 - (Nt/K)), sigma_proc);
}
'
model <- stan_model(model_code = stan_model)

source("1.5-compile-fit-function.R")
source("1.6-extract-function.R")

gpdd <- readRDS("gpdd-clean.rds")
id <- "logistic"

stan_dat <- "list(N = nrow(x) - 1, r_obs = x$r_obs[-1], nu_rate = 0.01, K_upper = 2, r_upper = 10, Nt = x$population_untrans_scaled[-N])"
stan_pars <- c("r", "K", "nu", "sigma_proc")

fit_gpdd_model(gpdd_dat = gpdd, model = model,
  sub_folder = id, stan_dat = stan_dat, pars = stan_pars)

out <- plyr::ldply(unique(gpdd$main_id), extract_model,
  sub_folder = id, get_phi = FALSE, type = "logistic")

saveRDS(out, file = paste0(id, "-hat.rds"))
