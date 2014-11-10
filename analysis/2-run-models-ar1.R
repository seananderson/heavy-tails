# This file runs the Gompertz models with an AR1 coefficient on the residuals

library("rstan")

# AR1 Gompertz that has been extensively tested
# coded similarly to the Stan 2.4.0 manual
stan_model <-
'data {
  int<lower=0> N; // rows of data
  vector[N] y; // vector to hold observations
  real<lower=0> nu_rate; // rate parameter for nu exponential prior
  real b_lower;
  real b_upper;
}
parameters {
  real lambda;
  real<lower=b_lower, upper=b_upper> b;
  real<lower=0> sigma_proc;
  real<lower=2> nu;
  real<lower=-1, upper=1> phi;
}
transformed parameters {
  vector[N] epsilon;    // error terms
  epsilon[1] <- 0;
  for (i in 2:N) {
    epsilon[i] <- y[i] - (lambda + b * y[i - 1])
                       - (phi * epsilon[i - 1]);
  }
}
model {
  nu ~ exponential(nu_rate);
  lambda ~ normal(0, 10);
  sigma_proc ~ cauchy(0, 2.5);
  phi ~ normal(0, 1);
  for (i in 2:N) {
    y[i] ~ student_t(nu,
                     lambda + b * y[i - 1]
                     + phi * epsilon[i - 1],
                     sigma_proc);
  }
}
'
model <- stan_model(model_code = stan_model)

source("1.5-compile-fit-function.R")
source("1.6-extract-function.R")

gpdd <- readRDS("gpdd-clean.rds")
id <- "gomp-ar1"

fit_gpdd_model(gpdd_dat = gpdd, model = model, sub_folder = id)

out <- plyr::ldply(unique(gpdd$main_id), extract_model,
  sub_folder = id, get_phi = TRUE, type = "gompertz")

saveRDS(out, file = paste0(id, "-hat.rds"))
