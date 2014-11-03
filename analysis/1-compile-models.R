# this file compiles the Stan models, nothing more

library(rstan)

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
}
model {
  nu ~ exponential(nu_rate);
  lambda ~ normal(0, 10);
  sigma_proc ~ cauchy(0, 2.5);
  for (i in 2:N) {
    y[i] ~ student_t(nu, lambda + b * y[i-1], sigma_proc);
  }
}
'

stan_gomp <- stan_model(model_code = stan_model)
saveRDS(stan_gomp, file = "stan-gomp.rds")

# with autocorrelation:

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
model {
  nu ~ exponential(nu_rate);
  lambda ~ normal(0, 10);
  sigma_proc ~ cauchy(0, 2.5);
  phi ~ normal(0, 1);
  for (i in 3:N) {
    y[i] ~ student_t(nu,
            (lambda + b * y[i-1]) + phi * (y[i-1] - (lambda + b * y[i-2])),
            sigma_proc);
  }
}
'
stan_gomp_ar1 <- stan_model(model_code = stan_model)
saveRDS(stan_gomp_ar1, file = "stan-gomp-ar1.rds")

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
stan_gomp_obs <- stan_model(model_code = stan_model)
saveRDS(stan_gomp_obs, file = "stan-gomp-obs.rds")

# a model that just estimates mu, sd, and nu from stationary t distribution:
stan_model <-
'data {
  int<lower=0> N; // rows of data
  vector[N] y; // vector to hold observations
  real<lower=0> nu_rate; // rate parameter for nu exponential prior
}
parameters {
  real<lower=0> sigma_proc;
  real<lower=2> nu;
  real mu;
}
model {
  nu ~ exponential(nu_rate);
  mu ~ normal(0, 10);
  sigma_proc ~ cauchy(0, 2.5);
  y ~ student_t(nu, mu, sigma_proc);
}
'

stan_t <- stan_model(model_code = stan_model)
saveRDS(stan_t, file = "stan-t.rds")

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

stan_logistic <- stan_model(model_code = stan_model)
saveRDS(stan_logistic, file = "stan-logistic.rds")

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
stan_gomp2_ar1 <- stan_model(model_code = stan_model)
saveRDS(stan_gomp2_ar1, file = "stan-gomp2-ar1.rds")

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

stan_rate <- stan_model(model_code = stan_model)
saveRDS(stan_rate, file = "stan-rate.rds")

# uniform prior on sigma:

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
  real<lower=0, upper=20> sigma_proc;
  real<lower=2> nu;
}
model {
  nu ~ exponential(nu_rate);
  lambda ~ normal(0, 10);
  for (i in 2:N) {
    y[i] ~ student_t(nu, lambda + b * y[i-1], sigma_proc);
  }
}
'

stan_gomp_uniform <- stan_model(model_code = stan_model)
saveRDS(stan_gomp_uniform, file = "stan-gomp-uniform.rds")
