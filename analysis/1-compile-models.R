# this file compiles the Stan models, nothing more

library(rstan)

# the main heavy-tailed model:

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

# the main model but with the Gamma prior:

stan_model <-
'data {
  int<lower=0> N; // rows of data
  vector[N] y; // vector to hold observations
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
  nu ~ gamma(2,0.1); // prior from https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
  lambda ~ normal(0, 10);
  sigma_proc ~ cauchy(0, 2.5);
  for (i in 2:N) {
    y[i] ~ student_t(nu, lambda + b * y[i-1], sigma_proc);
  }
}
'

stan_gomp_gamma <- stan_model(model_code = stan_model)
saveRDS(stan_gomp_gamma, file = "stan-gomp-gamma.rds")

# same base model but with skewed heavy-tailed process noise:

stan_model <-
'
functions {
  real skew_student_t_log(real y, real nu, real mu, real sigma, real skew) {
  real lp;
  if (skew <= 0)
    reject("Skew has to be positive. Found skew=", skew);
  if (sigma <= 0)
    reject("Scale has to be positive.  Found sigma=", sigma);
  lp <- log(skew) - log1p(square(skew));
  if (y < mu)
    return lp + student_t_log(y * skew, nu, mu * skew, sigma);
  else
    return lp + student_t_log(y / skew, nu, mu / skew, sigma);
  }
}
data {
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
  real log_skew;
}
model {
  nu ~ exponential(nu_rate);
  lambda ~ normal(0, 10);
  sigma_proc ~ cauchy(0, 2.5);
  log_skew ~ cauchy(0, 2.5);
  for (i in 2:N) {
    y[i] ~ skew_student_t(nu, lambda + b * y[i-1], sigma_proc, exp(log_skew));
  }
}
'
stan_gomp_skew <- stan_model(model_code = stan_model)
saveRDS(stan_gomp_skew, file = "stan-gomp-skew.rds")

# same base model but with skewed heavy-tailed process noise:
# (skew t as in Azzalini and Arellano-Valle and the sn R package)
# (quite inefficient to sample from)
stan_model <-
'
functions {
  real skew_student_t_log(real y, real nu, real mu, real sigma, real skew) {
   real z; real zc;
   if (sigma <= 0)
     reject("Scale has to be positive. Found sigma=", sigma);
   z <- (y - mu) / sigma;
   zc <- skew * z * sqrt((nu+1) / (nu + square(z)));
   return -log(sigma) + student_t_log(z, nu, 0, 1) + student_t_cdf_log(zc, nu + 1, 0, 1);
  }
}
data {
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
  real skew;
}
model {
  nu ~ exponential(nu_rate);
  lambda ~ normal(0, 10);
  sigma_proc ~ cauchy(0, 2.5);
  skew ~ cauchy(0, 5);
  for (i in 2:N) {
    y[i] ~ skew_student_t(nu, lambda + b * y[i-1], sigma_proc, skew);
  }
}
'
# stan_gomp_skew2 <- stan_model(model_code = stan_model)
# saveRDS(stan_gomp_skew2, file = "stan-gomp-skew2.rds")

# the base model but with normal process noise:

stan_model <-
'data {
int<lower=0> N; // rows of data
vector[N] y; // vector to hold observations
real b_lower;
real b_upper;
}
parameters {
real lambda;
real<lower=b_lower, upper=b_upper> b;
real<lower=0> sigma_proc;
}
model {
lambda ~ normal(0, 10);
sigma_proc ~ cauchy(0, 2.5);
for (i in 2:N) {
y[i] ~ normal(lambda + b * y[i-1], sigma_proc);
}
}
'
stan_gomp_normal <- stan_model(model_code = stan_model)
saveRDS(stan_gomp_normal, file = "stan-gomp-normal.rds")

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

# Logistic-ricker model, alternative parameterization
stan_model <-
'data {
  int<lower=0> N; // rows of data
  vector[N] y; // log(numbers) at time t
  real<lower=0> K_upper;
  real<lower=0> r_upper;
  real<lower=0> nu_rate; // rate parameter for nu exponential prior
}
parameters {
  real<lower=0, upper=K_upper> K;
  real<lower=0, upper=r_upper> r;
  real<lower=0> sigma_proc;
  real<lower=2> nu;
}
model {
  sigma_proc ~ cauchy(0, 2.5);
  nu ~ exponential(nu_rate);
  for (i in 2:N) {
    y[i] ~ student_t(nu, y[i-1] + r * (1 - exp(y[i-1]) / K), sigma_proc);
  }
}
'
stan_logistic2 <- stan_model(model_code = stan_model)
saveRDS(stan_logistic2, file = "stan-logistic2.rds")

## test
# x <- vector(length = 40)
# x[1] <- 2
# r <- 0.6
# K <- 20
# for (i in 2:40) {
#   x[i] <- x[i-1] + r * (1 - exp(x[i-1]) / K) + rnorm(1, 0, 0.2)
# }
# model <- readRDS("stan-logistic2.rds")
# sm <- sampling(model, data = list(N = length(x), y = x, K_upper = exp(max(x)) * 2,
#   r_upper = 10, nu_rate = 0.01))
##

# AR1 Gompertz that has been extensively tested
# coded similarly to the Stan 2.4.0 manual
stan_model <-
'data {
  int<lower=0> N; // rows of data
  vector[N] y; // vector to hold observations
  real<lower=0> nu_rate; // rate parameter for nu exponential prior
  real b_lower;
  real b_upper;
  real phi_sd;
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
  phi ~ normal(0, phi_sd);
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
# random walk with drift

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

# and a basic random walk model no drift

stan_model <-
  'data {
  int<lower=0> N;
  vector[N] r_obs;
  real<lower=0> nu_rate;
}
parameters {
  real<lower=0> sigma_proc;
  real<lower=2> nu;
}
model {
  nu ~ exponential(nu_rate);
  sigma_proc ~ cauchy(0, 2.5);
  r_obs ~ student_t(nu, 0, sigma_proc);
}
'

stan_rw <- stan_model(model_code = stan_model)
saveRDS(stan_rw, file = "stan-rw.rds")
