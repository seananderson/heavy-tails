// Beta regression with group-level intercepts for
// taxonomic class, order, and species
// Parameters and predictors are written out in full for clarity
// *This version only controls for time-series length as a predictor*
data {
  int<lower=0> N; // rows of data
  int<lower=0> n_class; // number of classes
  int<lower=0> n_order; // number of orders
  int<lower=0> n_sp; // number of species
  int<lower=1,upper=n_class> class_id[N];
  int<lower=1,upper=n_order> order_id[N];
  int<lower=1,upper=n_sp> sp_id[N];
  vector[N] x1; // predictor
  real<lower=0,upper=1> y[N]; // response
}
parameters {
  vector[n_class] a_class; // class-level deviates
  vector[n_order] a_order; // order-level deviates
  vector[n_sp] a_sp; // species-level deviates
  real b1; // coefficient
  real mu_a; // global intercept
  real<lower=0> sigma_a_order; // group-level standard deviations
  real<lower=0> sigma_a_class;
  real<lower=0> sigma_a_sp;
  real<lower=0> phi; // dispersion parameter
}
transformed parameters {
  vector[N] Xbeta; // linear predictor
  vector<lower=0, upper=1>[N] mu; // transformed linear predictor
  vector<lower=0>[N] A; // beta dist. parameter
  vector<lower=0>[N] B; // beta dist. parameter
  for (i in 1:N) {
    Xbeta[i] <- mu_a + a_class[class_id[i]]
                + a_order[order_id[i]]
                + a_sp[sp_id[i]]
                + b1 * x1[i];
    mu[i] <- inv_logit(Xbeta[i]);
  }
  A <- mu * phi;
  B <- (1.0 - mu) * phi;
}
model {
  // group-level intercept distributions:
  a_class ~ normal(0, sigma_a_class);
  a_order ~ normal(0, sigma_a_order);
  a_sp ~ normal(0, sigma_a_sp);
  // priors:
  mu_a ~ cauchy(0, 10);
  phi ~ cauchy(0, 10);
  b1 ~ cauchy(0, 2.5);
  sigma_a_class ~ cauchy(0, 2.5);
  sigma_a_order ~ cauchy(0, 2.5);
  sigma_a_sp ~ cauchy(0, 2.5);
  // likelihood:
  y ~ beta(A, B);
}
