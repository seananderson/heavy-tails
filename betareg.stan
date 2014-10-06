// beta regression with group-level intercepts for
// taxonomic class and order
data {
  int<lower=0> N;
  int<lower=0> n_order;
  int<lower=0> n_class;
  vector[N] x; // predictor
  int<lower=1,upper=n_order> order_id[N];
  int<lower=1,upper=n_class> class_id[N];
  real<lower=0,upper=1> y[N]; // response
}
parameters {
  vector[n_order] a_order;
  vector[n_class] a_class;
  real b;
  real mu_a;
  real<lower=0> sigma_a_order;
  real<lower=0> sigma_a_class;
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
                + b * x[i];
    mu[i] <- inv_logit(Xbeta[i]);
  }
  A <- mu * phi;
  B <- (1.0 - mu) * phi;
}
model {
  mu_a ~ normal(0, 10);
  a_class ~ normal(0, sigma_a_class);
  a_order ~ normal(0, sigma_a_order);
  phi ~ cauchy(0, 2.5);
  b ~ normal(0, 10);
  sigma_a_class ~ cauchy(0, 2.5);
  sigma_a_order ~ cauchy(0, 2.5);
  y ~ beta(A, B);
}
