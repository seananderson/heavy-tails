// beta regression with group-level intercepts for
// taxonomic class and order
data {
  int<lower=0> N;
  int<lower=0> n_order;
  int<lower=0> n_class;
  vector[N] x1; // predictor
  vector[N] x2; // predictor
  vector[N] x3; // predictor
  vector[N] x4; // predictor
  vector[N] x5; // predictor
  int<lower=1,upper=n_order> order_id[N];
  int<lower=1,upper=n_class> class_id[N];
  real<lower=0,upper=1> y[N]; // response
  vector<lower=0>[N] x1_sigma;
  vector<lower=0>[N] x2_sigma;
  vector<lower=0>[N] x3_sigma;
}
parameters {
  vector[n_order] a_order;
  vector[n_class] a_class;
  real b1;
  real b2;
  real b3;
  real b4;
  real b5;
  real mu_a;
  real<lower=0> sigma_a_order;
  real<lower=0> sigma_a_class;
  real<lower=0> phi; // dispersion parameter
  vector[N] x1_true;
  vector[N] x2_true;
  vector[N] x3_true;
}
transformed parameters {
  vector[N] Xbeta; // linear predictor
  vector<lower=0, upper=1>[N] mu; // transformed linear predictor
  vector<lower=0>[N] A; // beta dist. parameter
  vector<lower=0>[N] B; // beta dist. parameter
  for (i in 1:N) {
    Xbeta[i] <- mu_a + a_class[class_id[i]]
                + a_order[order_id[i]]
                + b1 * x1_true[i]
                + b2 * x2_true[i]
                + b3 * x3_true[i]
                + b4 * x4[i]
                + b5 * x5[i];
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
  x1 ~ normal(x1_true, x1_sigma);
  x2 ~ normal(x2_true, x2_sigma);
  x3 ~ normal(x3_true, x3_sigma);
  b1 ~ normal(0, 10);
  b2 ~ normal(0, 10);
  b3 ~ normal(0, 10);
  b4 ~ normal(0, 10);
  b5 ~ normal(0, 10);
  sigma_a_class ~ cauchy(0, 2.5);
  sigma_a_order ~ cauchy(0, 2.5);
  y ~ beta(A, B);
}
