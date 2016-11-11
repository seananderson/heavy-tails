// Beta regression with group-level intercepts for
// taxonomic class, order, and species
// Uncertainty specified around the predictors x1, x2, x3
data {
  int<lower=0> N; // rows of data
  int<lower=0> n_class; // number of classes
  int<lower=0> n_order; // number of orders
  int<lower=0> n_sp; // number of species
  int<lower=1,upper=n_class> class_id[N];
  int<lower=1,upper=n_order> order_id[N];
  int<lower=1,upper=n_sp> sp_id[N];
  vector[N] x1; // predictor
  vector[N] x2; // predictor
  vector[N] x3; // predictor
  vector[N] x4; // predictor
  vector[N] x5; // predictor
  real<lower=0,upper=1> y[N]; // response
  vector<lower=0>[N] x1_sigma; // sd around predictor
  vector<lower=0>[N] x2_sigma; // sd around predictor
  vector<lower=0>[N] x3_sigma; // sd around predictor
}
parameters {
  vector[n_class] a_class_z; // class-level deviates
  vector[n_order] a_order_z; // order-level deviates
  vector[n_sp] a_sp_z; // species-level deviates
  real b1; // coefficients
  real b2;
  real b3;
  real b4;
  real b5;
  real mu_a; // global intercept
  real<lower=0> sigma_a_order; // group-level standard deviations
  real<lower=0> sigma_a_class;
  real<lower=0> sigma_a_sp;
  real<lower=0> phi; // dispersion parameter
  vector[N] x1_true; // unobserved true predictor values
  vector[N] x2_true;
  vector[N] x3_true;
}
transformed parameters {
  vector[n_class] a_class; // class-level deviates
  vector[n_order] a_order; // order-level deviates
  vector[n_sp] a_sp; // species-level deviates
  vector[N] Xbeta; // linear predictor
  vector<lower=0, upper=1>[N] mu; // transformed linear predictor
  vector<lower=0>[N] A; // beta dist. parameter
  vector<lower=0>[N] B; // beta dist. parameter

  a_class = sigma_a_class * a_class_z;
  a_order = sigma_a_order * a_order_z;
  a_sp = sigma_a_sp * a_sp_z;

  for (i in 1:N) {
    Xbeta[i] = mu_a + a_class[class_id[i]]
                + a_order[order_id[i]]
                + a_sp[sp_id[i]]
                + b1 * x1_true[i]
                + b2 * x2_true[i]
                + b3 * x3_true[i]
                + b4 * x4[i]
                + b5 * x5[i];
    mu[i] = inv_logit(Xbeta[i]);
  }
  A = mu * phi;
  B = (1.0 - mu) * phi;
}
model {
  // group-level intercept distributions:
  a_class_z ~ normal(0, 1);
  a_order_z ~ normal(0, 1);
  a_sp_z ~ normal(0, 1);
  // measurement error:
  x1 ~ normal(x1_true, x1_sigma);
  x2 ~ normal(x2_true, x2_sigma);
  x3 ~ normal(x3_true, x3_sigma);
  // priors:
  mu_a ~ student_t(3, 0, 10);
  phi ~ student_t(3, 0, 25);
  b1 ~ student_t(3, 0, 2.5);
  b2 ~ student_t(3, 0, 2.5);
  b3 ~ student_t(3, 0, 2.5);
  b4 ~ student_t(3, 0, 2.5);
  b5 ~ student_t(3, 0, 2.5);
  sigma_a_class ~ student_t(3, 0, 2.5);
  sigma_a_order ~ student_t(3, 0, 2.5);
  sigma_a_sp ~ student_t(3, 0, 2.5);
  // likelihood:
  y ~ beta(A, B);
}
