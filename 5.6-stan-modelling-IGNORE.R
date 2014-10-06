# cut from 5.5-modelling.R

source("5-shape-data.R")

gomp_hat_base$log_Lifesp_scaled <- scale(log(gomp_hat_base$Lifesp))
gomp_hat_base$log_sigma_proc_50_scaled <- scale(log(gomp_hat_base$sigma_proc_50))
gomp_hat_base$b_50_scaled <- scale(gomp_hat_base$b_50)
gomp_hat_base$lambda_50_scaled <- scale(gomp_hat_base$lambda_50)
gomp_hat_base$log_dataset_length_scaled <- scale(log(gomp_hat_base$dataset_length))
gomp_hat_base$log_Len_scaled <- scale(log(gomp_hat_base$Len))
gomp_hat_base$heavy <- ifelse(gomp_hat_base$p10 > 0.5, 1, 0)

temp.dat <- gomp_hat_base[,c("log_Lifesp_scaled", "log_sigma_proc_50_scaled",
  "log_dataset_length_scaled", "b_50_scaled", "lambda_50_scaled",
  "taxonomic_class", "taxonomic_order", "heavy", "p10")]
temp.dat <- na.omit(temp.dat)
temp.dat$taxonomic_class <- as.factor(temp.dat$taxonomic_class)
temp.dat$taxonomic_order <- as.factor(temp.dat$taxonomic_order)

# now try with Stan:

# stan_model <-
# "data {
#   int<lower=1> N; // number of data points
#   //int<lower=1> n_class; // number of classes
#   int<lower=1> n_order; // number of orders
#   vector<lower=0,upper=1>[N] y; // response data
#   vector[N] x; // predictor vector
#   //int<lower=1,upper=n_class> class_id[N];
#   int<lower=1,upper=n_order> order_id[N];
# }
# parameters {
#   //vector[n_class] beta_class;
#   vector[n_order] beta_order;
#   //real mu_beta1;
#   //real<lower=0> tau_class; // sd on class-level intercept
#   real<lower=0> tau_order; // sd on order-level intercept
#   real mu_beta0; // global intercept
#   real beta1; // global slope
# }
# transformed parameters {
#   vector[N] y_hat;
#   for(i in 1:N) {
#     //y_hat[i] <- mu_beta0 + beta_order[order_id[i]] + x[i] * beta1;
#     y_hat[i] <- mu_beta0 + x[i] * beta1;
#   }
# }
# model {
#   mu_beta0 ~ normal(0, 10);
#   //tau_class ~ cauchy(0, 5);
#   tau_order ~ cauchy(0, 5);
# //  beta_class ~ normal(0, 10);
#   beta_order ~ normal(0, 10);
#   beta1 ~ normal(0, 10);
#   y ~ binomial_logit(y_hat);
# }
# "

stan_model <-
"data {
  int<lower=0> N;
  int<lower=0> n_order;
  vector[N] x;
  int<lower=1,upper=n_order> order[N];
  int<lower=0,upper=1> y[N];
}
parameters {
  vector[n_order] a_order;
  real b;
  real mu_a;
  real sigma_a;
}
transformed parameters {
  vector[N] y_hat;
  for (i in 1:N)
    y_hat[i] <- mu_a + a[order[i]] + b * x[i];
}
model {
  mu_a ~ normal(0, 10);
  a ~ normal(0, sigma_a);
  b ~ normal(0, 10);
  sigma_a ~ cauchy(0, 2.5);
  y ~ bernoulli_logit(y_hat);
}"

library(rstan)
if(!file.exists("stan_logistic_reg.rds")) {
  stan_logistic <- stan_model(model_code = stan_model)
  saveRDS(stan_logistic, file = "stan_logistic_reg.rds")
} else {
  stan_logistic <- readRDS("stan_logistic_reg.rds")
}

#m.admb.1 <- glmmadmb(heavy ~ log_Lifesp_scaled + log_sigma_proc_50_scaled +
  #log_dataset_length_scaled + b_50_scaled + lambda_50_scaled + (1 |
    #taxonomic_class / taxonomic_order), family = "binomial", data = temp.dat)

#d <- temp.dat[sample(1:nrow(temp.dat), 100), ]
d <- temp.dat
d$order_id <- as.numeric(d$taxonomic_order)
d$class_id <- as.numeric(d$taxonomic_class)

library(lme4)
m.lmer.1 <- glmer(heavy ~ log_dataset_length_scaled + (1 | taxonomic_order), family = "binomial", data = d)

m.stan.1 <- sampling(stan_logistic,
  data = list(
    N = nrow(d),
    n_order = max(d$order_id),
    x = as.numeric(d$log_dataset_length_scaled),
    order = d$order_id,
    y = d$heavy),
  pars = c("b", "mu_a", "sigma_a"),
  iter = 2000, chains = 3)
sink("m.stan.1.txt")
print(m.stan.1)
sink()
saveRDS(m.stan.1, file = "m.stan.1.rds")

# now with class-level predictor:

stan_model <-
"data {
  int<lower=0> N;
  int<lower=0> n_order;
  int<lower=0> n_class;
  vector[N] x;
  int<lower=1,upper=n_order> order_id[N];
  int<lower=1,upper=n_class> class_id[N];
  int<lower=0,upper=1> y[N];
  int<lower=1,upper=N> id[N];
}
parameters {
  vector[n_order] a_order;
  vector[n_class] a_class;
  vector[N] a_id;
  real b;
  real mu_a;
  real<lower=0> sigma_a_order;
  real<lower=0> sigma_a_class;
  real<lower=0> sigma_a_id;
}
transformed parameters {
  vector[N] y_hat;
  for (i in 1:N)
    y_hat[i] <- mu_a + a_class[class_id[i]]
                + a_order[order_id[i]]
                + a_id[id[i]] + b * x[i];
}
model {
  mu_a ~ normal(0, 10);
  a_class ~ normal(0, sigma_a_class);
  a_order ~ normal(0, sigma_a_order);
  a_id ~ normal(0, sigma_a_id);
  b ~ normal(0, 10);
  sigma_a_class ~ cauchy(0, 2.5);
  sigma_a_order ~ cauchy(0, 2.5);
  sigma_a_id ~ cauchy(0, 2.5);
  y ~ bernoulli_logit(y_hat);
}"

library(rstan)
if(!file.exists("stan_logistic_reg.rds")) {
  stan_logistic2 <- stan_model(model_code = stan_model)
  saveRDS(stan_logistic2, file = "stan_logistic_reg2.rds")
} else {
  stan_logistic2 <- readRDS("stan_logistic_reg2.rds")
}

m.lmer.2 <- glmer(heavy ~ log_dataset_length_scaled + (1 | taxonomic_class / taxonomic_order), family = "binomial", data = d)

m.stan.3 <- sampling(stan_logistic2,
  data = list(
    N = nrow(d),
    n_order = max(d$order_id),
    n_class = max(d$class_id),
    id = 1:nrow(d),
    x = as.numeric(d$log_dataset_length_scaled),
    order_id = d$order_id,
    class_id = d$class_id,
    y = d$heavy),
  pars = c("b", "mu_a", "sigma_a_class", "sigma_a_order", "sigma_a_id"),
  iter = 1000, chains = 3)
sink("m.stan.1.txt")
print(m.stan.1)
sink()
saveRDS(m.stan.1, file = "m.stan.1.rds")


stan_beta <- stan_model("betareg.stan")
saveRDS(stan_beta, file = "stan-beta.rds")

m.stan.beta <- sampling(stan_beta,
  data = list(
    N = nrow(d),
    n_order = max(d$order_id),
    n_class = max(d$class_id),
    x = as.numeric(d$log_dataset_length_scaled),
    order_id = d$order_id,
    class_id = d$class_id,
    y = d$p10),
  pars = c("b", "mu_a", "sigma_a_class", "sigma_a_order", "phi"),
  iter = 1000, chains = 3)
