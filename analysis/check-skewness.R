# a variety of scrap code from when I was testing the skewed model:

# for testing skew-t:
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
}
parameters {
  real<lower=0> sigma_proc;
  real<lower=2> nu;
  real log_skew;
  real mu;
}
model {
  nu ~ exponential(nu_rate);
  sigma_proc ~ cauchy(0, 2.5);
  mu ~ cauchy(0, 2.5);
  log_skew ~ cauchy(0, 2.5);
  for (i in 2:N) {
    y[i] ~ skew_student_t(nu, mu, sigma_proc, exp(log_skew));
  }
}
'
stan_skew <- stan_model(model_code = stan_model)
saveRDS(stan_skew, file = "stan-skew.rds")


gpdd <- readRDS("gpdd-clean.rds")
d <- subset(gpdd, main_id == "10127")
d <- subset(gpdd, main_id == "6528")
d <- subset(gpdd, main_id == "20579")

d <- subset(gpdd, main_id == "10040")

library(ggplot2)
ggplot(d, aes(sample_year, population_untransformed)) + geom_line()
d <- list(N = nrow(d), y = log(d$population_untransformed), nu_rate = 0.01,
  b_lower = -1, b_upper = 2)

#x <- skewt::rskt(300, df = 1e6, gamma = 1.0)
library(dplyr)
set.seed(1234)
out <- plyr::llply(seq_len(20), function(i) {
  x <- rnorm(40)
  plot(x)
  d <- list(N = length(x), y = x, nu_rate = 0.01)
  s <- sampling(stan_skew, data = d,
    pars = c("nu", "mu", "sigma_proc", "log_skew"),
    iter = 600, chains = 4, warmup = 300)
  s
  }, .parallel = FALSE)

out2 <- plyr::ldply(out, function(x) {
   e <- extract(x)$log_skew %>%
     quantile(probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) %>%
       t %>% as.data.frame
   e})
names(out2)<-c("ll", "l", "m", "u", "uu")
out2$iter <- 1:nrow(out2)

#ggplot(out2, aes(iter, log10(m))) + geom_pointrange(aes(ymax = log10(u), ymin = log10(l)), lwd = 0.8) + coord_flip() + geom_pointrange(aes(ymax = log10(uu), ymin = log10(ll)))  + geom_hline(yintecept = 0)
ggplot(out2, aes(iter, (m))) + geom_pointrange(aes(ymax = (u), ymin = (l)), lwd = 0.8) + coord_flip() + geom_pointrange(aes(ymax = (uu), ymin = (ll)))  + geom_hline(yintecept = 0)
x <- seq(-1, 1, length.out = 2000);plot(x, dcauchy(x, scale = 2.5), type = "l", ylim = c(0, 0.13))

# look at the skew-t distribution:
library(skewt)
params <- expand.grid(x = seq(-10, 10, length.out = 100), df = c(2, 5, 100), gamma = round(exp(c(-0.6, 0.3, 0, 0.3, 0.6)), 1))
o <- plyr::mdply(params, dskt)
ggplot(o, aes(x, V1)) + geom_line() + facet_grid(df~gamma) + ylab("Probability density") +
  xlab("x")

set.seed(3)
xx <- gompertz(100, 0.8, 0.9, proc_error = rnorm(100, 0, sd = 0.6), y1 = 8)
d <- list(N = length(xx), y = xx, nu_rate = 0.01,
  b_lower = -1, b_upper = 2)
plot(xx)

sm <- sampling(stan_gomp, data = d,
  pars = c("nu", "lambda", "b", "sigma_proc"),
  iter = 1000, chains = 4, warmup = 500)

sm_skew <- sampling(stan_gomp_skew, data = d,
  pars = c("nu", "lambda", "b", "sigma_proc", "log_skew"),
  iter = 1000, chains = 3, warmup = 500)

#sm_skew2 <- sampling(stan_gomp_skew2, data = d,
  #pars = c("nu", "lambda", "skew", "b", "sigma_proc"),
  #iter = 800, chains = 4, warmup = 400)

source("5-shape-data.R")

heavy <- filter(gomp_hat_base, nu_50 < 10) %>%
  select(main_id, common_name, p10)

check_t <- function(id) {
  d <- subset(gpdd, main_id == id)
  d <- list(N = nrow(d), y = log(d$population_untransformed), nu_rate = 0.01,
    b_lower = -1, b_upper = 2)
  sm_skew <- sampling(stan_gomp_skew, data = d,
    pars = c("nu", "lambda", "skew", "b", "sigma_proc"),
    iter = 3000, chains = 4, warmup = 1500)
  e <- extract(sm_skew)
  data.frame(
    main_id = id,
    b = e$b,
    lambda = e$lambda,
    sigma_proc = e$sigma_proc,
    skew = e$skew,
  nu_skew = e$nu,
  max_rhat = max(summary(sm_skew)$summary[, "Rhat"]),
  min_neff = min(summary(sm_skew)$summary[, "n_eff"]))
}

library(doParallel)
registerDoParallel(cores = 2)
out <- plyr::ldply(unique(heavy$main_id), function(i) check_t(i), .parallel = TRUE)

skews <- group_by(out, main_id) %>% summarise(median_skew = median(skew)) %>%
  arrange(median_skew) %>% mutate(order_ = 1:length(median_skew))
out2 <- inner_join(out, skews) %>% mutate(main_id = reorder(main_id, order_))

p <- ggplot(out2, aes(log10(nu_skew))) + geom_histogram() +
  facet_wrap(~main_id)  + geom_vline(xintercept = log10(10))
ggsave("skew-t-nu.pdf", width = 13, height = 9)
p <- ggplot(out2, aes(skew, group = main_id)) + geom_histogram(binwidth = 0.05) +
  facet_wrap(~main_id) + geom_vline(xintercept = 1, lty = 2) +
  scale_x_log10()
ggsave("skew-t-skew.pdf", width = 13, height = 9)

group_by(out2, main_id) %>% summarise(p10 = sum(nu_skew <= 10)/length(nu_skew), median_skew = median(skew))  %>% as.data.frame() %>% arrange(p10)




