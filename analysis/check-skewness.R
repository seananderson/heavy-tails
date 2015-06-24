# a variety of scrap code from when I was testing the skewed model:

# check if there's bias in the skewness parameter for a standard normal gomertz
# model:

library(rstan)
library(dplyr)
library(ggplot2)
library(Rcpp)
cppFunction("
  NumericVector gompertz(double n, double lambda, double b,
    NumericVector proc_error, double y1) {
    NumericVector y(n);
    y(0) = y1;
    for (int i = 1; i < n; ++i) {
      y(i) = lambda + b * y(i - 1) + proc_error(i - 1);
    }
    return y;
}")
stan_gomp_skew <- readRDS("stan-gomp-skew.rds")

test_skew_gompertz <- function(skew_param, N = 20) {
  set.seed(1234)
  out <- plyr::llply(seq_len(N), function(i) {
    message(i)
    xx <- gompertz(50, lambda = 1.1, b = 0.5,
      proc_error = skewt::rskt(50, df = 1e6, gamma = skew_param), y1 = 8)
    d <- list(N = length(xx), y = xx, nu_rate = 0.01,
      b_lower = -1, b_upper = 2)
    sm_skew <- sampling(stan_gomp_skew, data = d,
      pars = c("nu", "lambda", "b", "sigma_proc", "log_skew"),
      iter = 1500, chains = 4, warmup = 750)
    sm_skew})
  out2 <- plyr::ldply(out, function(x) {
     e <- extract(x)$log_skew %>%
       quantile(probs = c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)) %>%
         t %>% as.data.frame
     e})
  out2$iter <- 1:nrow(out2)
  out2
}

# no skew:
noskew <- test_skew_gompertz(skew_param = 1)
downskew <- test_skew_gompertz(skew_param = 0.6)
upskew <- test_skew_gompertz(skew_param = 1.6)

noskew$gamma_param <- 1
downskew$gamma_param <- 0.6
upskew$gamma_param <- 1.6

skew_sim <- bind_rows(noskew, downskew) %>% bind_rows(upskew) %>% as.data.frame
skew_sim$skew_param_label <- paste0("skew = ", skew_sim$gamma_param)

plain_theme <- theme_bw() + theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  strip.background = element_blank())

ticks <- c(0.2, 0.5, 1, 2, 5)

p <- ggplot(skew_sim, aes(iter, `50%`, group = gamma_param)) +
  geom_segment(
    aes(x = 1, xend = max(skew_sim$iter),
        y = log(gamma_param), yend = log(gamma_param)),
      col = "grey60", lwd = 0.3, lty = 2) +
  geom_segment(aes(yend = `97.5%`, y = `2.5%`, xend = iter), lwd = 0.3)  +
  geom_segment(aes(yend = `75%`, y = `25%`, xend = iter), lwd = 0.8) +
  geom_point() + theme_bw() +
  ylab("Skewness parameter") + xlab("Test iteration") +
  facet_wrap(~skew_param_label) + plain_theme +
  coord_cartesian(ylim = c(-2.0, 2.0)) +
  scale_y_continuous(breaks = log(ticks), labels = ticks)
  #print(p)
ggsave("skew-gompertz-simulation.pdf", width = 6, height = 2.6)

###############################


## # for testing skew-t:
## stan_model <-
## '
## functions {
##   real skew_student_t_log(real y, real nu, real mu, real sigma, real skew) {
##   real lp;
##   if (skew <= 0)
##     reject("Skew has to be positive. Found skew=", skew);
##   if (sigma <= 0)
##     reject("Scale has to be positive.  Found sigma=", sigma);
##   lp <- log(skew) - log1p(square(skew));
##   if (y < mu)
##     return lp + student_t_log(y * skew, nu, mu * skew, sigma);
##   else
##     return lp + student_t_log(y / skew, nu, mu / skew, sigma);
##   }
## }
## data {
##   int<lower=0> N; // rows of data
##   vector[N] y; // vector to hold observations
##   real<lower=0> nu_rate; // rate parameter for nu exponential prior
## }
## parameters {
##   real<lower=0> sigma_proc;
##   real<lower=2> nu;
##   real log_skew;
##   real mu;
## }
## model {
##   nu ~ exponential(nu_rate);
##   sigma_proc ~ cauchy(0, 2.5);
##   mu ~ cauchy(0, 2.5);
##   log_skew ~ cauchy(0, 2.5);
##   for (i in 2:N) {
##     y[i] ~ skew_student_t(nu, mu, sigma_proc, exp(log_skew));
##   }
## }
## '
## stan_skew <- stan_model(model_code = stan_model)
## saveRDS(stan_skew, file = "stan-skew.rds")
##
##
## gpdd <- readRDS("gpdd-clean.rds")
## d <- subset(gpdd, main_id == "10127")
## d <- subset(gpdd, main_id == "6528")
## d <- subset(gpdd, main_id == "20579")
##
## d <- subset(gpdd, main_id == "10040")
##
## library(ggplot2)
## ggplot(d, aes(sample_year, population_untransformed)) + geom_line()
## d <- list(N = nrow(d), y = log(d$population_untransformed), nu_rate = 0.01,
##   b_lower = -1, b_upper = 2)
##
## #x <- skewt::rskt(300, df = 1e6, gamma = 1.0)
## library(dplyr)
## set.seed(1234)
## out <- plyr::llply(seq_len(20), function(i) {
##   x <- rnorm(40)
##   plot(x)
##   d <- list(N = length(x), y = x, nu_rate = 0.01)
##   s <- sampling(stan_skew, data = d,
##     pars = c("nu", "mu", "sigma_proc", "log_skew"),
##     iter = 600, chains = 4, warmup = 300)
##   s
##   }, .parallel = FALSE)
##
## out2 <- plyr::ldply(out, function(x) {
##    e <- extract(x)$log_skew %>%
##      quantile(probs = c(0.025, 0.25, 0.5, 0.75, 0.975)) %>%
##        t %>% as.data.frame
##    e})
## names(out2)<-c("ll", "l", "m", "u", "uu")
## out2$iter <- 1:nrow(out2)
##
## #ggplot(out2, aes(iter, log10(m))) + geom_pointrange(aes(ymax = log10(u), ymin = log10(l)), lwd = 0.8) + coord_flip() + geom_pointrange(aes(ymax = log10(uu), ymin = log10(ll)))  + geom_hline(yintecept = 0)
## ggplot(out2, aes(iter, (m))) + geom_pointrange(aes(ymax = (u), ymin = (l)), lwd = 0.8) + coord_flip() + geom_pointrange(aes(ymax = (uu), ymin = (ll)))  + geom_hline(yintecept = 0)
## x <- seq(-1, 1, length.out = 2000);plot(x, dcauchy(x, scale = 2.5), type = "l", ylim = c(0, 0.13))
##
## # look at the skew-t distribution:
## library(skewt)
## params <- expand.grid(x = seq(-10, 10, length.out = 100), df = c(2, 5, 100), gamma = round(exp(c(-0.6, 0.3, 0, 0.3, 0.6)), 1))
## o <- plyr::mdply(params, dskt)
## ggplot(o, aes(x, V1)) + geom_line() + facet_grid(df~gamma) + ylab("Probability density") +
##   xlab("x")
##
## set.seed(3)
## xx <- gompertz(100, 0.8, 0.9, proc_error = rnorm(100, 0, sd = 0.6), y1 = 8)
## d <- list(N = length(xx), y = xx, nu_rate = 0.01,
##   b_lower = -1, b_upper = 2)
## plot(xx)
##
## sm <- sampling(stan_gomp, data = d,
##   pars = c("nu", "lambda", "b", "sigma_proc"),
##   iter = 1000, chains = 4, warmup = 500)
##
## sm_skew <- sampling(stan_gomp_skew, data = d,
##   pars = c("nu", "lambda", "b", "sigma_proc", "log_skew"),
##   iter = 1000, chains = 3, warmup = 500)
##
## #sm_skew2 <- sampling(stan_gomp_skew2, data = d,
##   #pars = c("nu", "lambda", "skew", "b", "sigma_proc"),
##   #iter = 800, chains = 4, warmup = 400)
##
## source("5-shape-data.R")
##
## heavy <- filter(gomp_hat_base, nu_50 < 10) %>%
##   select(main_id, common_name, p10)
##
## check_t <- function(id) {
##   d <- subset(gpdd, main_id == id)
##   d <- list(N = nrow(d), y = log(d$population_untransformed), nu_rate = 0.01,
##     b_lower = -1, b_upper = 2)
##   sm_skew <- sampling(stan_gomp_skew, data = d,
##     pars = c("nu", "lambda", "skew", "b", "sigma_proc"),
##     iter = 3000, chains = 4, warmup = 1500)
##   e <- extract(sm_skew)
##   data.frame(
##     main_id = id,
##     b = e$b,
##     lambda = e$lambda,
##     sigma_proc = e$sigma_proc,
##     skew = e$skew,
##   nu_skew = e$nu,
##   max_rhat = max(summary(sm_skew)$summary[, "Rhat"]),
##   min_neff = min(summary(sm_skew)$summary[, "n_eff"]))
## }
##
## library(doParallel)
## registerDoParallel(cores = 2)
## out <- plyr::ldply(unique(heavy$main_id), function(i) check_t(i), .parallel = TRUE)
##
## skews <- group_by(out, main_id) %>% summarise(median_skew = median(skew)) %>%
##   arrange(median_skew) %>% mutate(order_ = 1:length(median_skew))
## out2 <- inner_join(out, skews) %>% mutate(main_id = reorder(main_id, order_))
##
## p <- ggplot(out2, aes(log10(nu_skew))) + geom_histogram() +
##   facet_wrap(~main_id)  + geom_vline(xintercept = log10(10))
## ggsave("skew-t-nu.pdf", width = 13, height = 9)
## p <- ggplot(out2, aes(skew, group = main_id)) + geom_histogram(binwidth = 0.05) +
##   facet_wrap(~main_id) + geom_vline(xintercept = 1, lty = 2) +
##   scale_x_log10()
## ggsave("skew-t-skew.pdf", width = 13, height = 9)
##
## group_by(out2, main_id) %>% summarise(p10 = sum(nu_skew <= 10)/length(nu_skew), median_skew = median(skew))  %>% as.data.frame() %>% arrange(p10)
##
##
##
##
