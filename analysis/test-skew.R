
check <- function() {
  library(skewt)
  N <- 80
  lambda <- 0.6
  b <- 0.5
  y <- vector(length = N)
  # set.seed(1)
  proc_error <- skewt::rskt(N, df = 100, gamma = 1)
  for(i in 2:N) {
    y[i] <- lambda + b * y[i-1] + proc_error[i-1]
  }
  plot(y)

  sm <- readRDS("stan-gomp-skew.rds")

  library(rstan)
  # options(mc.cores = parallel::detectCores())
  options(mc.cores = 1)
  m <- sampling(sm, data = list(N = N, y = y, nu_rate = 0.01, b_lower = -1, b_upper = 2),
    iter = 400, control = list(adapt_delta = 0.90), chains = 2)

  library(dplyr)
  extract(m)$log_skew %>%
    exp() %>% quantile(probs = c(0.025, 0.5, 0.975)) %>%
    round(2) %>% t %>%
    as.data.frame()
}

d <- plyr::rdply(10, check())
library(ggplot2)
filter(d, .n < 11) %>%
  ggplot(aes(x = .n, y = `50%`, ymin = `2.5%`, ymax = `97.5%`)) +
  geom_pointrange() +
  geom_hline(yintercept = 1)
