# This file looks at quasi-extinction probability with heavy vs. normal tails.
#
# Do downwards heavy tails really matter for extinction risk?
#
# We'll write a little C++ function to simulate from a Gompertz model. We'll do
# this in C++ just to make everything fast:

library("dplyr")
library("Rcpp")
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

# And we'll write a function to simulate a time series, check how many
# generations go below some threshold, and return the data in a form we can use.
#
# For each, we'll be drawing all parameter values from the posterior and
# starting at the state the time series ended. Then we'll continue on for some
# given number of generations.

sim_gomp_extinction <- function(y1, abundance_ts = NULL, lambda = 1, b = 0.5, N = 10L,
  nu = 3, log_skew = 1, sigma_proc = 0.2, bootstrap_residuals = FALSE,
  model = c("heavy-skewd", "normal")) {

  y1 <- log(y1)
  if (bootstrap_residuals) {
    log_abundance_ts <- log(abundance_ts)
    res <- log_abundance_ts[-1] - (lambda + b * log_abundance_ts[-length(abundance_ts)])
    # bootstrap the process residuals:
    proc_error <- sample(res, size = N, replace = TRUE)
  } else {
    if (model[1] == "heavy-skewed") {
      proc_error <- skewt::rskt(N, df = nu, gamma = exp(log_skew)) * sigma_proc
    }
    if (model[1] == "normal") {
      proc_error <- rnorm(N, sd = sigma_proc)
    }
  }

  y <- gompertz(N, lambda, b, proc_error, y1)
  data.frame(log_abundance = y, year = seq_along(y))
}

# OK, let's try applying this to populations:

# some heavy tail IDs to look at:
# ids <- c(6528, 10127, 20579)

source("5-shape-data.R")
# heavy_pops <- dplyr::filter(gomp_hat_skew, log_skew_50 < 0.7, nu_50 < 30)
heavy_pops <- dplyr::filter(gomp_hat_base, p10 >= 0.5)

gpdd <- readRDS("gpdd-clean.rds") # to get the last abundance
skew_samples <- readRDS("skew_samples_heavy40000.rds")
normal_samples <- readRDS("normal_samples_heavy40000.rds")

library("doParallel")
registerDoParallel(cores = parallel::detectCores())

# join in the starting abundances first for speed:
last_abund <- group_by(gpdd, main_id) %>%
  summarise(last_abundance = population_untransformed[n()])
skew_samples <- inner_join(skew_samples, last_abund, by = "main_id")
normal_samples <- inner_join(normal_samples, last_abund, by = "main_id")

ww <- skew_samples %>% group_by(main_id) %>%
  mutate(heavyness = ifelse(median(nu) < 10, "heavy",
    ifelse(median(nu) < 70, "moderate", "normal"))) %>%
  group_by(main_id, heavyness) %>%
  summarise(ci_low_skew = quantile(log_skew, probs = 0.025),
    ci_up_skew = quantile(log_skew, probs = 0.975)) %>%
  mutate(skew_ci_contains_zero =
      ifelse(ci_low_skew < 0 & ci_up_skew > 0, TRUE, FALSE)) %>%
  as.data.frame()
ww <- select(ww, heavyness, skew_ci_contains_zero) %>% table
ww
ww / rowSums(ww) * 100

skew_samples_heavy <- filter(skew_samples, main_id %in% heavy_pops$main_id)
normal_samples_heavy <- filter(normal_samples, main_id %in% heavy_pops$main_id)

set.seed(1)
out_skew <- plyr::ldply(seq_len(5L), function(i) {
  temp <- plyr::ddply(skew_samples_heavy, c("main_id", "sample_ids"), function(x) {
      sim_gomp_extinction(y1 = x$last_abundance, lambda = x$lambda, b = x$b,
        N = 6L, nu = x$nu, sigma_proc = x$sigma_proc, log_skew = x$log_skew,
        model = "heavy-skewed")}, .parallel = FALSE, .progress = "text")
  temp$iteration_id <- i
  temp}, .parallel = FALSE)
saveRDS(out_skew, file = "project-skew-heavy.rds")

out_normal <- plyr::ldply(seq_len(5L), function(i) {
  temp <- plyr::ddply(normal_samples_heavy, c("main_id", "sample_ids"), function(x) {
      sim_gomp_extinction(y1 = x$last_abundance, lambda = x$lambda, b = x$b,
        N = 6L, nu = NULL, sigma_proc = x$sigma_proc, log_skew = NULL,
        model = "normal")}, .parallel = FALSE, .progress = "text")
  temp$iteration_id <- i
  temp}, .parallel = FALSE)
saveRDS(out_normal, file = "project-normal.rds")

out_skew$model <- "heavy-skew"
out_normal$model <- "normal"
out <- bind_rows(out_skew, out_normal)

q <- group_by(out, model, main_id, year) %>%
  mutate(abundance = exp(log_abundance)) %>%
  summarise(
    q_low1 = quantile(abundance, probs = 1/1000),
    q_low2 = quantile(abundance, probs = 1/500),
    q_low3 = quantile(abundance, probs = 1/200),
    q_low4 = quantile(abundance, probs = 1/100),
    q_low5 = quantile(abundance, probs = 1/20),
    med = quantile(abundance, probs = 0.5))

# make fake earlier data to merge in:
before_ts <- filter(gpdd, main_id %in% heavy_pops$main_id) %>%
  select(main_id, series_step, population_untransformed) %>%
  group_by(main_id) %>%
  mutate(year = series_step - max(series_step)) %>%
  select(-series_step) %>%
  rename(med = population_untransformed) %>%
  filter(year < 0) %>%
  mutate(year = year + 1) %>%
  as.data.frame()
before_ts <- bind_rows(
  data.frame(before_ts, model = "heavy-skew", stringsAsFactors = FALSE),
  data.frame(before_ts, model = "normal", stringsAsFactors = FALSE))
before_and_projections <- bind_rows(q, before_ts)

saveRDS(before_and_projections, file = "before_and_projections.rds")

library(ggplot2)
plain_theme <- theme_bw() + theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  strip.background = element_blank())
p <- ggplot(before_and_projections, aes(year, med, colour = model, fill = model)) +
  geom_line() + facet_wrap(~main_id, scales = "free") +
  geom_line(aes(y = q_low3), lty = 2) +
  scale_y_log10() + plain_theme + ylab("Abundance")
ggsave("heavy-skew-projections-ggplot.pdf", width = 16, height = 10)

qq <- before_and_projections %>% group_by(main_id, year) %>%
  summarise(r = q_low3[model == "heavy-skew"]/q_low3[model == "normal"]) %>%
  filter(year == 6)
plot(qq$r, log = "y", ylim = c(0.2, 5))
abline(h = 1)
median(1/qq$r)
plot(density(1/qq$r), xlim = c(-1, 4.5));abline(v = 1)
quantile(1/qq$r, probs = c(0, 0.25, 0.5, 0.75, 1)) %>%
  round(1)
saveRDS(qq, file = "skew-understimates.rds")

