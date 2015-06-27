# This file looks at quasi-extinction probability with heavy vs. normal tails.
#
# **I'm currently adapting this script to bring in the parameter values from the
# model fitting.**
#
# Do downwards heavy tails really matter for extinction risk?
#
# Two ways to go about this: repeatedly bootstrap from the observed process
# deviations or draw the deviations from a t distribution with given sigma and
# nu. The problem with the later is that the deviations won't be appropriately
# downwardly skewed. But the problem with the former is that I'm not sure this
# gets what we want out of the exercise. The assumed normal distributions might
# not in fact be normal. The model might just have fit poorly. Therefore, I'm
# going to start with the first option.
#
# We'll write a little C++ function to simulate from a Gompertz model. We'll do
# this in C++ just to make everything fast:

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

sim_gomp_extinction <- function(abundance_ts, lambda = 1, b = 0.5, N = 10L,
  nu = 3, log_skew = 1, sigma_proc = 0.2, bootstrap_residuals = FALSE,
  model = c("heavy-skewd", "normal")) {

  log_abundance_ts <- log(abundance_ts)
  y1 <- log_abundance_ts[length(abundance_ts)]

  if (bootstrap_residuals) {
    res <- log_abundance_ts[-1] - (lambda + b * log_abundance_ts[-length(abundance_ts)])
    # bootstrap the process residuals:
    proc_error <- sample(res, size = N, replace = TRUE)
  } else {
    if (model[1] == "heavy-skewed") {
      proc_error <- skewt::rskt(N, df = nu, gamma = exp(log_skew)) * sigma_proc
    } else {
      proc_error <- rnorm(N, sd = sigma_proc)
    }
  }

  y <- gompertz(N, lambda, b, proc_error, y1)
  data.frame(log_abundance = y, year = seq_along(y))
}

# OK, let's try applying this to populations:

# some heavy tail IDs to look at:
# ids <- c(6528, 10127, 20579)

gpdd <- readRDS("gpdd-clean.rds") # to get the last abundance
skew_samples <- readRDS("skew_samples.rds")
normal_samples <- readRDS("normal_samples.rds")

out_skew <- plyr::ddply(skew_samples, "main_id", function(x) {
  id_i <- x$main_id[1L]
  message(paste("main_id:", id_i))
  this_dat <- dplyr::filter(gpdd, main_id == id_i)
  out_inner <- plyr::ddply(x, c("main_id", "sample_ids"), function(xx) {
    sim_gomp_extinction(abundance_ts = this_dat$population_untransformed,
      lambda = xx$lambda[1L], b = xx$b[1L], N = 10L, nu = xx$nu[1L],
      sigma_proc = xx$sigma_proc[1L], log_skew = xx$log_skew[1L],
      model = "heavy")
  })
  out_inner
})

out_normal <- plyr::ddply(normal_samples, "main_id", function(x) {
  id_i <- x$main_id[1L]
  message(paste("main_id:", id_i))
  this_dat <- dplyr::filter(gpdd, main_id == id_i)
  out_inner <- plyr::ddply(x, c("main_id", "sample_ids"), function(xx) {
    sim_gomp_extinction(abundance_ts = this_dat$population_untransformed,
      lambda = xx$lambda[1L], b = xx$b[1L], N = 10L, nu = NULL,
      sigma_proc = xx$sigma_proc[1L], log_skew = NULL,
      model = "normal")
  })
  out_inner
})

# out_df2 <- mutate(out_df2, projection = ifelse(generation > 0, TRUE, FALSE))
# out_df2 <- filter(out_df2, exp(dat) < 10e5)
library("ggplot2")

cvar <- function(x) mean(x[x<=quantile(x, probs = 0.01)])
# group_by(out_df2, main_id) %>% summarise(cvar_ = cvar(exp(dat)))
#
# out_df3 <- out_df2 %>% group_by(main_id, generation, projection) %>%
#   summarise(
#     min_dat = min(dat, na.rm = TRUE),
#     q99 = quantile(dat, probs = 0.01, na.rm = TRUE)
#     )
#
# p <- ggplot(out_df3,
#   aes(generation, exp(min_dat), colour = projection)) +
#   geom_line() +
#   geom_line(aes(y = exp(q99)), lty = 2) +
#   geom_vline(xintercept = 0, lty = 2) +
#   geom_hline(yintercept = 10, lty = 2) +
#   xlim(-30, 10) +
#   theme_bw() + ylab("Abundance") + xlab("Year") +
#   scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
#     labels = scales::trans_format("log10", scales::math_format(10^.x)),
#     limits = c(1, 1e4)) +
#   annotation_logticks(sides = "l") +
#   guides(colour = guide_legend(override.aes= list(alpha = 1))) +
#   facet_wrap(~main_id)
# p
