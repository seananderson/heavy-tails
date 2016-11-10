# get the means and sds of predictors:
# for modelling later
gpdd <- readRDS("gpdd-clean.rds")
out <- plyr::ldply(unique(gpdd$main_id), function(x) {
  folder <- paste0("gomp-base/sm-", x, ".rds")

  library(rstan)
  sm <- readRDS(folder)

  samp_s <- extract(sm, pars = "sigma_proc")[[1]]
  samp_l <- extract(sm, pars = "lambda")[[1]]
  samp_b <- extract(sm, pars = "b")[[1]]

  mean_log_sigma_proc <- mean(log(samp_s))
  mean_lambda <- mean(samp_l)
  mean_b <- mean(samp_b)

  sd_log_sigma_proc <- sd(log(samp_s))
  sd_lambda <- sd(samp_l)
  sd_b <- sd(samp_b)

  data.frame(mean_log_sigma_proc = mean_log_sigma_proc,
    mean_lambda = mean_lambda,
    mean_b = mean_b,
    sd_log_sigma_proc = sd_log_sigma_proc,
    sd_lambda = sd_lambda,
    sd_b = sd_b, main_id = x)
  }, .progress = "text")
saveRDS(out, file = "gomp-base-mean-sd.rds")
