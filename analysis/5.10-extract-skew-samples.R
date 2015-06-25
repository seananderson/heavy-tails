library(dplyr)
f <- list.files("~/scratch/heavy/gomp-base-skew/", pattern = "*.rds")
set.seed(123456789)
skew_samples <- plyr::ldply(f, function(x) {
  N <- 1000L
  message(x)
  e <- readRDS(paste0("~/scratch/heavy/gomp-base-skew/", x)) %>% rstan::extract()

  n_samples_skew <- length(e$log_skew)
  if (n_samples_skew >= N) {
    s_skew <- base::sample(e$log_skew, size = N, replace = FALSE)
  } else {
    s_skew <- base::sample(e$log_skew, size = N, replace = TRUE)
  }

  n_samples_nu <- length(e$nu)
  if (n_samples_nu >= N) {
    s_nu <- base::sample(e$nu, size = N, replace = FALSE)
  } else {
    s_nu <- base::sample(e$nu, size = N, replace = TRUE)
  }

  id <- gsub("sm-([0-9]+)\\.rds", "\\1", x) %>% as.numeric

  data.frame(
    main_id        = rep(id, N),
    log_skew       = s_skew,
    nu             = s_nu,
    n_samples_skew = rep(n_samples_skew, N),
    n_samples_nu   = rep(n_samples_nu, N))
})
saveRDS(skew_samples, file = "skew_samples.rds")
