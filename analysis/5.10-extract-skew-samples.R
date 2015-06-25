library(dplyr)
f <- list.files("~/scratch/heavy/gomp-base-skew/", pattern = "*.rds")
skew_samples <- plyr::ldply(f, function(x) {
  N <- 1000L
  message(x)
  e <- readRDS(paste0("~/scratch/heavy/gomp-base-skew/", x)) %>% rstan::extract()
  n_samples_total <- length(e$log_skew)
  if (n_samples_total >= N) {
    s <- base::sample(e$log_skew, size = N, replace = FALSE)
  } else {
    s <- base::sample(e$log_skew, size = N, replace = TRUE)
  }
  id <- gsub("sm-([0-9]+)\\.rds", "\\1", x) %>% as.numeric
  data.frame(main_id = rep(id, length(s)), log_skew = s,
    n_samples_total = rep(n_samples_total, length(s)))
})
saveRDS(skew_samples, file = "skew_samples.rds")
