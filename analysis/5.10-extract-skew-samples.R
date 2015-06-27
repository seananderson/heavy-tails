sample_posterior <- function(path = ".", N = 1000L) {
  f <- list.files(path, pattern = "*.rds")
  samples <- plyr::ldply(f[1:3], function(x) {

    message(x)
    e <- rstan::extract(readRDS(paste0(path, x)))
    n_samples <- length(e$b) # check length of an arbitrary coefficient
    id <- as.numeric(gsub("sm-([0-9]+)\\.rds", "\\1", x))
    coefs <- names(e)
    s_ids <- sample(seq_len(n_samples), size = N, replace = FALSE)

    out <- data.frame(
      main_id    = rep(id, N),
      sample_ids = s_ids,
      sigma_proc = e$sigma_proc[s_ids],
      b          = e$b[s_ids],
      lambda     = e$lambda[s_ids],
      n_samples  = rep(n_samples, N))

    if ("nu" %in% coefs) {
      data.frame(out,
        log_skew = e$log_skew[s_ids],
        nu       = e$nu[s_ids])
    } else {
      out
    }
  })

  samples
}

set.seed(123)
skew_samples   <- sample_posterior("~/scratch/heavy/gomp-base-skew/")
normal_samples <- sample_posterior("~/scratch/heavy/gomp-base-normal/")
saveRDS(skew_samples,   file = "skew_samples.rds")
saveRDS(normal_samples, file = "normal_samples.rds")
