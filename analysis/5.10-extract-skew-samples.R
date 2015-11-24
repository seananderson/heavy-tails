sample_posterior <- function(path = ".", N = 5000L, files = NULL) {
  if (is.null(files)) {
    f <- list.files(path, pattern = "*.rds")
  } else {
    f <- files
  }
  samples <- plyr::ldply(f, function(x) {

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
      out <- data.frame(out, nu = e$nu[s_ids])
    }

    if ("log_skew" %in% coefs) {
      out <- data.frame(out, log_skew = e$log_skew[s_ids])
    }

    out
  })

  samples
}

set.seed(123)
skew_samples   <- sample_posterior("~/scratch/heavy/gomp-base-skew/", N = 5000L)
normal_samples <- sample_posterior("~/scratch/heavy/gomp-base-normal/", N = 5000L)
saveRDS(skew_samples,   file = "skew_samples.rds")
saveRDS(normal_samples, file = "normal_samples.rds")

source("5-shape-data.R")
h1 <- dplyr::filter(gomp_hat_skew, nu_50 < 30)
h2 <- dplyr::filter(gomp_hat_base, p10 >= 0.5)
heavy_ids <- union(h1$main_id, h2$main_id)
write.csv(heavy_ids, file = "heavy_ids.csv", row.names = FALSE)

heavy_ids <- read.csv("heavy_ids.csv")$x
set.seed(123)
skew_samples   <- sample_posterior("~/scratch/heavy/gomp-base-skew/", N = 40000L,
  files = paste0("sm-", heavy_ids, ".rds"))
normal_samples <- sample_posterior("~/scratch/heavy/gomp-base-normal/", N = 40000L,
  files = paste0("sm-", heavy_ids, ".rds"))
saveRDS(skew_samples,   file = "skew_samples_heavy40000.rds")
saveRDS(normal_samples, file = "normal_samples_heavy40000.rds")
