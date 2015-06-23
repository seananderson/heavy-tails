extract_model <- function(id, get_phi = TRUE,
  root_folder = "/global/scratch/anderson/heavy", sub_folder = "base",
  file_prefix = "sm",
  type = c("gompertz", "logistic", "rate"),
  get_skew = FALSE,
  get_nu = TRUE) {

  library("rstan")

  sm <- readRDS(paste0(root_folder, "/", sub_folder, "/", file_prefix, "-",
    id, ".rds"))

  samp <- extract(sm)
  q <- lapply(samp, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
  q <- do.call("rbind", q)

  get_q <- function(par) {
    x <- q[par, ]
    names(x) <- paste0(par, "_", names(x))
    as.data.frame(t(x))
  }

  if(type[1] == "rate") {
    lambda <- get_q("lambda")
    sigma_proc <- get_q("sigma_proc")
    out <- data.frame(lambda, nu, sigma_proc)
  }
  if(type[1] == "gompertz") {
    lambda <- get_q("lambda")
    b <- get_q("b")
    sigma_proc <- get_q("sigma_proc")
    out <- data.frame(lambda, nu, b, sigma_proc)
  }
  if(type[1] == "logistic") {
    r <- get_q("r")
    K <- get_q("K")
    sigma_proc <- get_q("sigma_proc")
    out <- data.frame(r, nu, K, sigma_proc)
  }

  sm_summ <- summary(sm)$summary
  max_rhat <- max(sm_summ[, "Rhat"])
  min_neff <- min(sm_summ[, "n_eff"])

  if(get_phi) {
    phi <- get_q("phi")
  }
  if(get_skew) {
    skew <- get_q("log_skew")
  }
  if(get_nu) {
    nu <- get_q("nu")
    nu_samples <- samp$nu
    p10 <- length(nu_samples[nu_samples < 10]) / length(nu_samples)
    p20 <- length(nu_samples[nu_samples < 20]) / length(nu_samples)
    nu_rhat <- sm_summ["nu", "Rhat"]
    nu_neff <- sm_summ["nu", "n_eff"]
  } else {
    p10 <- NA
    p20 <- NA
    nu_rhat <- NA
    nu_neff <- NA
  }

  iter <- sm@stan_args[[1]]$iter

  if(get_phi) {
    out <- data.frame(out, phi)
  }
  if(get_skew) {
    out <- data.frame(out, skew)
  }

  out <- data.frame(out, p10 = p10, p20 = p20, max_rhat = max_rhat,
    min_neff = min_neff, nu_rhat = nu_rhat, nu_neff = nu_neff, iter = iter)

  out$main_id <- id
  names(out) <- sub("\\.", "", names(out))
  out

}
