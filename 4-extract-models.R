# this file extracts posterior and diagnostic data from all the model runs

gpdd <- readRDS("gpdd-clean.rds")

library(rstan)
library(dplyr)

extract_model <- function(id, get_phi = FALSE,
  folder = "/global/scratch/anderson/heavy/", subfolder = "base",
  type = "gompertz") {

  sm <- readRDS(paste0(folder, subfolder, "/sm-", id, ".rds"))
  samp <- extract(sm)
  q <- lapply(samp, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
  q <- do.call("rbind", q)

  get_q <- function(par) {
    x <- q[par, ]
    names(x) <- paste0(par, "_", names(x))
    as.data.frame(t(x))
  }

  if(type[1] == "gompertz") {
    lambda <- get_q("lambda")
    nu <- get_q("nu")
    b <- get_q("b")
    sigma_proc <- get_q("sigma_proc")
    sigma_proc <- get_q("sigma_proc")
    out <- data.frame(lambda, nu, b, sigma_proc)
  }
  if(type[1] == "logistic") {
    r <- get_q("r")
    nu <- get_q("nu")
    K <- get_q("K")
    sigma_proc <- get_q("sigma_proc")
    out <- data.frame(r, nu, K, sigma_proc)
  }

  if(get_phi) {
    phi <- get_q("phi")
  }

  nu_samples <- samp$nu
  p10 <- length(nu_samples[nu_samples < 10]) / length(nu_samples)
  p20 <- length(nu_samples[nu_samples < 20]) / length(nu_samples)

  sm_summ <- summary(sm)$summary
  max_rhat <- max(sm_summ[, "Rhat"])
  min_neff <- min(sm_summ[, "n_eff"])
  nu_rhat <- sm_summ["nu", "Rhat"]
  nu_neff <- sm_summ["nu", "n_eff"]

  iter <- sm@stan_args[[1]]$iter


  if(get_phi) {
    out <- data.frame(out, phi)
  }

  out <- data.frame(out, p10 = p10, p20 = p20, max_rhat = max_rhat,
    min_neff = min_neff, nu_rhat = nu_rhat, nu_neff = nu_neff, iter = iter)

  out$main_id <- id
  names(out) <- sub("\\.", "", names(out))
  out

}

ids <- unique(gpdd$main_id)
#mammals_ids <- gpdd %>% filter(taxonomic_class == "Mammalia") %>%
  #group_by(main_id) %>% summarise(id = main_id[1]) %>% select(main_id)

#gomp_hat_base <- plyr::ldply(ids, extract_model, subfolder = "base-5.0",
  #get_phi = FALSE)
gomp_hat_ar1 <- plyr::ldply(ids, extract_model, subfolder = "ar1-5.0",
  get_phi = TRUE)
gomp_hat_logistic <- plyr::ldply(ids, extract_model, type = "logistic",
  subfolder = "logistic-5.0", get_phi = FALSE)
#gomp_hat_ar1_obs <- plyr::ldply(mammals_ids$main_id, extract_model,
#  subfolder = "ar1-obs0.2", get_phi = TRUE)
#saveRDS(gomp_hat_base, file = "gomp_hat_base.rds")
saveRDS(gomp_hat_ar1, file = "gomp_hat_ar1-5.0.rds")
saveRDS(gomp_hat_logistic, file = "gomp_hat_logistic-5.0.rds")
#saveRDS(gomp_hat_ar1_obs, file = "gomp_hat_ar1_obs0.2.rds")
