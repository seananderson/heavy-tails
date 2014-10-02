# This file runs the models with a logistic growth model

source("1.5-compile-fit-function.R")
source("1.6-extract-function.R")

model <- readRDS("stan-logistic.rds")
gpdd <- readRDS("gpdd-clean.rds")
id <- "logistic"

stan_dat <- "list(N = nrow(x) - 1, r_obs = x$r_obs[-1], nu_rate = 0.01, K_upper = 2, r_upper = 10, Nt = x$population_untrans_scaled[-N])"
stan_pars <- c("r", "K", "nu", "sigma_proc")

fit_gpdd_model(gpdd_dat = gpdd, model = model,
  sub_folder = id, stan_dat = stan_dat, pars = stan_pars)

out <- plyr::ldply(unique(gpdd$main_id), extract_model,
  sub_folder = id, get_phi = FALSE, type = "logistic")

saveRDS(out, file = paste0(id, "-hat.rds"))


# stan_logistic <- readRDS("stan-logistic.rds")
# gpdd <- readRDS("gpdd-clean.rds")
# library(rstan)
#
# main_id_vec <- unique(gpdd$main_id)
#
# out_gomp <- plyr::d_ply(subset(gpdd, main_id %in% main_id_vec), "main_id",
#   function(x) {
#
#     max_rhat <- 999
#     min_neff <- 0
#     iterations<- 2000
#     warmup <- 1000
#     file_base <- "/global/scratch/anderson/heavy/logistic-5.0/sm-"
#
#     if(!file.exists(paste0(file_base, unique(x$main_id), ".rds"))) {
#       print("already done")
#       print(unique(x$main_id))
#       # TODO scale the data to max of 1 and calculate the log-diff rate column here
#       x$r_obs <- c(NA, diff(log(x$population_untransformed)))
#       x$population_untrans_scaled <-
#         x$population_untransformed / max(x$population_untransformed)
#       N <- nrow(x) - 1
#
#       while((max_rhat > 1.05 | min_neff < 200) & iterations <= 4000) {
#         sm <- sampling(stan_logistic,
#           data = list(N = N, r_obs = x$r_obs[-1], K_upper = 2, r_upper = 20,
#             nu_rate = 0.01, Nt = x$population_untrans_scaled[-N]),
#           pars = c("r", "K", "nu", "sigma_proc"), iter = iterations,
#           chains = 4, warmup = warmup)
#
#         # check:
#         max_rhat <- max(summary(sm)$summary[, "Rhat"])
#         min_neff <- min(summary(sm)$summary[, "n_eff"])
#
#         warmup <- warmup * 2
#         iterations <- iterations * 2
#       }
#       saveRDS(sm, file = paste0(file_base, unique(x$main_id), ".rds"))
#       sink(paste0(file_base, unique(x$main_id), ".txt"))
#       print(sm)
#       sink()
#     }
#   })
