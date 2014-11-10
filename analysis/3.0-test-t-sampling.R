# this file tests the combination of N samples and degree of heavy tailedness
# in terms of the frequency that nu is estimated correctly
# this is part 1 of 2 simulation testing checks

library("rstan")
stan_t <- readRDS("stan-t.rds")

# take a basic t- dist
# and try and re-capture it with our prior at different nu values
# use a ton of data (2000 points, 500 points, 100 points, 50 points)
# show how probability of noticing goes down
set.seed(123)

N_vec <- c(1600, 800, 400, 200, 100, 50, 25)
sigma_proc <- 1
nu_vec <- c(1, 3, 5, 10, 1e6)
reps <- 1:20
show_plot <- FALSE

# test:
#N_vec <- 50
#nu_vec <- 100
#reps <- 1:2

out <- plyr::ldply(reps, function(k) {
  out_N <- plyr::ldply(seq_along(nu_vec), function(i) {

    y <- metRology::rt.scaled(max(N_vec), df = nu_vec[i], mean = 0, sd = sigma_proc)

    if(show_plot) {
      plot(y, ylim = c(-6, 6), xaxs = "i", log = "x")
      abline(v = N_vec, col = "red", lwd = 2)
      rect(1, -100, max(N_vec) + 10, qnorm(0.001), col = "#00000030", border = NA)
      rect(1, 100, max(N_vec) + 10, qnorm(0.999), col = "#00000030", border = NA)
      #rect(1, -100, max(N_vec) + 10, qnorm(0.01), col = "#00000030", border = NA)
      #rect(1, 100, max(N_vec) + 10, qnorm(0.99), col = "#00000030", border = NA)
    }

    out_nu <-  plyr::ldply(seq_along(N_vec), function(j) {
      sm <- sampling(stan_t,
        data = list(N = N_vec[j], y = y[1:N_vec[j]], nu_rate = 0.05),
        pars = c("nu", "sigma_proc"), iter = 2000,
        chains = 4, warmup = 1000)
      max_rhat <- max(summary(sm)$summary[, "Rhat"])
      min_neff <- min(summary(sm)$summary[, "n_eff"])
      e <- extract(sm, pars = "nu")[[1]]
      e_sig <- extract(sm, pars = "sigma_proc")[[1]]
      data.frame(med_nu = median(e), l_nu = quantile(e, probs = 0.1),
        u_nu = quantile(e, probs = 0.9),
        l_sigma = quantile(e_sig, probs = 0.1),
        u_sigma = quantile(e_sig, probs = 0.9),
        med_sigma = median(e_sig),
        max_rhat = max_rhat, min_neff = min_neff,
        nu_true = nu_vec[i], N = N_vec[j])
})
    out_nu
})
  out_N$iter <- k
  out_N
})

saveRDS(out, file = "sample-t-sim-check.rds")
