# Try estimating nu with a weakly informative prior in Stan
# See how the ability to detect heavy tails deteriorates with sample size,
# observation error and the degree of heavy-tailedness.

# truncated distribution function from:
# http://www.jstatsoft.org/v16/c02/paper
dtrunc <- function(x, spec, a = -Inf, b = Inf, ...) {
  tt <- rep(0, length(x))
  g <- get(paste("d", spec, sep = ""), mode = "function")
  G <- get(paste("p", spec, sep = ""), mode = "function")
  tt[x>=a & x<=b] <- g(x[x>=a&x<=b], ...)/(G(b, ...) - G(a, ...))
  tt
}

# Check priors:
plot_prior <- function(x, y, xlab, xlim = NULL, add = FALSE, lty = 1,
  col = "black", log = "", ylim = NULL, ylab = "", ...) {
  if(is.null(ylim)) {
    ylim <- c(0, max(y)*1.04)
  }
  if(!add) {
    plot(x, y, type = "l", lty = lty, xlab = xlab, ylab = ylab,
      ylim = ylim, yaxs = "i", xaxs = "i", xlim = xlim, col = col,
      log = log, ...)
  } else {
    lines(x, y, lty = lty, col = col)
  }
}

 pdf("priors-gomp-base.pdf", width = 6, height = 6)
# par(mfrow = c(2, 2), cex = 0.8, las = 1, mar = c(4.5, 3.5, 0, 0),
#   oma = c(1, 1, 1, 1))
par(mfrow = c(2, 2), mar = c(3,3,0,0), oma = c(.5, 3, .5, .5),
  tck = -0.02, mgp = c(1.5, 0.4, 0), col.axis = "grey25", col = "grey25", las = 1)
x <- seq(-9, 9, length.out = 200)
plot_prior(x, dnorm(x, 0, 10), expression(lambda))
x <- seq(0, 6, length.out = 200)
plot_prior(x, dcauchy(x, 0, 2.5), expression(sigma[proc]))
x <- seq(2, 500, length.out = 2000)
plot_prior(x, dtrunc(x, "exp", a = 2, b = Inf, rate = 0.02), expression(nu),
  col = "grey70", log = "", ylim = NULL, lty = "93")
plot_prior(x, dtrunc(x, "exp", a = 2, b = Inf, rate = 0.01),
  expression(nu), add = TRUE, lty = 1, col = "black")
plot_prior(x, dtrunc(x, "exp", a = 2, b = Inf, rate = 0.005),
  expression(nu), add = TRUE, lty = "93", col = "grey20")

# text(1.7, 1.2/100, "Base prior", pos = 4, cex = 0.9)
# text(1, 1.2/100, "Stonger prior", cex = 0.9, pos = 4)
# text(1, 0.2/100, "Weaker\nprior", pos = 4, cex = 0.9)

TeachingDemos::subplot({
  x <- seq(2, 500, length.out = 2000)
  plot_prior(log(x), dtrunc(x, "exp", a = 2, b = Inf, rate = 0.02), "",
    col = "grey70", log = "", ylim = NULL, lty = "93", yaxt = "n", xaxt = "n",
    ylab = "")
  plot_prior(log(x), dtrunc(x, "exp", a = 2, b = Inf, rate = 0.01),
    expression(nu), add = TRUE, lty = 1, col = "black")
  plot_prior(log(x), dtrunc(x, "exp", a = 2, b = Inf, rate = 0.005),
    expression(nu), add = TRUE, lty = "93", col = "grey20")
   axis(2, at = c(0, 0.01, 0.02))
   axis(1, labels = c(2, 5, 20, 200), at = log(c(2, 5, 20, 200)))
  text(0.5, 1.2/100, "Base prior", pos = 4, cex = 0.9)
  text(0.5, 1.8/100, "Stonger prior", cex = 0.9, pos = 4)
  text(0.5, 0.2/100, "Weaker\nprior", pos = 4, cex = 0.9)

},
  x = c(200, 500),
  y = c(0.01, 0.02))

x <- seq(-1.1, 1.1, length.out = 200)
plot_prior(x, dtrunc(x, "norm", a = -1, b = 1, mean = 0, sd = 1),
  expression(phi), xlim = c(-1.1, 1.1))
mtext("Probability density", side = 2, outer = TRUE, line = 1, las = 0, cex = 0.9)
dev.off()

# How much probability mass is below 10 and above 2 in the prior?
# x <- rexp(2e6, rate = 0.01)
#hist(x)
# round(length(x[x<10 & x >2])/length(x), 2)
# round(length(x[x>10])/length(x), 2)
# round(length(x[x<50 & x >2])/length(x), 2)
# round(length(x[x>50])/length(x), 2)

# # Compile Stan model:
# mstan2.obs.nu <-
# "data {
#   int<lower=0> N; // rows of data
#   vector[N] y; // vector to hold observations
#   real<lower=0> sigma_obs; // observation error standard deviation
#   real nu_rate; // rate parameter for nu exponential prior
# }
# parameters {
#   vector[N] t; // state vector of `true` values
#   real a; // intercept (mu)
#   real<lower=2> nu; // student-t degrees of freedom
#   real<lower=0> sigma_proc; // process noise standard deviation
# }
# model {
#   // priors:
#   a ~ normal(0, 5);
#   sigma_proc ~ cauchy(0, 5);
#   nu ~ exponential(nu_rate);
#   t ~ student_t(nu, a, sigma_proc);
#   y ~ normal(t, sigma_obs);
# }
# // generated quantities {
# //   vector[N] log_lik;
# //   real dev; //deviance
# //   dev <- 0.0;
# //   for (i in 1:N) {
# //     log_lik[i] <-
# //       student_t_log(t[i], nu, a, sigma_proc) +
# //       normal_log(o[i], t[i], sigma_obs);
# //     dev <- dev + (-2) * log_lik[i];
# //   }
# //}
# "
# write(mstan2.obs.nu, file = "mstan2.obs.nu.stan")

## #library(rstan)
## #mstan2.obs.nu.c <- stan_model("mstan2.obs.nu.stan")
##
## Nx <- 50
## sigma_proc <- 0.2
## iter <- 1000
## chains <- 4
## warmup <- 500
## nu <- 3
## sigma_obs_vec <- c(0.01, 0.1, 0.2)
## dat <- list()
## dat_obs <- list()
## sm <- list()
## sigma_obs <- 0.2
## for(i in seq_along(sigma_obs_vec)) {
##   set.seed(123)
##   dat[[i]] <- metRology::rt.scaled(Nx, df = nu, mean = 0, sd = sigma_proc)
##   dat_obs[[i]] <- rnorm(Nx, mean = dat[[i]], sd = sigma_obs_vec[i])
##   sm[[i]] <- sampling(mstan2.obs.nu.c,
##     data = list(N = Nx, y = dat_obs[[i]], nu_rate = 0.005,
##       sigma_obs = sigma_obs),
##     pars = c("a", "sigma_proc", "nu"), iter = iter,
##     chains = chains, warmup = warmup)
## }
##
## # NOTE if nu is detectable as very small than it is informative and overpowers the prior; if the data are not informative about nu (because nu is not detectably small) then the prior will show through; in the second case, nu could really be anything over 10 or 20 and without a very large quantity of data, the data will never be informative about that.
##
## pdf("nu-estimation-example.pdf", width = 6.5, height = 5.0)
## par(mfrow = c(3, 2), cex = 0.8)
## par(xpd = FALSE, mar = c(0, 3, 0, 0), oma = c(4, 0, 1.0, 1.0),
##   mgp = c(2, 0.5, 0), tck = -0.05)
## for(i in 1:3) {
##   plot(dat[[i]], col = "grey40", type = "l", xaxt = "n", yaxt = "n", ylab = "")
##   if(i == 3) axis(1)
##   if(i == 2) mtext("Abundance", side = 2, line = 1.0)
##   points(dat_obs[[i]], pch = 21, col = "grey20", bg = "grey60", cex = 0.9)
##   nu.samples <- extract(sm[[i]], pars = "nu")[[1]]
##   if(i == 3) mtext("Time", side = 1, line = 2.3)
##   h <- hist(nu.samples, breaks = seq(2, max(nu.samples)+20,6),
##     xlim = c(2, 150), xlab = "nu hat", xaxt = "n", main = "", col = "white",
##     xaxs = "i", yaxs = "i", plot = FALSE, warn.unused = FALSE)
##   xlim <- c(2, 150)
##   ylim  <- c(0, 1500)
##   plot(1, 1, type = "n", xlim = xlim, yaxs = "i", xaxs = "i",
##     xlab = quote(widehat(nu)), xaxt = "n", ylim = ylim, yaxt = "n", ylab = "")
##   if(i == 3) axis(1, at = c(2, seq(50, 200, 50)))
##   cols <- rev(RColorBrewer::brewer.pal(9, "YlOrRd"))
##   locs <- c(2, seq(10, 100, 3))
##   for(j in 1:9) {
##     rect(locs[j], 0, locs[j+1], max(ylim)*1.2, border = NA,
##       col = paste0(cols[j], "80"))
##   }
##   for(j in 1:100) {
##     rect(h$breaks[j], 0, h$breaks[j+1], h$counts[j], border = "black",
##       col = c("white", rep("white", 99))[j], lwd = 1)
##   }
##
##   # Show median and interquartile range?
##   #segments(quantile(nu.samples, probs = 0.25)[[1]],
##   # max(h$counts) * 0.5, quantile(nu.samples, probs = 0.75)[[1]],
##   # max(h$counts) * 0.5)
##   #points(median(nu.samples), max(h$counts) * 0.5)
##   #points(mean(nu.samples), max(h$counts) * 0.5, pch = 4)
##
##   par(new = TRUE)
##   y <- seq(2, max(xlim), length.out = 300)
##   d <- dexp(y, 0.005)
##   plot(y, d, type = "l", col = "grey60", axes = FALSE, xlim = xlim,
##     xlab = "", ylab ="", lwd = 2, ylim = c(0, max(d)*8), xaxs = "i", lty = 2)
##   par(new = FALSE)
##   mtext(paste0(sprintf("%.2f", round(length(nu.samples[nu.samples < 10]) /
##       length(nu.samples), 2)), " probability nu < 10"), line = -2,
##     col = "grey30", cex = 0.9, adj = 0.6)
##   box()
##   if(i == 1)
##     mtext("Prior", side = 1, line = -2.0, adj = 0.9, col = "grey30", cex = 0.9)
##   if(i == 3) mtext(quote(widehat(nu)), side = 1, line = 2.3)
##   if(i == 2) mtext("Frequency", side = 2, line = 1.0)
## }
## dev.off()
##
## # Checks:
## rstan::traceplot(sm[[1]], inc_warmup = F)
## rstan::traceplot(sm[[2]], inc_warmup = F)
## rstan::traceplot(sm[[3]], inc_warmup = F)
## # check Rhat:
## summary(sm[[1]])$summary["lp__", "Rhat"]
## summary(sm[[2]])$summary["lp__", "Rhat"]
## summary(sm[[3]])$summary["lp__", "Rhat"]
## max(summary(sm[[1]])$summary[, "Rhat"])
##
## summary(sm[[1]])$summary["lp__", "n_eff"]
## min(summary(sm[[1]])$summary[, "n_eff"])
## min(summary(sm[[2]])$summary[, "n_eff"])
## min(summary(sm[[3]])$summary[, "n_eff"])
##
## pairs(sm[[1]])
## pairs(sm[[2]])
## pairs(sm[[3]])
