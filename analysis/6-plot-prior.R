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
  col = "black", log = "", ylim = NULL, ylab = "", label = "", ...) {
  if(is.null(ylim)) {
    ylim <- c(0, max(y)*1.04)
  }
  if(!add) {
    plot(x, y, type = "l", lty = lty, xlab = xlab, ylab = ylab,
      ylim = ylim, yaxs = "i", xaxs = "i", xlim = xlim, col = col,
      log = log, ...)
    mtext(label, side = 3, line = -1.2, font = 2, adj = 0.05)
  } else {
    lines(x, y, lty = lty, col = col)
  }
}

pdf("priors-gomp-base.pdf", width = 6, height = 5.2)
par(mfrow = c(2, 2), mar = c(3,3,0,0), oma = c(.5, 3, .5, .5),
  tck = -0.02, mgp = c(1.5, 0.4, 0), col.axis = "grey25", col = "grey25", las = 1)
x <- seq(-9, 9, length.out = 200)
plot_prior(x, dnorm(x, 0, 10), expression(lambda), label = "(a)")
x <- seq(0, 6, length.out = 200)
plot_prior(x, dcauchy(x, 0, 2.5), expression(sigma[proc]), ylim = c(0, 0.16), label = "(b)")
x <- seq(2, 500, length.out = 2000)
plot_prior(x, dtrunc(x, "exp", a = 2, b = Inf, rate = 0.02), expression(nu),
  col = "grey70", log = "", ylim = NULL, lty = "93", label = "(c)")
plot_prior(x, dtrunc(x, "exp", a = 2, b = Inf, rate = 0.01),
  expression(nu), add = TRUE, lty = 1, col = "black")
plot_prior(x, dtrunc(x, "exp", a = 2, b = Inf, rate = 0.005),
  expression(nu), add = TRUE, lty = "93", col = "grey20")

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
  text(0.5, 0.4/100, "Weaker\nprior", pos = 4, cex = 0.9)

},
  x = c(200, 500),
  y = c(0.01, 0.02))

x <- seq(-1.1, 1.1, length.out = 200)
plot_prior(x, dtrunc(x, "norm", a = -1, b = 1, mean = 0, sd = 0.5),
  expression(phi), xlim = c(-1.1, 1.1), label = "(d)")
mtext("Probability density", side = 2, outer = TRUE, line = 1, las = 0, cex = 0.9)
dev.off()
