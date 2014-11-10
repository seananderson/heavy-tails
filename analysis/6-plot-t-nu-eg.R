# Check how student-t changes with nu
# to justify tailing off prior around 100 or 200
#pdf("t-nu.pdf", width = 5, height = 4.5)

library("rstan")

bs_cols <- c(rev(RColorBrewer::brewer.pal(5, "YlOrRd"))[-5], "grey30")

add_label <- function(xfrac = 0.1, yfrac = 0.1, label = "", pos = 4, ...) {
  u <- par("usr")
  x <- u[1] + xfrac * (u[2] - u[1])
  y <- u[4] - yfrac * (u[4] - u[3])
  text(x, y, label, pos = pos, ...)
}

box.col <- "grey20"
pdf("t-nu-eg2.pdf", width = 7.5, height = 4)
par(col.axis = "grey20", col = "grey20")
l <- rbind(c(1, 1, 1, 3, 3, 6, 6), c(1, 1, 1, 4,4, 7,7), c(2,2,2, 5,5, 8,8))
layout(l)
cols <- c(rev(RColorBrewer::brewer.pal(5, "YlOrRd"))[-5], "grey30")
lwd <- 2.0
par(mar = c(0, 2, 0, 0), oma = c(3.5, 2.5, 1.5, 0.5), las = 1)
par(xpd = FALSE)
par(cex = 0.8, tck = -0.025, mgp = c(2, 0.6, 0))
x <- seq(-9.5, 9.5, length.out = 500)
plot(x, dt(x, 2), type = "l", xlab = "",
  ylab = "", yaxs = "i", lty = 1, , lwd = lwd,
  xlim = c(-8, 8), ylim = c(0, 0.42), col = cols[1], xaxt = "n", yaxt = "n")
text(-7.5, 0.393, quote((a)))
mtext("Student-t distribution", side = 3, line = 0, cex = 0.9)
axis(2, at = c(0, 0.2, 0.4), col = box.col)
lines(x, dt(x, 3), lty = 1, col = cols[2], lwd = lwd)
lines(x, dt(x, 5), lty = 1, col = cols[3], lwd = lwd)
lines(x, dt(x, 20), lty = 1, col = cols[4], lwd = lwd)
lines(x, dt(x, 1e9), lty = 1, col = cols[5], lwd = lwd)
legend("topright", legend = c(expression(nu==2), expression(nu==3),
  expression(nu==5), expression(nu==20), expression(Normal)),
  col = cols, bty = "n", lty = c(1, 1, 1, 1, 1),
  lwd = c(lwd, lwd, lwd, lwd, lwd))
# zoomed panel:
plot(x, dt(x, 2), type = "l", xlab = "Standard deviation",
  ylab = "", yaxs = "i", lty = 1, , lwd = lwd,
  xlim = c(-8, 8), ylim = c(0, 0.025), col = cols[1], yaxt = "n")
text(-7.5, 0.0218, quote((b)))
axis(2, at = c(0, 0.02), col = box.col)
lines(x, dt(x, 3), lty = 1, col = cols[2], lwd = lwd)
lines(x, dt(x, 5), lty = 1, col = cols[3], lwd = lwd)
lines(x, dt(x, 20), lty = 1, col = cols[4], lwd = lwd)
lines(x, dt(x, 1e6), lty = 1, col = cols[5], lwd = lwd)
par(xpd = NA)
mtext("Probability density", side = 2, outer = TRUE, line = 1.2, las = 0,
  cex = 0.9)
mtext("Standard deviations", side = 1, outer = FALSE, line = 2, las = 0,
  cex = 0.9)
par(xpd = FALSE)

load("nu_effective_seeds.rda")

make_ts <- function(seed, nu, ylim = c(-6.5, 9), bs_col = "red") {
  N <- 50
  set.seed(seed)
  proc_error <- metRology::rt.scaled(N, df = nu, mean = 0, sd = 0.65)
  lambda <- 1.5
  y1 <- 3
  b = 0.5
  y <- rep(NA, N)
  y[1] <- y1
  for(i in 2:N) {
    y[i] <- lambda + b * y[i-1] + proc_error[i]
  }
  black_swan_i <- seq_along(proc_error)[pnorm(proc_error) < 0.001]
  col <- rep("grey50", N)
  cex = rep(0.7, N)
  col[black_swan_i] <- bs_col
  cex[black_swan_i] <- 1
  pch = rep(21, N)
  pch[black_swan_i] <- 19
  plot(y, type = "l", ylim = ylim, ylab = "", xlab = "", yaxt = "n",
    xaxt = "n", cex = 0.8, col = "grey20", xaxs = "i")
  points(seq_along(y), y, pch = pch, cex = cex, col = col, bg = "grey70")
  return(list(y = y, nu = nu))
}

seed <- nu_5_seeds_N50$seeds[19]
y1 <- make_ts(seed = seed, nu = 3, bs_col = bs_cols[2])
text(1, 7.1, quote((c)~nu==3), pos = 4)
mtext("Simulated population dynamics", side = 3, line = 0, cex = 0.9)
text(13, -3.9, "Black swans", pos = 4)
y2 <- make_ts(seed = seed, nu = 5, bs_col = bs_cols[3])
mtext("log Abundance", side = 2, line = 0.5, las = 0, cex = 0.9)
#text(24, -2.5, "Gompertz\npopulation\nmodel", pos = 4)
text(1, 7.1, quote((d)~nu==5), pos = 4)
y3 <- make_ts(seed = seed, nu = 1e9, bs_col = "darkgrey")
axis(1, at = c(1, seq(10, 50, 10)), col = box.col)
text(1, 7.1, quote((e)~Normal), pos = 4)
mtext("Time step", side = 1, line = 2, las = 0, cex = 0.9)

if(!file.exists("illustrate-nu-sampling.rda")) {
  stan_gomp <- readRDS("stan-gomp.rds")
  sm1 <- sampling(stan_gomp,
    data = list(N = 50, y = y1$y, nu_rate = 0.01,
      b_lower = -1, b_upper = 2),
    pars = c("lambda", "sigma_proc", "nu", "b"), iter = 2000,
    chains = 4, warmup = 1000)

  sm2 <- sampling(stan_gomp,
    data = list(N = 50, y = y2$y, nu_rate = 0.01,
      b_lower = -1, b_upper = 2),
    pars = c("lambda", "sigma_proc", "nu", "b"), iter = 2000,
    chains = 4, warmup = 1000)

  sm3 <- sampling(stan_gomp,
    data = list(N = 50, y = y3$y, nu_rate = 0.01,
      b_lower = -1, b_upper = 2),
    pars = c("lambda", "sigma_proc", "nu", "b"), iter = 2000,
    chains = 4, warmup = 1000)
  save(sm1, sm2, sm3, file = "illustrate-nu-sampling.rda")
} else {
  load("illustrate-nu-sampling.rda")
}
e1 <- extract(sm1, pars = "nu")[[1]]
e2 <- extract(sm2, pars = "nu")[[1]]
e3 <- extract(sm3, pars = "nu")[[1]]

make_hist2 <- function(x) {
  xlim <- c(2, 88)
  h <- hist(x, breaks = seq(2, max(x)+20,3),
    xlim = xlim, xlab = "", xaxt = "n", main = "", col = "white",
    xaxs = "i", yaxs = "i", plot = FALSE, warn.unused = FALSE)
  ylim  <- c(0,3200)
  plot(1, 1, type = "n", xlim = xlim, yaxs = "i", xaxs = "i",
    xlab = quote(widehat(nu)), xaxt = "n", ylim = ylim, yaxt = "n", ylab = "")
  cols <- rev(RColorBrewer::brewer.pal(9, "YlOrRd"))
  locs <- c(2, seq(10, 100, 1.5))

bs_df <- data.frame(nu = c(2, 8, 20, 30),
   col = c(bs_cols[1], bs_cols[3], bs_cols[4], "#FFFFFF"),
  stringsAsFactors = FALSE)

bs_bg <- data.frame(nu = min(bs_df$nu):(max(bs_df$nu)-1),
  col = "colour", stringsAsFactors = FALSE)

finterp <- function(cols, num) {
  f <- colorRampPalette(cols)
  f(num)
}

bs_bg$col <- unlist(sapply(1:(nrow(bs_df) - 1), function(i)
  finterp(bs_df$col[i:(i+1)], bs_df$nu[i+1] - bs_df$nu[i])))

for(j in 1:(nrow(bs_bg) - 1)) {
  rect(bs_bg$nu[j], 0, bs_bg$nu[j+1], max(ylim)*1.2, border = NA, lwd = 0,
    col = paste0(bs_bg$col[j], "80"))
}

#   cols_df <- data.frame(nu = c(1.9, 2.9, 4.9, 15.9, 25),
#     col = c(bs_cols[-length(bs_cols)], "grey30"),
#     stringsAsFactors = FALSE)
#   h_cols <- cols_df$col[findInterval(h$breaks, cols_df$nu)]

  for(j in seq_along(h$breaks)) {
    rect(h$breaks[j], 0, h$breaks[j+1], h$counts[j], border = "grey10",
#       col = h_cols[j], lwd = 1)
       col = "white", lwd = 1)
  }
  segments(quantile(x, probs = 0.25), ylim[2] / 2,
    quantile(x, probs = 0.75), ylim[2] / 2, col = "grey10")
  points(median(x), ylim[2] / 2, pch = 21, col = "grey10", bg = "grey70",
    lwd = 0.7)
  z <- seq(xlim[1], xlim[2], length.out = 400)
  par(new = TRUE)
  plot(z, dexp(z, 0.01), type = "l", ylim = c(0, 0.025), xlab = "", ylab = "",
    axes = FALSE, xlim = xlim, xaxs = "i", lty = 2)
  par(new = FALSE)
}
make_hist2(e1)
text(2.5, 0.022, quote((f)~Heavy~tails), pos = 4)
text(30, 0.009, "Prior", pos = 4)
text(8, 0.0015, "Posterior samples", pos = 4)
mtext("Model fits", side = 3, line = 0, cex = 0.9)
make_hist2(e2)
text(2.5, 0.022, quote((g)~Slightly~heavy~tails), pos = 4)
mtext("Samples or probability", side = 2, line = 0.5, las = 0, cex = 0.9)
make_hist2(e3)
axis(1, at = c(2, seq(20, 80, 20)), col = box.col)
mtext(quote(widehat(nu)), side = 1, line = 2, las = 0, cex = 0.9)
text(37, 0.015, "Median + IQR", pos = 4)
text(2.5, 0.022, quote((h)~Normal~tails), pos = 4)
dev.off()
