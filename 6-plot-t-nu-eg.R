# Check how student-t changes with nu
# to justify tailing off prior around 100 or 200
pdf("t-nu.pdf", width = 5, height = 4.5)
cols <- c(rev(RColorBrewer::brewer.pal(5, "YlOrRd"))[-5], "grey30")
lwd <- 2.0
par(mfrow = c(2, 1))
par(mar = c(0, 4, 0, 1), oma = c(4, 0, 1, 0), las = 1)
par(xpd = FALSE)
par(cex = 0.8, tck = -0.025, mgp = c(2, 0.6, 0))
x <- seq(-9.5, 9.5, length.out = 500)
plot(x, dt(x, 2), type = "l", xlab = "Standard deviation",
  ylab = "", yaxs = "i", lty = 1, , lwd = lwd,
  xlim = c(-8, 8), ylim = c(0, 0.42), col = cols[1], xaxt = "n", yaxt = "n")
axis(2, at = c(0, 0.2, 0.4))
lines(x, dt(x, 3), lty = 1, col = cols[2], lwd = lwd)
lines(x, dt(x, 10), lty = 1, col = cols[3], lwd = lwd)
lines(x, dt(x, 50), lty = 1, col = cols[4], lwd = lwd)
lines(x, dt(x, 1e9), lty = 1, col = cols[5], lwd = lwd)
legend("topright", legend = c(expression(nu==2), expression(nu==3),
  expression(nu==10), expression(nu==50), expression(nu==infinity~(Normal))),
  col = cols, bty = "n", lty = c(1, 1, 1, 1, 1),
  lwd = c(lwd, lwd, lwd, lwd, lwd))
# zoomed panel:
plot(x, dt(x, 2), type = "l", xlab = "Standard deviation",
  ylab = "", yaxs = "i", lty = 1, , lwd = lwd,
  xlim = c(-8, 8), ylim = c(0, 0.025), col = cols[1], yaxt = "n")
axis(2, at = c(0, 0.02))
lines(x, dt(x, 3), lty = 1, col = cols[2], lwd = lwd)
lines(x, dt(x, 10), lty = 1, col = cols[3], lwd = lwd)
lines(x, dt(x, 50), lty = 1, col = cols[4], lwd = lwd)
lines(x, dt(x, 1e6), lty = 1, col = cols[5], lwd = lwd)
par(xpd = NA)
mtext("Probability density", side = 2, outer = TRUE, line = -1.2, las = 0)
mtext("Standard deviations", side = 1, outer = FALSE, line = 2.5, las = 0)
dev.off()
