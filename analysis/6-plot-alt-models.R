cols_df <- data.frame(col =
    c(RColorBrewer::brewer.pal(4, "Set3")[c(3, 4, 1, 2)],
      rep("#FFFFFF", 3)),
  taxonomic_class = c("Aves", "Mammalia", "Insecta", "Osteichthyes",
    "Chondrichtyhes", "Crustacea", "Gastropoda"),
  stringsAsFactors = FALSE)

comp_panel <- function(xdat, ydat, xlab, ylab) {

  not_converged <- which(ydat$max_rhat > 1.05)
  if(length(not_converged) > 0) {
    xdat <- xdat[-not_converged, ]
    ydat <- ydat[-not_converged, ]
  }
  # colours:

  xdat <- inner_join(xdat, lookup[,c("main_id", "taxonomic_class")])
  xdat <- inner_join(xdat, cols_df)

  set.seed(1)
  scramble <- sample(1:nrow(xdat))
  xdat <- xdat[scramble, ]
  ydat <- ydat[scramble, ]

  par(xpd = NA)
  plot(1/xdat$nu_50, 1/ydat$nu_50, yaxs = "i", xaxs = "i",
    xlim = c(0, 0.5), ylim = c(0, 0.5), xlab = xlab, ylab = ylab,
    axes = FALSE, pch = 21, bg = xdat$col, fg = "grey40",
    type = "n")
  rect(0.1, 0, 0.5, 0.1, col = "grey90", border = "grey60")
  rect(0, 0.1, 0.1, 0.5, col = "grey90", border = "grey60")
  segments(1/xdat$nu_25, 1/ydat$nu_50, 1/xdat$nu_75, 1/ydat$nu_50,
    col = "#00000035", lwd = 0.6)
  segments(1/xdat$nu_50, 1/ydat$nu_25, 1/xdat$nu_50, 1/ydat$nu_75,
    col = "#00000035", lwd = 0.6)
  points(1/xdat$nu_50, 1/ydat$nu_50, pch = 21, bg = xdat$col,
    fg = "grey40", cex = 1.1)
  box()
  ticks <- c(2, 3, 5, 10, Inf)
  axis(1, at = 1/ticks, labels = ticks)
  axis(2, at = 1/ticks, labels = ticks)
  par(xpd = FALSE)
  abline(a = 0, b = 1)
  par(xpd = NA)
}

pdf("gomp-comparison.pdf", width = 7, height = 6.8)
par(mfrow = c(2, 2), mar = c(3,3,0,0), oma = c(.5, .5, 3.5, .5),
  tck = -0.02, mgp = c(1.5, 0.4, 0), col.axis = "grey25", col = "grey25", las = 1)
par(cex = 0.9)
comp_panel(gomp_hat_base, gomp_hat_logistic, quote(Gompertz~widehat(nu)), quote(Ricker-logistic~widehat(nu)))
legend(0.3, 0.63, legend = cols_df$taxonomic_class[1:4], fill = cols_df$col[1:4], bty = "n", ncol = 2)
comp_panel(gomp_hat_base, gomp_hat_ar1, quote(Gompertz~widehat(nu)), quote(Gompertz~AR1~widehat(nu)))
comp_panel(gomp_hat_base, gomp_hat_rate, quote(Gompertz~widehat(nu)), quote(Rate~only~widehat(nu)))
comp_panel(gomp_hat_base, gomp_hat_obs_0.2, quote(Gompertz~widehat(nu)), quote(Gompertz~obs.~error~widehat(nu)))
dev.off()

# pdf("gomp-prior-comparison.pdf", width = 7, height = 6.8)
pdf("gomp-prior-comparison.pdf", width = 7.2, height = 4)
par(mfrow = c(1, 2), mar = c(3,3,0,0), oma = c(.5, .5, 3.5, .5),
  tck = -0.02, mgp = c(1.5, 0.4, 0), col.axis = "grey25", col = "grey25", las = 1)
par(cex = 0.9)
comp_panel(gomp_hat_base, gomp_hat_weaker, quote(Gompertz~widehat(nu)~(base~model)), quote(Gompertz~widehat(nu)~(weaker~prior)))
comp_panel(gomp_hat_base, gomp_hat_stronger, quote(Gompertz~widehat(nu)~(base~model)), quote(Gompertz~widehat(nu)~(stronger~prior)))

# legend(0.3, 0.63, legend = cols_df$taxonomic_class[1:4], fill = cols_df$col[1:4], bty = "n", ncol = 2)
# comp_panel(gomp_hat_base, gomp_hat_ar1, quote(Gompertz~widehat(nu)), quote(Gompertz~AR1~widehat(nu)))
# comp_panel(gomp_hat_base, gomp_hat_rate, quote(Gompertz~widehat(nu)), quote(Rate~only~widehat(nu)))
# comp_panel(gomp_hat_base, gomp_hat_obs_0.2, quote(Gompertz~widehat(nu)), quote(Gompertz~obs.~error~widehat(nu)))
dev.off()
