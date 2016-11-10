source("5-shape-data.R")

gomp_hat_base <- arrange(gomp_hat_base, main_id)
gomp_hat_obs_0.2 <- arrange(gomp_hat_obs_0.2, main_id)
gomp_hat_logistic <- arrange(gomp_hat_logistic, main_id)
gomp_hat_rate <- arrange(gomp_hat_rate, main_id)
gomp_hat_rw <- arrange(gomp_hat_rw, main_id)
gomp_hat_ar1 <- arrange(gomp_hat_ar1, main_id)
gomp_hat_gamma <- arrange(gomp_hat_gamma, main_id)

heavy_ids <- filter(gomp_hat_base, p10 >= 0.50)$main_id

check_nc_heavy <- function(dat, heavy_ids) {
  nc <- filter(dat, max_rhat > 1.05 | min_neff < 200)$main_id
  x <- nc[nc %in% heavy_ids]
  if (length(x) > 0) warning(paste(x, collapse = ", "))
}
check_nc_heavy(gomp_hat_obs_0.2, heavy_ids)
check_nc_heavy(gomp_hat_logistic, heavy_ids)
check_nc_heavy(gomp_hat_rate, heavy_ids)
check_nc_heavy(gomp_hat_rw, heavy_ids)
check_nc_heavy(gomp_hat_ar1, heavy_ids)

cols_df <- data.frame(col =
    c(RColorBrewer::brewer.pal(4, "Set3")[c(3, 4, 1, 2)],
      rep("#FFFFFF", 3)),
  taxonomic_class = c("Aves", "Mammalia", "Insecta", "Osteichthyes",
    "Chondrichtyhes", "Crustacea", "Gastropoda"),
  stringsAsFactors = FALSE)

comp_panel <- function(xdat, ydat, xlab, ylab, label = "") {

  not_converged <- which(ydat$max_rhat > 1.05)
  if (length(not_converged) > 0) {
    message(paste(paste(not_converged, collapse = ", "),
      "didn't converge and will not be plotted"))
  }
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
  mtext(label, side = 3, line = -1.2, font = 2, adj = 0.04, cex = 1.2)
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
comp_panel(gomp_hat_base, gomp_hat_logistic, quote(Gompertz~widehat(nu)), quote(Ricker-logistic~widehat(nu)), label = "A")
legend(0.3, 0.63, legend = cols_df$taxonomic_class[1:4], fill = cols_df$col[1:4], bty = "n", ncol = 2)
comp_panel(gomp_hat_base, gomp_hat_ar1, quote(Gompertz~widehat(nu)), quote(Gompertz~AR1~widehat(nu)), label = "B")
comp_panel(gomp_hat_base, gomp_hat_rate, quote(Gompertz~widehat(nu)), quote(Rate~only~widehat(nu)), label = "C")
comp_panel(gomp_hat_base, gomp_hat_obs_0.2, quote(Gompertz~widehat(nu)), quote(Gompertz~obs.~error~widehat(nu)), label = "D")
dev.off()

# check:
# comp_panel(gomp_hat_base, gomp_hat_rw, quote(Gompertz~widehat(nu)), quote(Random~walk~widehat(nu)), label = "Z")

# pdf("gomp-prior-comparison.pdf", width = 7, height = 6.8)
pdf("gomp-prior-comparison.pdf", width = 7.2, height = 4)
par(mfrow = c(1, 2), mar = c(3,3,0,0), oma = c(.5, .5, 3.5, .5),
  tck = -0.02, mgp = c(1.5, 0.4, 0), col.axis = "grey25", col = "grey25", las = 1)
par(cex = 0.9)
comp_panel(gomp_hat_base, gomp_hat_weaker, quote(Gompertz~widehat(nu)~(base~model)), quote(Gompertz~widehat(nu)~(weaker~prior)), label = "A")
comp_panel(gomp_hat_base, gomp_hat_stronger, quote(Gompertz~widehat(nu)~(base~model)), quote(Gompertz~widehat(nu)~(stronger~prior)), label = "B")
# comp_panel(gomp_hat_base, gomp_hat_gamma, quote(Gompertz~widehat(nu)~(base~model)), quote(Gompertz~widehat(nu)~("Gamma"~prior)), label = "C")

# legend(0.3, 0.63, legend = cols_df$taxonomic_class[1:4], fill = cols_df$col[1:4], bty = "n", ncol = 2)
# comp_panel(gomp_hat_base, gomp_hat_ar1, quote(Gompertz~widehat(nu)), quote(Gompertz~AR1~widehat(nu)))
# comp_panel(gomp_hat_base, gomp_hat_rate, quote(Gompertz~widehat(nu)), quote(Rate~only~widehat(nu)))
# comp_panel(gomp_hat_base, gomp_hat_obs_0.2, quote(Gompertz~widehat(nu)), quote(Gompertz~obs.~error~widehat(nu)))
dev.off()

d <- inner_join(select(gomp_hat_base, main_id, p10),
  select(gomp_hat_gamma, main_id, p10) %>% rename(p10_gamma = p10))

ggplot(d, aes(p10, p10_gamma)) + geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  geom_vline(xintercept = 0.5) +
  geom_hline(yintercept = 0.5) +
  geom_smooth(se = F)

filter(d, p10 > 0.5, p10_gamma < 0.5)
filter(d, p10 > 0.5, p10_gamma > 0.5)

filter(d, p10 > 0.5, p10_gamma < 0.25)
filter(d, p10 > 0.5, p10_gamma > 0.25)

x <- seq(2.01, 100, length.out = 200)
plot(x, dgamma(x, shape = 2, scale = 1/0.1), type = "l", col = "black", ylab = "")
par(new = TRUE)
plot(x, dexp(x, 0.01), type = "l")

r <- rgamma(1e7, shape = 2, scale = 1/0.1)
r <- r[r>2]
round(sum(r<10)/length(r), 2) * 100
round(median(r), 2)

r <- rexp(1e7, 0.01)
r <- r[r>2]
round(sum(r<10)/length(r), 3) * 100
round(median(r), 2)

d2 <- inner_join(select(gomp_hat_base, main_id, nu_50),
  select(gomp_hat_gamma, main_id, nu_50) %>% rename(nu_50_gamma = nu_50))

ggplot(d2, aes(nu_50, nu_50_gamma)) + geom_point() +
  geom_abline(intercept = 0, slope = 1) +
  geom_vline(xintercept = 10) +
  geom_hline(yintercept = 10) +
  geom_smooth(se = F)

cor(d$p10, d$p10_gamma, method = "spearman")
