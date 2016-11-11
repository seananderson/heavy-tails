# a pretty version of nu coeffs and p(nu < 10) for a main paper figure

library("grImport")

classes <- c("Aves","Mammalia", "Insecta", "Osteichthyes")
cols <- c(rev(RColorBrewer::brewer.pal(5, "YlOrRd"))[-c(4, 5)])
pics <- c("silhouettes/aves.eps.xml", "silhouettes/rabbit.eps.xml",
  "silhouettes/fly.eps.xml", "silhouettes/chinook.eps.xml")
pic_width <- c(9, 9, 9, 17)
pic_height <- c(.10, .10, .10, .08)

add_label <- function(xfrac = -0.02, yfrac = -0.04, label = "", pos = 4, ...) {
  u <- par("usr")
  x <- u[1] + xfrac * (u[2] - u[1])
  y <- u[4] - yfrac * (u[4] - u[3])
  text(x, y, label, pos = pos, ...)
}

set.seed(1) # for jittering
pdf("nu-coefs-2.pdf", width = 3.8, height = 9.2)

par(mfrow = c(4, 1), mar = c(1.5,0,0,0), oma = c(1.8, 3.5, 1, .8), tck = -0.02, mgp = c(2, 0.5, 0), col.axis = "grey25", col = "grey25")
par(cex = 0.9)

for(i in 1:length(classes)) {

  xaxt  <- ifelse(i == 4, "s", "n")

  x <- subset(gomp_hat_base, taxonomic_class == classes[i])

  ticks <- c(2, 3, 5, 10, Inf)

  plot(1, 1, type = "n", xlim = c(-0.2, 100.2), ylim = c(0, 0.5), xlab = "", ylab = "", yaxt = "n", xaxs = "i", yaxs = "i", xaxt = xaxt, las = 1)
  #axis(2, at = seq(0, 100, 50), las = 1)
  par(xpd = NA)
  add_label(label = LETTERS[i], font = 2, cex = 1.2)
  add_label(label = classes[i], xfrac = 0.03)
  if(i == 1) {
    mtext("Heavy", side = 2, line = 2.5, adj = 0.90, col = "grey55")
    mtext("tails", side = 2, line = 1.5, adj = 0.85, col = "grey55")
    mtext("Normal", side = 2, line = 2.5, adj = -0.09, col = "grey55")
    mtext("tails", side = 2, line = 1.5, adj = 0, col = "grey55")
  }
  par(xpd = FALSE)
  axis(2, las = 1, at = 1/ticks, labels = ticks)

  cols_df_coefs <- data.frame(
    col = c(rev(RColorBrewer::brewer.pal(6, "YlOrRd"))[-6], "grey50"),
    cuts = c(1.9, 3, 5, 10, 15, 25), stringsAsFactors = FALSE)
  x$coef_col <- as.character(cols_df_coefs$col[findInterval(x$nu_50, cols_df_coefs$cuts)])

  segments(x$sort_id_perc, 1/x$nu_5, x$sort_id_perc, 1/x$nu_95, lwd = 1.1, col = "grey60")
  segments(x$sort_id_perc, 1/x$nu_25, x$sort_id_perc, 1/x$nu_75, lwd = 1.3, col = "grey30")
  points(x$sort_id_perc, 1/x$nu_50, pch = 21, col = "grey10", lwd = 0.5, bg = x$coef_col, cex = 0.85)

  par(xpd = NA)
  p <- readPicture(pics[i])
  if(i != 4) {
    p@paths$path@rgb <- "grey37"
  } else {
    p@paths[[2]]@rgb <- "grey37"  # different vector structure
  }
  grImport::picture(p, 12, 0.37, 12 + pic_width[i], 0.37 + pic_height[i])
  par(xpd = FALSE)

  if(unique(x$taxonomic_class) == "Aves") {
    cutoff <- 4
  } else {
    cutoff <- 0
  }
  x <- x %>% group_by(taxonomic_order) %>%
    mutate(n_pops = length(p10)) %>%
    filter(n_pops >= cutoff)

  p <- x %>% group_by(taxonomic_order) %>% summarise(med_p10 = mean(p10)) %>%
    arrange(med_p10) %>% mutate(p10_order = 1:length(med_p10))
  x <- dplyr::inner_join(x, p, by = "taxonomic_order")

  # in the following, the first 5 are for the heavy points
  # and the last one is for the normal points
  cols_df <- data.frame(
    col = rev((c(rev(RColorBrewer::brewer.pal(6, "YlOrRd"))[-6], "#00000007"))),
    bg = rev(c(rep("grey40", 5), "#33333350")),
    cex = rev(c(rep(0.75, 5), 0.65)),
    cuts = rev(c(0.9, 0.8, 0.7, 0.6, 0.5, 0)), stringsAsFactors = FALSE)

  x$p_col <- as.character(cols_df$col[findInterval(x$p10, cols_df$cuts)])
  x$bg <- as.character(cols_df$bg[findInterval(x$p10, cols_df$cuts)])
  x$cex <- as.numeric(cols_df$cex[findInterval(x$p10, cols_df$cuts)])

  TeachingDemos::subplot({
    plot(x$p10, x$p10_order, type = "n", xlab = "", ylab = "", las = 1,
      bty = "n", yaxt = "n", xaxt = "n", axes = FALSE, xlim = c(0, 1), xaxs = "i")
    axis(2, at = p$p10_order, labels = p$taxonomic_order, las = 1,
      cex = 0.8, cex.axis= 0.8, lwd = 0, col.axis = "grey30", line = -0.36)
    axis(1, at = c(0, 0.5, 1), col = "grey30", lwd = 0.7, tck = -0.04,
      mgp = c(2, 0.25, 0), col.axis = "grey30", cex.axis = 0.8)
    mtext(expression(Pr(nu<10)), side = 1, line = 1.4, cex = 0.8, col = "grey30")
    par(xpd = FALSE)
    abline(h = as.numeric(as.factor(x$taxonomic_order)), col = "grey93", lwd = 0.7, lty = 1)
    par(xpd = NA)
    points(x$p10, jitter(x$p10_order, amount = 0.25), pch = 21, bg = x$p_col,
      col = x$bg, cex = x$cex, lwd = 0.8)
    par(xpd = FALSE)
  },
    x = c(60, 95),
    y = c(0.25, 0.48))
}

par(xpd = NA)
mtext(quote(t~distribution~degrees~of~freedom~(widehat(nu))), outer = TRUE, side = 2, line = 2.2)
mtext("Percentage of populations", outer = TRUE, side = 1, line = 0.4)
par(xpd = FALSE)

dev.off()
