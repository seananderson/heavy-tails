# a pretty version of the correlates of nu plot for the main paper
# need to run 5.8... first

# randomly sort the rows:
set.seed(1)
gomp_hat_base_corr <- gomp_hat_base[sample(1:nrow(gomp_hat_base)), ]

gomp_hat_base_corr$col <- NULL
cols_df <-
  data.frame(col = c(RColorBrewer::brewer.pal(4, "Set3")[c(3, 4, 1, 2)],
    rep("#FFFFFF", 3)),
  taxonomic_class = c("Aves", "Mammalia", "Insecta", "Osteichthyes",
    "Chondrichtyhes", "Crustacea", "Gastropoda"),
  stringsAsFactors = FALSE)
gomp_hat_base_corr <- plyr::join(gomp_hat_base_corr, cols_df)

add_label <- function(xfrac = 0.03, yfrac = 0.04, label = "", pos = 4, ...) {
  u <- par("usr")
  x <- u[1] + xfrac * (u[2] - u[1])
  y <- u[4] - yfrac * (u[4] - u[3])
  text(x, y, label, pos = pos, ...)
}

make_panel <- function(x, xl = NULL, xu = NULL, log = "", yaxis = FALSE,
  label = "",
  xaxis_ticks = NULL) {

  i <<- i + 1
  y <- gomp_hat_base_corr$p10
  class_col <- gomp_hat_base_corr$col

  y <- y[!is.na(x)]
  class_col <- class_col[!is.na(x)]
  x <- x[!is.na(x)]

  plot(1, 1, xlim = range(x), ylim = c(0, 1.02), type = "n", xlab = "",
    ylab = "", log = log, yaxs = "i", axes = FALSE)

  if(!is.null(xl) & !is.null(xu)) {
    segments(xl, y, xu, y, col = paste0(class_col, "50"), lwd = 0.9)
  }

  points(x, y, pch = 21, col = "#66666690", lwd = 0.8,
    bg = paste0(class_col, ""), cex = 1)

  if(i == 4) {
    with(pred_dataset_length,
      lines(dataset_length, X50., lwd = 2.8, col = "#00000072"))
    with(pred_dataset_length,
      polygon(c(dataset_length, rev(dataset_length)),
        c(X5., rev(X95.)), border = FALSE, col = "#00000025"))
  }
  if(i == 1) {
    with(pred_sigma,
      lines(sigma_proc_mean, X50., lwd = 2.8, col = "#00000072"))
    with(pred_sigma,
      polygon(c(sigma_proc_mean, rev(sigma_proc_mean)),
        c(X5., rev(X95.)), border = FALSE, col = "#00000025"))
  }

  if(is.null(xaxis_ticks))
    axis(1)
  else
    axis(1, at = xaxis_ticks)
  if(yaxis) {
    axis(2, las = 1, at = seq(0, 1, 0.5))
  }
  box()
  par(xpd = NA)
  lab <- letters[i]
# legend("topleft", legend = substitute(paste("(", lab, ") ", label)),
  legend("topleft", legend = substitute(paste("(", lab, ")")),
    bty = "n", inset = c(-0.1, -0.03), cex = 0.9)
  mtext(substitute(label), side = 1, line = 1.5)
  par(xpd = FALSE)
}

pdf("correlates-p10.pdf", width = 7, height = 4.6)
par(mfrow = c(2, 3), mar = c(3.2,0,0,0), oma = c(0, 4.2,1.1, 5.4),
  tck = -0.02, mgp = c(2, 0.3, 0), col.axis = "grey25", col = "grey25")
par(cex = 0.9)

i <<- 0

with(gomp_hat_base_corr,
  make_panel(x = sigma_proc_50, xl = sigma_proc_25, xu = sigma_proc_75,
    log = "x", yaxis = TRUE, label = paste("Gompertz ", sigma),
    xaxis_ticks = c(0.05, 0.2, 1, 5)))

mtext("Heavy tails", side = 2, line = 1.8, adj = 1.15, col = "grey55",
  cex = 0.8)
mtext("Normal tails", side = 2, line = 1.8, adj = -0.22, col = "grey55",
  cex = 0.8)

with(gomp_hat_base_corr,
  make_panel(x = lambda_50, xl = lambda_25, xu = lambda_75,
    label = paste("Gompertz ", lambda)))

with(gomp_hat_base_corr,
  make_panel(x = b_50, xl = b_25, xu = b_75, label = paste("Gompertz ", b)))

add_sil <- function(x, y, file, tax_class, width_mult = 1, height_mult = 1) {
  #x <- x+5
  y <- y + 0.1
  bg <- subset(cols_df, taxonomic_class == tax_class)$col
  library("grImport")
  par(xpd = NA)
  p <- readPicture(paste0("silhouettes/", file, ".eps.xml"))
  p@paths$path@rgb <- "grey37"
  if(file == "chinook")
    p@paths[[2]]@rgb <- "grey37"  # different vector structure
  grImport::picture(p, x, y, x + 0.28*width_mult, y+0.22*height_mult)
  rect(x-0.2, y+0.06, x-0.2 + 0.14, y + 0.13, col = bg, lwd = 0.8, border = "grey50")
  text(x-0.315, y-0.05, tax_class, col = "grey40", cex = 0.9, pos=4)
  par(xpd = FALSE)
}

add_sil(1.65, 0.3, "aves", tax_class = "Aves")
add_sil(1.65, -0.1, "rabbit", tax_class = "Mammalia")
add_sil(1.65, -0.5, "fly", tax_class = "Insecta")
add_sil(1.65, -0.9, "chinook", tax_class = "Osteichthyes", width_mult = 1.6,
  height_mult = 0.8)

with(gomp_hat_base_corr,
  make_panel(x = dataset_length, log = "x", yaxis = TRUE,
    label = "Time-series length",
    xaxis_ticks = c(20, 30, 50, 100)))
mtext("Heavy tails", side = 2, line = 1.8, adj = 1.15, col = "grey55",
  cex = 0.8)
mtext("Normal tails", side = 2, line = 1.8, adj = -0.22, col = "grey55",
  cex = 0.8)

with(gomp_hat_base_corr, make_panel(x = Len, log = "x",
  label = "Body length (mm)",
  xaxis_ticks = c(2, 5, 20, 100, 500)))
with(gomp_hat_base_corr, make_panel(x = Lifesp/12, log = "x",
  label = "Lifespan (years)", xaxis_ticks = c(0.2, 1, 5, 20)))

# mtext("Possible covariate value", side = 1, outer = TRUE, line = 0,
#   cex = 0.9)
mtext(quote(Pr(nu<10)), side = 2,
  outer = TRUE, line = 2.8, cex = 0.9, adj = 0.6)

dev.off()
