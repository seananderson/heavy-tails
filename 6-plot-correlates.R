# a pretty version of the correlates of nu plot for the main paper

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
  yl <- gomp_hat_base_corr$nu_25
  yu <- gomp_hat_base_corr$nu_75
  y <- gomp_hat_base_corr$nu_50
  class_col <- gomp_hat_base_corr$col

  y <- y[!is.na(x)]
  yl <- yl[!is.na(x)]
  yu <- yu[!is.na(x)]
  class_col <- class_col[!is.na(x)]
  x <- x[!is.na(x)]

  plot(1, 1, xlim = range(x), ylim = c(0, 0.5), type = "n", xlab = "",
    ylab = "", log = log, yaxs = "i", axes = FALSE)

  if(!is.null(xl) & !is.null(xu))
    segments(xl, 1/y, xu, 1/y, col = paste0(class_col, "50"), lwd = 0.9)

  segments(x, 1/yl, x, 1/yu, col = paste0(class_col, "50"), lwd = 0.9)
  points(x, 1/y, pch = 21, col = "#66666690", lwd = 0.8,
    bg = paste0(class_col, ""), cex = 1)
  if(is.null(xaxis_ticks))
    axis(1)
  else
    axis(1, at = xaxis_ticks)
  if(yaxis)
    axis(2, las = 1, at = 1/c(2, 3, 5, 10, Inf), labels = c(2, 3, 5, 10, Inf))
  box()
  par(xpd = NA)
  lab <- letters[i]
  legend("topleft", legend = substitute(paste("(", lab, ") ", label)),
    bty = "n", inset = c(-0.1, -0.04), cex = 0.9)
  par(xpd = FALSE)
}

pdf("correlates.pdf", width = 7, height = 4.55)
par(mfrow = c(2, 3), mar = c(2.5,0,0,0), oma = c(.8, 4.2,1.1, 5.4),
  tck = -0.02, mgp = c(2, 0.5, 0), col.axis = "grey25", col = "grey25")
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
# text(-0.5, 1/10, "Density\ndependent", cex = 0.8, col = "grey50")

add_sil <- function(x, y, file, tax_class, width_mult = 1, height_mult = 1) {
  #x <- x+5
  y <- y + -0.1
  bg <- subset(cols_df, taxonomic_class == tax_class)$col
  library(grImport)
  par(xpd = NA)
  p <- readPicture(paste0("silhouettes/", file, ".eps.xml"))
  p@paths$path@rgb <- "grey37"
  if(file == "chinook")
    p@paths[[2]]@rgb <- "grey37"  # different vector structure
  grImport::picture(p, x, y, x + 0.28*width_mult, y+0.1*height_mult)
  rect(x-0.2, y+0.05, x-0.2 + 0.14, y + 0.08, col = bg, lwd = 0.8, border = "grey50")
  text(x-0.315, y-0.025, tax_class, col = "grey40", cex = 0.9, pos=4)
  par(xpd = FALSE)
}

#add_sil(-12, 0.54, "grey-heron", tax_class = "Aves")
#add_sil(-2, 0.54, "rabbit", tax_class = "Mammalia")
#add_sil(8, 0.54, "fly", tax_class = "Insecta")
#add_sil(17.5, 0.54, "chinook", tax_class = "Osteichthyes", width_mult = 1.7, height_mult = 0.8)
add_sil(1.65, 0.3, "grey-heron", tax_class = "Aves")
add_sil(1.65, 0.1, "rabbit", tax_class = "Mammalia")
add_sil(1.65, -0.1, "fly", tax_class = "Insecta")
add_sil(1.65, -0.3, "chinook", tax_class = "Osteichthyes", width_mult = 1.7, height_mult = 0.8)

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

mtext("Possible correlate value", side = 1, outer = TRUE, line = -0.5,
  cex = 0.9)
mtext(quote(t~distribution~degrees~of~freedom~(widehat(nu))), side = 2,
  outer = TRUE, line = 2.8, cex = 0.9, adj = 0.6)


dev.off()
