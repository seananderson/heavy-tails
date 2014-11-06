# Need to run 5.9-order-level-posteriors.R first

# a_df_select <- filter(a_df, n_pops > 3)

source("add_phylopic.R") # modifed from Scott's rphylopic package

class_cols <- lu[, c("taxonomic_class", "col")]
class_cols <- class_cols[!duplicated(class_cols), ]

a_df <- plyr::join(a_df, class_cols, by = "taxonomic_class")

# get silhouette images:
or <- read.table("orders.csv", stringsAsFactors = F, header = T)
or <- mutate(or,
  hash = gsub("http://phylopic.org/image/([0-9a-z-]+)/", "\\1", url),
  svg_url = paste0("http://phylopic.org/assets/images/submissions/",
    hash, ".svg"),
  png_url = paste0("http://phylopic.org/assets/images/submissions/",
    hash, ".256.png"),
  scaling_factor = rep(1, nrow(or)))
or[or$taxonomic_order == "Gadiformes", "scaling_factor"] <- 0.6
or[or$taxonomic_order == "Salmoniformes", "scaling_factor"] <- 0.6
or[or$taxonomic_order == "Perciformes", "scaling_factor"] <- 0.6
or[or$taxonomic_order == "Pleuronectiformes", "scaling_factor"] <- 0.8

if(any(!file.exists(paste0("silhouettes/", or$taxonomic_order, ".png")))) {
  for(i in 1:nrow(or)) {
    download.file(or$png_url[i],
      destfile = paste0("silhouettes/", or$taxonomic_order[i], ".png"))
  }
}

op <- vector(mode = "list", length = nrow(a_df))
for(i in seq_along(op)) {
  op[[i]]$taxonomic_order <- a_df$taxonomic_order[i]
  op[[i]]$n_pops <- a_df$n_pops[i]
  op[[i]]$order_id <- a_df$order_id[i]
  op[[i]]$class_id <- a_df$class_id[i]
  op[[i]]$taxonomic_class <- a_df$taxonomic_class[i]
  op[[i]]$col <- a_df$col[i]
  op[[i]]$img <- png::readPNG(paste0("silhouettes/", or$taxonomic_order[i], ".png"))
}

for(i in seq_along(op)) {
  op[[i]]$post <- plogis(mu_a +
        a_order[, op[[i]]$order_id] +
        a_class[, op[[i]]$class_id])
  op[[i]]$dens <- density(op[[i]]$post,
    from = quantile(op[[i]]$post, probs = 0.001)[[1]],
    to = quantile(op[[i]]$post, probs = 0.999)[[1]])
  op[[i]]$med_post <- median(op[[i]]$post)
  op[[i]]$med_dens_height <-
    op[[i]]$dens$y[max(which(op[[i]]$dens$x < op[[i]]$med_post))]
}

x <- rexp(2e6, 0.01)
x <- x[x > 2]
prior_p10 <- length(x[x < 10])/length(x)

pdf("order-posteriors.pdf", width = 3.5, height = 4.5)
par(mfrow = c(1, 1), mar = c(2.7,8,0,0), oma = c(0.2, 0.2, 1.1, 0.8),
  tck = -0.02, mgp = c(2, 0.5, 0), col.axis = "grey25", col = "grey25")
par(cex = 0.8)

xlim <- c(-.025, 0.29)
plot(1, 1, xlim = xlim, ylim = c(length(op), 1), type = "n",
  ylab = "", xlab = "", axes = FALSE, xaxs = "i")
abline(v = prior_p10, lty = 2, col = "grey40", lwd = 0.6)
scaling_factor <- 55
for(i in seq_along(op)) {
  segments(xlim[1]+0.05, i, min(op[[i]]$dens$x), i, col = "grey90")
  segments(max(op[[i]]$dens$x), i, xlim[2], i, col = "grey90")
  polygon(c(op[[i]]$dens$x, rev(op[[i]]$dens$x)),
    i + c(op[[i]]$dens$y/scaling_factor, -rev(op[[i]]$dens$y/scaling_factor)),
    border = "grey50", lwd = 0.5, col = "white")
  polygon(c(op[[i]]$dens$x, rev(op[[i]]$dens$x)),
    i + c(op[[i]]$dens$y/scaling_factor, -rev(op[[i]]$dens$y/scaling_factor)),
    border = "grey50", lwd = 1.5, col = paste0(op[[i]]$col, "90"))
  segments(
    op[[i]]$med_post,
    i - op[[i]]$med_dens_height/scaling_factor,
    op[[i]]$med_post,
    i + op[[i]]$med_dens_height/scaling_factor,
    col = "grey50", lwd = 1.5)
  par(xpd = NA)
  add_phylopic(op[[i]]$img, alpha = 1, x = 0.004, y = i,
    ysize = 0.9 * or$scaling_factor[i], xy_ratio = 35, color = "grey45")
  par(xpd = FALSE)
}

axis(2, at = seq_along(op),
  labels = as.character(unlist(lapply(op, function(x) x$taxonomic_order))),
  las = 1, lwd = 0, line = -0.6)
axis(1, at = seq(0, 0.3, 0.1))
mtext(quote(Pr(nu<10)), side = 1, line = 2, cex = 0.85)
# box()
dev.off()
