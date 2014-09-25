# a pretty version of nu coeffs and p(nu < 10) for a main paper figure

classes <- c("Aves","Mammalia", "Insecta", "Osteichthyes")
cols <- c(rev(RColorBrewer::brewer.pal(5, "YlOrRd"))[-c(4, 5)])

add_label <- function(xfrac = -0.02, yfrac = -0.04, label = "", pos = 4, ...) {
  u <- par("usr")
  x <- u[1] + xfrac * (u[2] - u[1])
  y <- u[4] - yfrac * (u[4] - u[3])
  text(x, y, label, pos = pos, ...)
}

set.seed(1) # for jittering
pdf("nu-coefs-2.pdf", width = 3.8, height = 9.5)

par(mfrow = c(4, 1), mar = c(1.3,0,0,0), oma = c(2.2, 3.5, 1, .8), tck = -0.02, mgp = c(2, 0.5, 0), col.lab = "grey20", col.axis = "grey20", col = "grey20")
par(cex = 0.9)

for(i in 1:length(classes)) {

xaxt  <- ifelse(i == 4, "s", "n")

x <- subset(gomp_hat_base, taxonomic_class == classes[i])

ticks <- c(2, 3, 5, 10, Inf)


plot(1, 1, type = "n", xlim = c(-0.2, 100.2), ylim = c(0, 0.5), xlab = "", ylab = "", yaxt = "n", xaxs = "i", yaxs = "i", xaxt = xaxt, las = 1)
#axis(2, at = seq(0, 100, 50), las = 1)
par(xpd = NA)
add_label(label = paste0("(", letters[i], ") ", classes[i]))
if(i == 1) {
  mtext("Heavy tails", side = 2, line = 1.8, adj = 0.90, col = "grey55")
  mtext("Normal tails", side = 2, line = 1.8, adj = -0.09, col = "grey55")
}
par(xpd = FALSE)
axis(2, las = 1, at = 1/ticks, labels = ticks)

  #locs <- 1/ticks
  #for(j in 1:length(cols)) {
    #rect(0, locs[j], 100, locs[j+1], border = NA,
      #col = paste0(cols[j], "35"))
  #}

#segments(1/x$nu_5, x$sort_id_perc, 1/x$nu_95, x$sort_id_perc, lwd = 1.1, col = "grey60")
#segments(1/x$nu_25, x$sort_id_perc, 1/x$nu_75, x$sort_id_perc, lwd = 1.1, col = "grey30")
#points(1/x$nu_50, x$sort_id_perc, pch = 21, col = "white", lwd = 0.5, bg = "white", cex = 0.4)

segments(x$sort_id_perc, 1/x$nu_5, x$sort_id_perc, 1/x$nu_95, lwd = 1.1, col = "grey60")
segments(x$sort_id_perc, 1/x$nu_25, x$sort_id_perc, 1/x$nu_75, lwd = 1.1, col = "grey30")
points(x$sort_id_perc, 1/x$nu_50, pch = 21, col = "grey10", lwd = 0.5, bg = "black", cex = 0.6)



 #x <- filter(x, !taxonomic_order %in% c("Cuculiformes", "Gruiformes", "Pteroclidiformes", "Procellariiformes", "Gaviiformes"))
if(unique(x$taxonomic_class) == "Aves") {
  cutoff <- 4
} else {
  cutoff <- 1
}
 x <- x %>% group_by(taxonomic_order) %>%
  mutate(n_pops = length(p10)) %>%
  filter(n_pops >= cutoff)

  #filter(!taxonomic_order %in% c("Cuculiformes", "Gruiformes", "Pteroclidiformes", "Procellariiformes", "Gaviiformes"))

 p <- x %>% group_by(taxonomic_order) %>% summarise(med_p10 = mean(p10)) %>%
  arrange(med_p10) %>% mutate(p10_order = 1:length(med_p10))
 x <- plyr::join(x, p, by = "taxonomic_order")
 TeachingDemos::subplot({
 plot(x$p10, x$p10_order, type = "n", xlab = "", ylab = "", las = 1, bty = "n", yaxt = "n", xaxt = "n", axes = FALSE)
 axis(2, at = p$p10_order, labels = p$taxonomic_order, las = 1, cex = 0.8, cex.axis= 0.8, lwd = 0, col.axis = "grey30")
 axis(1, at = c(0, 0.5, 1), col = "grey30", lwd = 0.7, tck = -0.04, mgp = c(2, 0.25, 0), col.axis = "grey30", cex.axis = 0.8)
 mtext(expression(p(nu<10)), side = 1, line = 1.4, cex = 0.8, col = "grey30")
 abline(h = as.numeric(as.factor(x$taxonomic_order)), col = "grey93", lwd = 0.7, lty = 1)
   par(xpd = NA)
 points(x$p10, jitter(x$p10_order, amount = 0.25), pch = 21, bg = "grey80", col = "grey40", cex = 0.7)
   par(xpd = FALSE)
 },
 x = c(60, 95),
 y = c(0.25, 0.48))
}

par(xpd = NA)
#axis(1, at = 1/ticks, labels = ticks, las = 1)
mtext(quote(t~distribution~degrees~of~freedom~(widehat(nu))), outer = TRUE, side = 2, line = 2.2)
mtext("Percentage of populations", outer = TRUE, side = 1, line = 1)
par(xpd = FALSE)



dev.off()


#
#
# par(mfrow = c(2, 1))
# plot(3, 3)
# subplot({plot(1, 1)}, c(0.5, 3), c(0.5, 3))
# plot(1, 1)
# subplot(plot(1, 1), c(0.5, 1), c(0.5, 1))
