# this file produces some values for inclusion of the text of the associated
# paper
#
library(dplyr)

write_tex <- function(x, macro, ...) {
  out <- paste0("\\newcommand{\\", macro, "}{", x, "}")
  cat(out, file = zz)
  cat("\n", file = zz)
  #writeLines(out, con = "values.tex", ...)
}

# what's the median and mean nu in the exponential prior?
#x <- rexp(1e7, 0.01)
#median(x)
#mean(x)

# how much was imputed?
gpdd <- readRDS("gpdd-clean.rds")
perc_imputed <- gpdd %>% group_by(main_id) %>% summarise(has_imputed = ifelse(max(interpolated) == 1, TRUE, FALSE)) %>% summarise(percent_imputed = sum(has_imputed)/length(has_imputed)) %>% round(2) * 100
perc_imputed_pops <- as.character(perc_imputed$percent_imputed)

perc_imputed_points <- round(sum(gpdd$interpolated)/nrow(gpdd) * 100, 1)

fig2_pop_n <- filter(gpdd, taxonomic_class %in% c("Aves", "Mammalia", "Insecta", "Osteichthyes")) %>% summarise(n = length(unique(main_id))) %>% summarise(N = sum(n))
fig2_pop_n <- fig2_pop_n$N

NPops <- length(unique(gpdd$main_id))
NOrders <- length(unique(gpdd$taxonomic_order))
NClasses <- length(unique(gpdd$taxonomic_class))


# p_base <- round(length(which(gomp_hat_base$p10 > 0.5)) / NPops * 100)
# p_obs <- round(length(which(gomp_hat_obs_0.2$p10 > 0.5)) / NPops * 100)
# p_logistic <- round(length(which(gomp_hat_logistic$p10 > 0.5)) / NPops * 100)
# p_rate <- round(length(which(gomp_hat_rate$p10 > 0.5)) / NPops * 100)
# MinPNu10 <- min(p_base, p_obs, p_logistic)
# MaxPNu10 <- max(p_base, p_obs, p_logistic)

cols_include <- c("taxonomic_class", "main_id", "p10", "nu_50", "type",
  "max_rhat", "min_neff")
gomp_hat_base$ type <- "base"
gomp_hat_ar1$type <- "ar1"
gomp_hat_logistic$type <- "logistic"
gomp_hat_obs_0.2$type <- "obs"

gtemp <- rbind(
  gomp_hat_base[,cols_include],
  gomp_hat_ar1[,cols_include],
  gomp_hat_logistic[,cols_include],
  gomp_hat_obs_0.2[,cols_include])

# gtemp %>% filter(max_rhat < 1.05) %>% group_by(taxonomic_class, type) %>% summarise(n = n(), h =  length(which(p10 > 0.5)), p = round(100 * h / n))

pheavy_class <- gtemp %>% filter(max_rhat < 1.05) %>% group_by(taxonomic_class, type) %>% summarise(n = n(), h =  length(which(p10 > 0.5)), p = round(100 * h / n)) %>% group_by(taxonomic_class) %>% summarise(min_p = min(p), max_p = max(p)) %>% filter(min_p > 0)



pheavy_overall <- gtemp %>% group_by(type) %>% summarise(n = n(), h = length(which(p10 > 0.5)), p = round(100 * h / n))

temp <- gpdd %>% group_by(main_id) %>% summarise(assume_log10 = assumed_log10[1]) %>% summarise(total_assumed_log10 = sum(assume_log10))
total_assumed_log10 <- temp$total_assumed_log10

zz <- file("values.tex", "w")
write_tex(perc_imputed_pops, "percImputedPops")
write_tex(perc_imputed_points, "percImputedPoints")
write_tex(fig2_pop_n, "nuCoefPopN")
write_tex(total_assumed_log10, "totalAssumedLog")
for(cl in unique(pheavy_class$taxonomic_class)) {
  write_tex(paste0(
    pheavy_class[pheavy_class$taxonomic_class == cl, "min_p"],
    "--",
    pheavy_class[pheavy_class$taxonomic_class == cl, "max_p"]),
    paste0(cl, "RangePerc"))
}

write_tex(min(pheavy_overall$p), "overallMinPerc")
write_tex(max(pheavy_overall$p), "overallMaxPerc")
write_tex(NPops, "NPops")
write_tex(NOrders, "NOrders")
write_tex(NClasses, "NClasses")
close(zz)
