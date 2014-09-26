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


temp <- gpdd %>% group_by(main_id) %>% summarise(assume_log10 = assumed_log10[1]) %>% summarise(total_assumed_log10 = sum(assume_log10))
total_assumed_log10 <- temp$total_assumed_log10

zz <- file("values.tex", "w")
write_tex(perc_imputed_pops, "percImputedPops")
write_tex(perc_imputed_points, "percImputedPoints")
write_tex(fig2_pop_n, "nuCoefPopN")
write_tex(total_assumed_log10, "totalAssumedLog")
close(zz)
