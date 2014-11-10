# massive time-series plot:
library("ggplot2")

make_big_ts_plot <- function(dat_, id, width, height) {
  p <- ggplot(dat_, aes(series_step, log10(population_untransformed))) +
      #colour = taxonomic_class)) +
  geom_point(pch = 20, col = "grey70") + geom_line() +
  facet_wrap(~label, scales = "free") +
  geom_point(data = subset(dat_, interpolated == TRUE),
    aes(series_step, log10(population_untransformed)), colour = "red", pch = 20, cex = 3.5) +
  geom_point(data = subset(dat_, zero_sub == TRUE),
    aes(series_step, log10(population_untransformed)), colour = "blue", pch = 20, cex = 3.5) +
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank()) +
  xlab("Time series step") + ylab("log10(Abundance)")
  ggsave(paste0("all-clean-ts-", id, ".pdf"), width = width, height = height,
    limitsize = FALSE)
}
make_big_ts_plot(gpdd[gpdd$taxonomic_class == "Insecta", ], id = "insects",
  width = 24, height = 20)
make_big_ts_plot(gpdd[gpdd$taxonomic_class == "Mammalia", ], id = "mammals",
  width = 24, height = 20)
make_big_ts_plot(gpdd[gpdd$taxonomic_class == "Aves", ], id = "birds",
  width = 24, height = 20)
make_big_ts_plot(gpdd[gpdd$taxonomic_class %in%
  c("Gastropoda", "Crustacea", "Chondrichtyhes", "Osteichthyes"), ],
  id = "fishes-others", width = 24, height = 20)
