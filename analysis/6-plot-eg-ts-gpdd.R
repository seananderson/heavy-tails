heavy <- filter(gomp_hat_base, p10 > .5) %>%
#   filter(fur == FALSE) %>%
  select(main_id, common_name, p10, p20, nu_5, nu_25, nu_50, nu_75, nu_95)

#gpdd_heavy <- filter(gpdd, main_id %in% heavy$main_id)
#gpdd_heavy <- filter(gpdd_heavy, !grepl("harvest", ref))
#gpdd_heavy <- filter(gpdd_heavy, !grepl("commerce", ref))


heavy <-
  plyr::join(heavy, as.data.frame(gpdd)) %>%
  arrange(main_id, p10) %>%
  mutate(id_name = paste(round(p10, 2), common_name, main_id))

library(ggplot2)
theme_set(theme_bw())

p <- ggplot(heavy, aes(series_step, population_untransformed, colour = taxonomic_class)) + geom_line() + facet_wrap(~id_name, scales = "free")
ggsave("ts-gpdd-heavy-eg-base.pdf", width = 15, height = 12)

p <- ggplot(heavy, aes(decimal_year_begin, log10(population_untransformed), colour = taxonomic_class)) + geom_line() + facet_wrap(~id_name, scales = "free", ncol = 4) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_point() +
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank())
ggsave("ts-gpdd-heavy-eg-log10-base.pdf", width = 15, height = 17)

heavy_table <- plyr::ddply(heavy, "main_id", function(x) x[1,]) %>%
  arrange(id_name)

write.csv(heavy_table, file = "heavy.csv", row.names = FALSE)

p<-ggplot(subset(gpdd, assumed_log10 == TRUE), aes(series_step, log10(population_untransformed))) + geom_point() + geom_line() + facet_wrap(~label, scales = "free_x") + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
ggsave("log10-assumed.pdf", width = 13, height = 13)

zz <- file("mainids.tex", "w")
cat(sort(unique(gpdd$main_id)), file = zz)
close(zz)
