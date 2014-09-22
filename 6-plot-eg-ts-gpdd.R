heavy <- filter(gomp_hat_base, p10 > .5) %>%
#   filter(fur == FALSE) %>%
  select(main_id, common_name, p10, p20, nu_5, nu_25, nu_50, nu_75, nu_95)

#gpdd_heavy <- filter(gpdd, main_id %in% heavy$main_id)
#gpdd_heavy <- filter(gpdd_heavy, !grepl("harvest", ref))
#gpdd_heavy <- filter(gpdd_heavy, !grepl("commerce", ref))


heavy <-
  plyr::join(heavy, gpdd) %>%
  arrange(main_id, p10) %>%
  mutate(id_name = paste(round(p10, 2), common_name, main_id))

library(ggplot2)
p <- ggplot(heavy, aes(series_step, population_untransformed, colour = taxonomic_class)) + geom_line() + facet_wrap(~id_name, scales = "free")
ggsave("ts-gpdd-heavy-eg-base.pdf", width = 15, height = 12)

p <- ggplot(heavy, aes(series_step, log10(population_untransformed), colour = taxonomic_class)) + geom_line() + facet_wrap(~id_name, scales = "free", ncol = 4)
ggsave("ts-gpdd-heavy-eg-log10-base.pdf", width = 15, height = 17)

heavy_table <- plyr::ddply(heavy, "main_id", function(x) x[1,]) %>%
  arrange(id_name)

write.csv(heavy_table, file = "heavy.csv", row.names = FALSE)
