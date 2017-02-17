# Make a table for the SOM with citations and verified (or not) causes

library("dplyr")
library("xtable")

heavy_refs <- read.csv("heavy-filled.csv", stringsAsFactors = FALSE,
  strip.white = TRUE)[, c("citation", "common_name", "known_cause", "reasons",
    "data_correct", "taxon_name", "main_id", "exact_name")]

source("5-shape-data.R")

heavy <- filter(gomp_hat_base, p10 > .6) %>%
  select(main_id, common_name, p10, nu_50, nu_5, nu_95, taxon_name)

heavy <- inner_join(select(heavy, -taxon_name, -common_name), heavy_refs)

heavy <- mutate(heavy,
  population = paste0(common_name, ", \\textit{", taxon_name, "}, ",
    exact_name)) %>%
  select(-common_name, -taxon_name, -exact_name) %>%
  arrange(nu_50) %>%
  mutate(p10 = sprintf("%.2f", round(p10, 2))) %>%
  mutate(nu_50 = sprintf("%.0f", round(nu_50, 0))) %>%
  mutate(nu_5 = sprintf("%.0f", round(nu_5, 0))) %>%
  mutate(nu_95 = sprintf("%.0f", round(nu_95, 0))) %>%
  mutate(spark = paste0("\\includegraphics[width=1.7cm]{../analysis/sparks/", main_id, ".pdf}")) %>%
  mutate(nu_hat = paste0(nu_50, " (", nu_5, "--", nu_95, ")")) %>%
  select(spark, population, main_id, citation, reasons, p10) %>%
  rename("Time series" = spark, "Population" = population, "ID" = main_id,
    "Ref" = citation,
    "Description" = reasons, "Pr($\\nu < 10$)" = p10)#, "$\\widehat{\\nu}$" = nu_hat)

heavy$Ref <- sub("citep", "cite", heavy$Ref)
heavy$Ref <- sub("citet", "cite", heavy$Ref)

print.xtable(xtable(heavy,
    caption = ""),
  include.rownames = FALSE, file = "cause-table.tex",
  booktabs = TRUE,  caption.placement = "top", size = "footnotesize",
  sanitize.text.function = identity, only.contents = TRUE, timestamp = NULL)
