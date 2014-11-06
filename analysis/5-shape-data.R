# this file reads in the posterior data, filter it a bit, and joins it with
# the original gpdd and mammal life history databases
# it then does some reshaping for plotting purposes

library(dplyr)

gpdd <- readRDS("gpdd-clean.rds")

gomp_hat_base <- readRDS("gomp-base-hat.rds")
gomp_hat_logistic <- readRDS("logistic-hat.rds")
gomp_hat_ar1 <- readRDS("gomp-ar1-hat.rds")
gomp_hat_obs_0.2 <- readRDS("gomp-obs-0.2-hat.rds")
gomp_hat_rate <- readRDS("rate-hat.rds")

gomp_hat_weaker <- readRDS("gomp-base-weaker-hat.rds")

brook <- read.csv("brook-etal.csv", stringsAsFactors = FALSE)

gpdd$dataset_length <- NULL
gpdd.temp <- gpdd %>% group_by(main_id) %>% summarise(dataset_length = n())
gpdd <- left_join(gpdd, gpdd.temp)

# back to reshaping
# add pgr to gpdd:
gpdd <- gpdd %>% group_by(main_id) %>%
  mutate(pgr = c(diff(log(population_untransformed)), NA))

# make some joins for looking at the effect of observation error
# and estimating AR1 parameter:
temp <- gomp_hat_base
names(temp) <- paste0(names(temp), "_base")
temp <- plyr::rename(temp, c("main_id_base" = "main_id"))
ar1_vs_base <- left_join(gomp_hat_ar1, temp)
rm(temp)

# join in taxonomic data from the gpdd:
lookup <- gpdd[,c("main_id", "common_name", "taxonomic_class",
  "taxon_name", "taxonomic_order", "taxonomic_family", "dataset_length")]
lookup <- lookup[!duplicated(lookup), ]

gomp_hat_base <- inner_join(gomp_hat_base, lookup, by = "main_id")
gomp_hat_ar1 <- inner_join(gomp_hat_ar1, lookup, by = "main_id")
gomp_hat_obs_0.2 <- inner_join(gomp_hat_obs_0.2, lookup, by = "main_id")
gomp_hat_logistic <- inner_join(gomp_hat_logistic, lookup, by = "main_id")
gomp_hat_weaker <- inner_join(gomp_hat_weaker, lookup, by = "main_id")
ar1_vs_base <- inner_join(ar1_vs_base, lookup, by = "main_id")

# for the paper:
# (max_rhat <- max(gomp_hat_base$max_rhat))
# (max_nu_rhat <- max(gomp_hat_base$nu_rhat))
# (min_neff <- min(gomp_hat_base$min_neff))
# (min_nu_neff <- min(gomp_hat_base$nu_neff))
#
# (max_rhat <- max(gomp_hat_ar1$max_rhat))
# (max_nu_rhat <- max(gomp_hat_ar1$nu_rhat))
# (min_neff <- min(gomp_hat_ar1$min_neff))
# (min_nu_neff <- min(gomp_hat_ar1$nu_neff))

gomp_hat_base <- gomp_hat_base %>%
  arrange(taxonomic_class, nu_50) %>%
  group_by(taxonomic_class) %>%
  mutate(sort_id = seq_along(main_id),
    sort_id_perc = 100 * (sort_id / max(sort_id)))

gomp_hat_base <- mutate(gomp_hat_base, log_sigma_proc_50 = log(sigma_proc_50))

# bring in Brook et al. 2006 dataset of characteristics
brook <- transform(brook, taxon_name = paste(Genus, Species))
gomp_hat_base <- left_join(gomp_hat_base, brook[,c("taxon_name", "Mass",
    "Len", "MinAge", "Lifesp", "BS", "EF", "HI", "RA")])
gomp_hat_base$log_Mass <- log(gomp_hat_base$Mass)
# gomp_hat_base$Mass <- NULL
gomp_hat_base$log_Len <- log(gomp_hat_base$Len)
gomp_hat_base$log10_Len <- log10(gomp_hat_base$Len)
# gomp_hat_base$Len <- NULL
gomp_hat_base$log10_MinAge <- log10(gomp_hat_base$MinAge)
gomp_hat_base$log10_Lifesp <- log10(gomp_hat_base$Lifesp)

# make a long data version for ggplot:
library(reshape2)
gomp_hat_base_long <- melt(gomp_hat_base,
  id.vars = c("main_id", "taxonomic_class", "taxonomic_order", "taxon_name"),
  measure.vars = c("log_sigma_proc_50", "dataset_length", "b_50", "lambda_50"))

gomp_hat_base_long_brooks <- melt(gomp_hat_base,
  id.vars = c("main_id", "taxonomic_class", "taxonomic_order", "taxon_name"),
  measure.vars = c("log_Mass", "log_Len", "MinAge", "log10_MinAge", "Lifesp",
    "BS", "EF", "HI", "RA"))
gomp_hat_base_long_brooks$variable <-
  plyr::revalue(gomp_hat_base_long_brooks$variable,
  c("BS" = "Body size index",
    "EF" = "Ecological flexibility index",
    "HI" = "Human impact index",
    "RA" = "Geographic range index"))

gomp_hat_base_long_p <- left_join(gomp_hat_base_long,
  gomp_hat_ar1[,c("main_id", "p20", "p10")])

gomp_hat_base <- as.data.frame(gomp_hat_base)
gomp_hat_base_long <- left_join(gomp_hat_base_long,
  gomp_hat_base[,c("main_id", "nu_50", "nu_25", "nu_75")])
gomp_hat_base_long_brooks <- left_join(gomp_hat_base_long_brooks,
  gomp_hat_base[,c("main_id", "nu_50", "nu_25", "nu_75")])

# reorder gomp_hat_ar1 for taxonomic order dot plotting:
temp <- gomp_hat_base %>% group_by(taxonomic_order) %>%
  summarise(median_p20 = median(p20)) %>%
  arrange(desc(median_p20)) %>%
  mutate(order_p20 = 1:length(median_p20))

refs <- gpdd[,c("main_id", "ref")]
refs <- refs[!duplicated(refs), ]

gomp_hat_ar1 <- left_join(gomp_hat_ar1, temp)
gomp_hat_ar1 <- left_join(gomp_hat_ar1, refs, type = "left")
rm(temp)
