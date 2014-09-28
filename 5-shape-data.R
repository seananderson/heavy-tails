# this file reads in the posterior data, filter it a bit, and joins it with
# the original gpdd and mammal life history databases
# it then does some reshaping for plotting purposes

library(dplyr)

gpdd <- readRDS("gpdd-clean.rds")
gomp_hat_base <- readRDS("gomp_hat_base.rds")
gomp_hat_logistic <- readRDS("gomp_hat_logistic-5.0.rds")
gomp_hat_ar1 <- readRDS("gomp_hat_ar1-5.0.rds")
#gomp_hat_ar1_obs <- readRDS("gomp_hat_ar1_obs0.2.rds")
mammals <- readRDS("mammals.rds")
brook <- read.csv("brook-etal.csv", stringsAsFactors = FALSE)

gomp_hat_ar1 <- subset(gomp_hat_ar1,
  main_id %in% gpdd$main_id)
gomp_hat_base <- subset(gomp_hat_base,
  main_id %in% gpdd$main_id)
#gomp_hat_ar1_obs <- subset(gomp_hat_ar1_obs,
  #main_id %in% gpdd$main_id)

gpdd$dataset_length <- NULL
gomp_hat_base$dataset_length <- NULL
gpdd.temp <- gpdd %>% group_by(main_id) %>% summarise(dataset_length = n())
gomp_hat_base <- plyr::join(gomp_hat_base, gpdd.temp)
gpdd <- plyr::join(gpdd, gpdd.temp)

# back to reshaping
# add pgr to gpdd:
gpdd <- gpdd %>% group_by(main_id) %>%
  mutate(pgr = c(diff(log(population_untransformed)), NA))

# make some joings for looking at the effect of observation error
# and estimating AR1 parameter:
temp <- gomp_hat_base
names(temp) <- paste0(names(temp), "_base")
temp <- plyr::rename(temp, c("main_id_base" = "main_id"))
ar1_vs_base <- plyr::join(gomp_hat_ar1, temp)
rm(temp)
#
# temp <- gomp_hat_ar1
# names(temp) <- paste0(names(temp), "_ar1_no_obs")
# temp <- plyr::rename(temp, c("main_id_ar1_no_obs" = "main_id"))
# ar1_no_obs_vs_obs <- plyr::join(gomp_hat_ar1_obs, temp)
# rm(temp)

# join in taxonomic data from the gpdd:
lookup <- gpdd[,c("main_id", "common_name", "taxonomic_class", "taxon_name",
  "taxonomic_order", "taxonomic_family", "dataset_length")]
lookup <- lookup[!duplicated(lookup), ]

gomp_hat_base <- plyr::join(gomp_hat_base, lookup)
gomp_hat_ar1 <- plyr::join(gomp_hat_ar1, lookup)
ar1_vs_base <- plyr::join(ar1_vs_base, lookup)
#ar1_no_obs_vs_obs <- plyr::join(ar1_no_obs_vs_obs, lookup)

# for the paper:
(max_rhat <- max(gomp_hat_base$max_rhat))
(max_nu_rhat <- max(gomp_hat_base$nu_rhat))
(min_neff <- min(gomp_hat_base$min_neff))
(min_nu_neff <- min(gomp_hat_base$nu_neff))

(max_rhat <- max(gomp_hat_ar1$max_rhat))
(max_nu_rhat <- max(gomp_hat_ar1$nu_rhat))
(min_neff <- min(gomp_hat_ar1$min_neff))
(min_nu_neff <- min(gomp_hat_ar1$nu_neff))
#
# (max_rhat <- max(gomp_hat_ar1_obs$max_rhat))
# (max_nu_rhat <- max(gomp_hat_ar1_obs$nu_rhat))
# (min_neff <- min(gomp_hat_ar1_obs$min_neff))
# (min_nu_neff <- min(gomp_hat_ar1_obs$nu_neff))

# sort nicely for plotting:
gomp_hat_ar1 <- gomp_hat_ar1 %>%
  arrange(taxonomic_class, nu_50) %>%
  group_by(taxonomic_class) %>%
  mutate(sort_id = seq_along(main_id),
    sort_id_perc = 100 * (sort_id / max(sort_id)))

gomp_hat_base <- gomp_hat_base %>%
  arrange(taxonomic_class, nu_50) %>%
  group_by(taxonomic_class) %>%
  mutate(sort_id = seq_along(main_id),
    sort_id_perc = 100 * (sort_id / max(sort_id)))



# bring in PanTHERIA:
mammals <- plyr::rename(mammals, c("binomial" = "taxon_name"))

gomp_hat_ar1 <- plyr::join(gomp_hat_ar1, mammals[,c("taxon_name", "adult_body_mass_g", "adult_head_body_len_mm", "sexual_maturity_age_d", "home_range__indiv_km2", "population_density_n.km2", "trophic_level")])

gomp_hat_ar1 <- mutate(gomp_hat_ar1, log_adult_body_mass_g = log(adult_body_mass_g), log_sexual_maturity_age_d = log(sexual_maturity_age_d), log_sigma_proc_50 = log(sigma_proc_50))

gomp_hat_base <- plyr::join(gomp_hat_base, mammals[,c("taxon_name", "adult_body_mass_g", "adult_head_body_len_mm", "sexual_maturity_age_d", "home_range__indiv_km2", "population_density_n.km2", "trophic_level")])

gomp_hat_base <- mutate(gomp_hat_base, log_adult_body_mass_g = log(adult_body_mass_g), log_sexual_maturity_age_d = log(sexual_maturity_age_d), log_sigma_proc_50 = log(sigma_proc_50))

# bring in Brook et al. 2006 dataset of characteristics
brook <- transform(brook, taxon_name = paste(Genus, Species))
gomp_hat_base <- plyr::join(gomp_hat_base, brook[,c("taxon_name", "Mass",
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
gomp_hat_ar1_long <- melt(gomp_hat_ar1, id.vars = c("main_id", "taxonomic_class", "taxonomic_order"), measure.vars = c("phi_50", "log_sigma_proc_50", "dataset_length", "b_50", "lambda_50", "log_adult_body_mass_g"))

gomp_hat_base_long <- melt(gomp_hat_base, id.vars = c("main_id", "taxonomic_class", "taxonomic_order", "taxon_name"), measure.vars = c("log_sigma_proc_50", "dataset_length", "b_50", "lambda_50", "log_adult_body_mass_g"))

gomp_hat_base_long_brooks <- melt(gomp_hat_base, id.vars = c("main_id", "taxonomic_class", "taxonomic_order", "taxon_name"), measure.vars = c("log_Mass", "log_Len", "MinAge", "log10_MinAge", "Lifesp", "BS", "EF", "HI", "RA"))
gomp_hat_base_long_brooks$variable <-
  plyr::revalue(gomp_hat_base_long_brooks$variable,
  c("BS" = "Body size index",
    "EF" = "Ecological flexibility index",
    "HI" = "Human impact index",
    "RA" = "Geographic range index"))


gomp_hat_ar1_long_p <- plyr::join(gomp_hat_ar1_long, gomp_hat_ar1[,c("main_id", "p20", "p10")])
gomp_hat_base_long_p <- plyr::join(gomp_hat_base_long, gomp_hat_ar1[,c("main_id", "p20", "p10")])

gomp_hat_ar1_long <- plyr::join(gomp_hat_ar1_long, gomp_hat_ar1[,c("main_id", "nu_50", "nu_25", "nu_75")])
gomp_hat_base_long <- plyr::join(gomp_hat_base_long, gomp_hat_base[,c("main_id", "nu_50", "nu_25", "nu_75")])
gomp_hat_base_long_brooks <- plyr::join(gomp_hat_base_long_brooks, gomp_hat_base[,c("main_id", "nu_50", "nu_25", "nu_75")])

# reorder gomp_hat_ar1 for taxonomic order dot plotting:
temp <- gomp_hat_ar1 %>% group_by(taxonomic_order) %>%
  summarise(median_p20 = median(p20)) %>%
  arrange(desc(median_p20)) %>%
  mutate(order_p20 = 1:length(median_p20))

refs <- gpdd[,c("main_id", "ref")]
refs <- refs[!duplicated(refs), ]

gomp_hat_ar1 <- plyr::join(gomp_hat_ar1, temp)
gomp_hat_ar1 <- plyr::join(gomp_hat_ar1, refs, type = "left")
rm(temp)


# clean_heavy <- clean_heavy %>%
#   arrange(taxonomic_class, nu_50) %>%
#   group_by(taxonomic_class) %>%
#   mutate(sort_id = seq_along(main_id),
#     sort_id_perc = 100 * (sort_id / max(sort_id)))


# bring together the estimates from the Ricker-logistic and the Gompertz for
# comparison:
names(gomp_hat_logistic) <- paste0("logistic_", names(gomp_hat_logistic))
gomp_hat_logistic <- plyr::rename(gomp_hat_logistic, c("logistic_main_id" = "main_id"))
gomp_vs_logistic <- plyr::join(gomp_hat_base, gomp_hat_logistic)

