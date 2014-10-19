# system("extract_gpdd.sh")

dat <- read.csv("data.csv", stringsAsFactors = FALSE,
  strip.white = TRUE)
taxon <- read.csv("taxon.csv", stringsAsFactors = FALSE,
  strip.white = TRUE)
datasource <- read.csv("datasource.csv", stringsAsFactors = FALSE,
  strip.white = TRUE)
version <- read.csv("version.csv", stringsAsFactors = FALSE,
  strip.white = TRUE)
timeperiod <- read.csv("timeperiod.csv", stringsAsFactors = FALSE,
  strip.white = TRUE)
restricted_datasets <- read.csv("restricted_datasets.csv",
  stringsAsFactors = FALSE, strip.white = TRUE)
location <- read.csv("location.csv", stringsAsFactors = FALSE,
  strip.white = TRUE)
biotope <- read.csv("biotope.csv", stringsAsFactors = FALSE,
  strip.white = TRUE)
main <- read.csv("main.csv", stringsAsFactors = FALSE,
  strip.white = TRUE)

underscore_names <- function(x) {
  y <- gsub("([a-z])([A-Z])", "\\1_\\2", names(x))
  tolower(y)
}

names(dat) <- underscore_names(dat)
names(taxon) <- underscore_names(taxon)
names(datasource) <- underscore_names(datasource)
names(version) <- underscore_names(version)
names(timeperiod) <- underscore_names(timeperiod)
names(restricted_datasets) <- underscore_names(restricted_datasets)
names(location) <- underscore_names(location)
names(biotope) <- underscore_names(biotope)
names(main) <- underscore_names(main)

dat$decimal_year_begin <- as.numeric(dat$decimal_year_begin)
dat$decimal_year_end <- as.numeric(dat$decimal_year_end)
dat$population <- as.numeric(dat$population)
dat$population_untransformed <- as.numeric(dat$population_untransformed)

location$longitude_degrees <- as.numeric(location$longitude_degrees)
location$longitude_minutes <- as.numeric(location$longitude_minutes)
location$latitude_degrees <- as.numeric(location$latitude_degrees)
location$latitude_minutes <- as.numeric(location$latitude_minutes)
location$lat_dd <- as.numeric(location$lat_dd)
location$long_dd <- as.numeric(location$long_dd)
location$area <- as.numeric(location$area)

main$dataset_length <- as.numeric(main$dataset_length)
main$sibly_fitted_theta <- as.numeric(main$sibly_fitted_theta)
main$sibly_theta_cilower <- as.numeric(main$sibly_theta_cilower)
main$sibly_theta_ciupper <- as.numeric(main$sibly_theta_ciupper)
main$reliability <- as.numeric(main$reliability)
main$sampling_frequency <- as.numeric(main$sampling_frequency)
main$sibly_return_rate <- as.numeric(main$sibly_return_rate)
main$sibly_carrying_capacity <- as.numeric(main$sibly_carrying_capacity)
main <- plyr::rename(main, c("data_source_id" = "datasource_id"))

datasource <- dplyr::mutate(datasource, ref = paste(author, year, title, reference),
  data_notes = notes)

gpd <- plyr::join(dat, main[,c("main_id", "taxon_id", "location_id",
  "sampling_units", "sampling_protocol", "sampling_effort", "reliability",
  "sampling_frequency", "dataset_length", "notes", "datasource_id")], by = "main_id")

gpd <- plyr::join(gpd, taxon[,c("taxon_id", "taxon_name", "common_name",
  "taxonomic_phylum", "taxonomic_class", "taxonomic_order",
  "taxonomic_family")], by = "taxon_id")

gpd <- plyr::join(gpd, location[,c("location_id", "exact_name",
  "long_dd", "lat_dd")], by = "location_id")

gpd <- plyr::join(gpd, datasource[,c("datasource_id", "ref", "data_medium",
    "data_notes")], by = "datasource_id")

saveRDS(gpd, file = "../gpdd.rds")
