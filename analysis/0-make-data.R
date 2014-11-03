# - this file makes the filtered GPDD dataset
# - some of the gpdd processing happens in another folder (`gpdd`)
#   but gets sourced from here
# - note that there's a shell script in the `gpdd` folder to turn the original
#   Access database files into the `.csv` files
# - I'm not running it here because I haven't put the appropriate binary in R's
#   path

library(dplyr)

# GPDD:
setwd("gpdd")
source("make-gpdd.R")
setwd("..")

# Data filtering (TODO note)
gpdd <- readRDS("gpdd.rds")
#
# for checking later:
gpdd$population_untransformed_original <- gpdd$population_untransformed

# There are a couple that are out of order
# this tripped me up before:
gpdd <- arrange(gpdd, main_id, series_step)

gpdd <- subset(gpdd, sampling_protocol != "Harvest")

# hidden harvests: the fur trade stats:
gpdd <- plyr::ddply(gpdd, "main_id", function(y) {
  fur <- grepl("harvest", y$ref[1])
  if(!fur) y
})
gpdd <- plyr::ddply(gpdd, "main_id", function(y) {
  fur <- grepl("commerce", y$ref[1])
  if(!fur) y
})

# now a big one: some have missing time steps
# first go through and expand missing rows to NAs:
gpdd <- plyr::ddply(gpdd, "main_id", function(y) {
  steps_df <- data.frame(series_step =
      min(y$series_step, na.rm = TRUE):max(y$series_step, na.rm = TRUE),
    stringsAsFactors = FALSE, main_id = y$main_id[1])
  out <- plyr::join(steps_df, y, by = c("main_id", "series_step"))
  out <- dplyr::arrange(out, series_step)
  as.data.frame(out)
})

# now deal with hidden logged values
# assume those with < 0 are logged:
# (some say, some don't)
gpdd <- plyr::ddply(gpdd, "main_id", function(x) {
  x$assumed_log10 <- FALSE
  if(min(x$population_untransformed, na.rm = TRUE) < 0) {
    x$population_untransformed <- 10^x$population_untransformed
    x$assumed_log10 <- TRUE
  }
  x
})

# substitute single zeros for the lowest value
# as in Brook et al. 2006
# (slightly different since they did the first 0 no matter what)
# I'm only subbing out zeros if they are surrounded by non-zeros

gpdd <-
  plyr::ddply(gpdd, "main_id", function(y) {
  y$zero_sub <- FALSE
  for(i in 2:(nrow(y)-1)) {
    if(!is.na(y$population_untransformed[i]) &
        !is.na(y$population_untransformed[i-1]) &
        !is.na(y$population_untransformed[i+1])) {
      if(y$population_untransformed[i] == 0 &
          y$population_untransformed[i-1] != 0 &
          y$population_untransformed[i+1] != 0) {
        lowest_non_zero <-
          min(y$population_untransformed[-which(y$population_untransformed == 0)],
            na.rm = TRUE)
        y$population_untransformed[i] <- lowest_non_zero
        y$decimal_year_begin[i] <- mean(c(y$decimal_year_begin[i-1],
          y$decimal_year_begin[i + 1]))
        y$decimal_year_end[i] <- mean(c(y$decimal_year_end[i-1],
          y$decimal_year_end[i + 1]))
        y$sample_year[i] <- mean(c(y$sample_year[i-1], y$sample_year[i + 1]))
        copy_cols <- c("time_period_id", "taxon_id", "location_id",
          "sampling_units", "sampling_protocol", "sampling_effort",
          "sampling_frequency", "dataset_length", "notes", "datasource_id",
          "taxon_name", "common_name", "taxonomic_phylum", "taxonomic_class",
          "taxonomic_order", "taxonomic_family", "exact_name", "long_dd",
          "lat_dd", "ref", "data_medium", "data_notes")
        y[i, copy_cols] <- y[i-1, copy_cols]
         y$zero_sub[i] <- TRUE
      }
    }
  }
  y
})

# now turn the remaining zeros to NAs for our windowing function
# below
gpdd$population_untransformed[gpdd$population_untransformed == 0] <- NA

# find those with uneven sampling intervals:
multi_sample_diff_ids <-
  gpdd %>% group_by(main_id) %>%
  mutate(sample_diff = decimal_year_end - decimal_year_begin) %>%
  group_by(main_id) %>%
  summarize(n_sample_diff = length(as.numeric(na.omit(unique(sample_diff))))) %>%
  filter(n_sample_diff > 1)
gpdd <- subset(gpdd, !main_id %in% multi_sample_diff_ids$main_id)

saveRDS(gpdd, file = "gpdd-temp.rds") # save time later if needed

# remove those with 4 or more identical values in a row:
check_identical_window <- function(x) {
  ident_window <- FALSE
  x <- as.numeric(na.omit(x))
  if(length(x) > 5) {
    for(i in 4:length(x)) {
      if(length(unique(x[(i-3):i])) == 1)
        ident_window <- TRUE
    }
  } else {
    ident_window <- TRUE # too short to check, so remove now
  }
  ident_window
}

gpdd <- gpdd %>% group_by(main_id) %>%
  mutate(identical_window = check_identical_window(population)) %>%
  filter(identical_window == FALSE)
gpdd$identical_window <- NULL

# as in Brook et al. 2006, remove those without at least 4 unique values
gpdd <- plyr::ddply(gpdd, "main_id", function(x) {
  if(length(as.numeric(na.omit(unique(x$population_untransformed)))) >= 4)
    x
})

# do the interpolation:
gpdd <- plyr::ddply(gpdd, "main_id", function(y) {
  y$interpolated <- FALSE
  interpolated_count <- 0
  for(i in 2:(nrow(y)-1)) {
    if(is.na(y$population_untransformed[i])
      & !is.na(y$population_untransformed[i-1])
      & !is.na(y$population_untransformed[i+1])) {
      y$population_untransformed[i] <- exp(mean(c(log(y$population_untransformed[i-1]),
        log(y$population_untransformed[i + 1]))))
      y$decimal_year_begin[i] <- mean(c(y$decimal_year_begin[i-1],
          y$decimal_year_begin[i + 1]))
      y$decimal_year_end[i] <- mean(c(y$decimal_year_end[i-1],
          y$decimal_year_end[i + 1]))
      y$sample_year[i] <- mean(c(y$sample_year[i-1], y$sample_year[i + 1]))
      interpolated_count <- interpolated_count + 1
      copy_cols <- c("time_period_id", "taxon_id", "location_id",
        "sampling_units", "sampling_protocol", "sampling_effort",
        "sampling_frequency", "dataset_length", "notes", "datasource_id",
        "taxon_name", "common_name", "taxonomic_phylum", "taxonomic_class",
        "taxonomic_order", "taxonomic_family", "exact_name", "long_dd",
        "lat_dd", "ref", "data_medium", "data_notes")
      y[i, copy_cols] <- y[i-1, copy_cols]
      y$interpolated[i] <- TRUE
    }
  }
  return(y)
  })


# now go through and for those that still have NAs, take the longest unbroken
# window if it's 25 or over steps
#
# For testing:
#  junk <- data.frame(main_id = rep("a", 11),
#    population = c(4, 4, 3, 2, 3, 3, 4, 4, 3, 4, 4))
#
#  plyr::ddply(junk, "main_id", function(y) {
#    rl <- rle(is.na(y$population))
#    rl <- data.frame(lengths = rl$lengths, missing = rl$values)
#    rl <- rbind(data.frame(lengths = 0, missing = FALSE), rl)
#
#    long_index <- (1:nrow(rl))[rl$missing == FALSE & (rl$lengths == max(rl$lengths[rl$missing == FALSE]))]
#    long_index <- long_index[length(long_index)] # in case there are multiple, take newest
#
#    start_i <- sum(rl$lengths[1:(long_index - 1)]) + 1
#    end_i <- sum(rl$lengths[1:(long_index)])
#    y[start_i:end_i, ]
#    })

gpdd <- plyr::ddply(gpdd, "main_id", function(y) {
  rl <- rle(is.na(y$population_untransformed))
  rl <- data.frame(lengths = rl$lengths, missing = rl$values)
  rl <- rbind(data.frame(lengths = 0, missing = FALSE), rl)

  long_index <- (1:nrow(rl))[rl$missing == FALSE &
    (rl$lengths == max(rl$lengths[rl$missing == FALSE]))]
  long_index <- long_index[length(long_index)] # in case there are multiple, take newest

  start_i <- sum(rl$lengths[1:(long_index - 1)]) + 1
  end_i <- sum(rl$lengths[1:(long_index)])
  message(paste0("Keeping rows ", start_i, ":", end_i))
  return(y[start_i:end_i, ])
  })

# And now remove those that have fewer than N points:
gpdd <- plyr::ddply(gpdd, "main_id", function(y) {
  if(nrow(y) >= 20) y
})

# a bit of manual tuning to note in ms:
# Remove this heron population, it is duplicated:
# main_id = 20531 and 10139
# jj <- subset(gpdd, main_id %in% c(20531, 10139))[,c("main_id", "series_step", "population")]
# ggplot(jj, aes(series_step, population, colour = as.factor(main_id))) + geom_point()
gpdd <- subset(gpdd, main_id != 20531)
gpdd <- mutate(gpdd, label = paste(main_id, common_name))

# looks like this value should be interpolated not zero subbed:
# main ID: 20659
gpdd$population_untransformed[gpdd$data_id == 1022378] <-
  exp(mean(c(log(7643), log(7117))))
gpdd$zero_sub[gpdd$data_id == 1022378] <- FALSE
gpdd$interpolated[gpdd$data_id == 1022378] <- TRUE

gpdd %>% group_by(taxonomic_class, main_id) %>% summarise(n = n()) %>%
  summarise(n = n())

# now focus the taxonomy slightly:
gpdd <- subset(gpdd, taxonomic_class != "Angiospermopsida (Dicotyledoneae)")
gpdd <- subset(gpdd, taxonomic_class != "Angiospermopsida (Monocotyledonae)")
gpdd <- subset(gpdd, taxonomic_class != "Bacillariophyceae")
gpdd <- subset(gpdd, taxonomic_class != "Unknown")

gpdd <- subset(gpdd, !grepl("^[A-Za-z]+ $", taxon_name))
gpdd <- subset(gpdd, !grepl("^[A-Za-z]+$", taxon_name))
gpdd <- subset(gpdd, !grepl("Unknown", taxon_name))

gpdd %>% group_by(main_id, common_name) %>%
  summarise(rel = reliability[1]) %>%
  arrange(desc(rel)) %>%
  group_by(rel) %>%
  summarise(n = n()) %>% as.data.frame

# and remove those that we zero-subbed more than N times:
# this doesn't remove many, but it removes a few silly ones
# with a ton of zero subs
gpdd <- plyr::ddply(gpdd, "main_id", function(x) {
  if(sum(x$zero_sub) <= 5)
    x
})

# and remove those coded by GPDD as reliability == 1
# mostly hidden fur records
gpdd <- plyr::ddply(gpdd, "main_id", function(x) {
  rel <- as.numeric(na.omit(x$reliability)) # some were turned to NAs
  if(rel[1] != 1)
    x
})

stat_table <-
  gpdd %>% group_by(taxonomic_class, main_id) %>%
  summarise(x = n(), x_zero_sub = sum(zero_sub),
    x_impute = sum(interpolated), x_sp = taxon_name[1],
    x_order = taxonomic_order[1]) %>%
  summarise(
    n = n(),
    n_orders = length(unique(x_order)),
    n_species = length(unique(x_sp)),
    #min_N = min(x),
    median_N = median(x),
    max_N = max(x),
    n_interpolated = sum(x_impute),
    n_zero_sub = sum(x_zero_sub)
    ) %>% arrange(desc(n)) %>% as.data.frame
write.csv(stat_table, file = "stat_table.csv", row.names = FALSE)

stat_table$max_N <- NULL

names(stat_table) <- c("Taxonomic class", "Populations", "Orders", "Species",
  "Median length", "Interpolated pts", "Zeros pts")
library(xtable)
print.xtable(xtable(stat_table,
    caption = "Summary statistics for the filtered Global Population Dynamics Database time series arranged by taxonomic class. Columns are: number of populations, number of taxonomic orders, numbers of species, median time series length, total number of interpolated time steps, and total number of substituted zeros."),
  include.rownames = FALSE, file = "stat-table.tex",
  booktabs = TRUE,  caption.placement = "top")

saveRDS(gpdd, file = "gpdd-clean.rds")

# massive ts plot:
library(ggplot2)
p <- ggplot(gpdd, aes(series_step, log10(population_untransformed),
  colour = taxonomic_class)) +
  geom_point() + geom_line() +
  facet_wrap(~label, scales = "free") +
  geom_point(data = subset(gpdd, interpolated == TRUE),
    aes(series_step, log10(population_untransformed)), colour = "black", pch = 21) +
  geom_point(data = subset(gpdd, zero_sub == TRUE),
    aes(series_step, log10(population_untransformed)), colour = "black", pch = 20) +
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    strip.background = element_blank()) +
  xlab("Time series step") + ylab("log10(Abundance)")
ggsave("all-clean-ts-3.pdf", width = 50, height = 50, limitsize = FALSE)

p <- ggplot(subset(gpdd, assumed_log10 == TRUE),
  aes(series_step, log10(population_untransformed))) +
  geom_point() + geom_line() +
  facet_wrap(~label, scales = "free")
ggsave("log10-assumed.pdf", width = 20, height = 20)
