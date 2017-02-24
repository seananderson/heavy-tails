# Make little time series sparklines (copyright Tufte) of the heavy-tailed
# populations for a table

source("5-shape-data.R")
dir.create("sparks")

heavy <- filter(gomp_hat_base, nu_50 < 20) %>%
  select(main_id, common_name, p10)

heavy <- plyr::join(heavy, as.data.frame(gpdd)) %>%
  arrange(main_id, p10) %>%
  mutate(id_name = paste(round(p10, 2), common_name, main_id))

get_gomp_res <- function(id_show) {
  p <- subset(gomp_hat_base, main_id == id_show)
  pop <- subset(gpdd, main_id == id_show)$population_untransformed
  yr <- subset(gpdd, main_id == id_show)$sample_year
  res <- rep(NA, length(pop))
  for(i in 2:length(pop)) {
    res[i] <- log(pop)[i] - (p$lambda_50 + p$b_50 * log(pop[i-1]))
  }
  l <- qnorm(0.0001, 0, sd = p$sigma_proc_50)
  u <- qnorm(0.9999, 0, sd = p$sigma_proc_50)
  return(list(res = res, l = l, u = u, yr = yr))
}

heavy_res <- plyr::ddply(heavy, "main_id", function(x) {
  qq <- get_gomp_res(x$main_id[1])
  with(qq, data.frame(res, l, u, yr))
})
saveRDS(heavy_res, file = "heavy_residuals.rds")

make_spark <- function(x) {
  res <- get_gomp_res(x$main_id[1])
  bsw_l <- as.numeric(na.omit(seq_along(res$res)[res$res<res$l]))
  bsw_u <- as.numeric(na.omit(seq_along(res$res)[res$res>res$u]))
  pdf(paste0("sparks/", x$main_id[1], ".pdf"), width = 2.5, height = 1.2)
  par(mar = c(0.5,0.2,0.5,0.2), oma = c(0,0,0.1,0), cex = 0.8)
  par(xpd = NA)
  with(x, plot(seq_along(population_untransformed), population_untransformed, type = "l", axes = FALSE,
      xlab = "", ylab = "", yaxs = "i", xaxs = "i", log = "y"))

# check for interpolation and blacks swan association:
  if(length(bsw_l) > 0)
    points(bsw_l, x$population_untransformed[bsw_l], col = "red", pch = 20, cex = 3)
  if(length(bsw_u) > 0)
    points(bsw_u, x$population_untransformed[bsw_u], col = "blue", pch = 20, cex = 3)
  dev.off()
}

plyr::d_ply(heavy, "main_id", function(y) {
  make_spark(y)
})

# output .csv to get started
h <- heavy %>% plyr::ddply("main_id", function(x) x[1,]) %>% arrange(desc(p10))
h <- h[,c("main_id", "p10", "common_name", "taxon_name", "exact_name", "ref", "notes", "data_notes")]

write.csv(h, file = "heavy.csv")

heavy_res$bs <- ifelse(heavy_res$res < heavy_res$l | heavy_res$res > heavy_res$u,
  TRUE, FALSE)
b <- filter(heavy_res, bs == TRUE, yr > 0)
b <- b %>% filter(res < l) %>% mutate(bs_frac = res / l)
ggplot(b, aes(yr, bs_frac)) + geom_point()

yrs <- gpdd %>% group_by(sample_year) %>% summarise(n = n()) %>%
  filter(sample_year > 0)
b_yrs <- b %>% group_by(yr) %>% summarise(nbs = n()) %>%
  rename(sample_year = yr) %>% inner_join(yrs) %>%
  mutate(perc_bs = nbs / n)
ggplot(b_yrs, aes(sample_year, perc_bs)) + geom_point()

