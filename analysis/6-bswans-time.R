# check - is there a trend through time or an alignment of black swans?

library("dplyr")
source("5-shape-data.R")

heavy <- filter(gomp_hat_base, nu_50 < 20) %>%
  select(main_id, common_name, p10)

heavy <- plyr::join(heavy, as.data.frame(gpdd)) %>%
  arrange(main_id, p10) %>%
  mutate(id_name = paste(round(p10, 2), common_name, main_id))

get_gomp_res <- function(id_show) {
  p <- subset(gomp_hat_base, main_id == id_show)
  pop <- subset(gpdd, main_id == id_show)$population_untransformed
  res <- rep(NA, length(pop))
  for(i in 2:length(pop)) {
    res[i] <- log(pop)[i] - (p$lambda_50 + p$b_50 * log(pop[i-1]))
      #p$phi_50 * (p$lambda_50 + p$b_50 * log(pop[i-2]))
  }
  l <- qnorm(0.0001, 0, sd = p$sigma_proc_50)
  u <- qnorm(0.9999, 0, sd = p$sigma_proc_50)
  return(list(res = res, l = l, u = u))
}

heavy_res <- plyr::ddply(heavy, "main_id", function(x) {
  qq <- get_gomp_res(x$main_id[1])
  with(qq, data.frame(res, l, u, sample_year = x$sample_year))
})
#saveRDS(heavy_res, file = "heavy_residuals.rds")

heavy_res <- heavy_res %>% group_by(main_id) %>%
  mutate(
    bs = ifelse(res < l | res > u, TRUE, FALSE),
    bs_type = ifelse(bs & res > u, "up", ifelse(bs & res < l, "down", "none"))) %>%
  filter(sample_year > 0) %>%
  as.data.frame

xx <- heavy_res %>%
  mutate(decade = paste0(substr(sample_year, 1, 3), "0")) %>%
  group_by(decade) %>%
  summarise(bs_count = sum(bs, na.rm = TRUE), samples = length(bs)) %>%
  mutate(bs_perc = bs_count / samples) %>%
  as.data.frame
plot(xx$decade, xx$bs_perc)

xx <- plyr::ldply(
  seq(min(heavy_res$sample_year),
    max(heavy_res$sample_year) - 10),
  function(i) {
out <-  filter(heavy_res, sample_year > i & sample_year < i + 10) %>%
  summarise(bs_count = sum(bs, na.rm = TRUE), samples = length(bs)) %>%
  mutate(bs_perc = bs_count / samples, year = i)
    })
plot(xx$year, xx$bs_perc)


