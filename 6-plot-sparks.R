# make little time series sparklines (copyright Tufte) of the heavy-tailed
# populations for a table


#source("5-shape-data.R")

heavy <- filter(gomp_hat_base, nu_50 < 20) %>%
  select(main_id, common_name, p10)

heavy <-
  plyr::join(heavy, gpdd) %>%
  arrange(main_id, p10) %>%
  mutate(id_name = paste(round(p10, 2), common_name, main_id))

get_gomp_res <- function(id_show) {
  p <- subset(gomp_hat_base, main_id == id_show)
  pop <- subset(gpdd, main_id == id_show)$population_untransformed
  res <- rep(NA, length(pop))
  for(i in 3:length(pop)) {
    res[i] <- log(pop)[i] - (p$lambda_50 + p$b_50 * log(pop[i-1]))
      #p$phi_50 * (p$lambda_50 + p$b_50 * log(pop[i-2]))
  }
  l <- qnorm(0.0001, 0, sd = p$sigma_proc_50)
  u <- qnorm(0.9999, 0, sd = p$sigma_proc_50)
  return(list(res = res, l = l, u = u))
}

make_spark <- function(x) {
  res <- get_gomp_res(x$main_id[1])
  bsw_l <- as.numeric(na.omit(seq_along(res$res)[res$res<res$l]))
  bsw_u <- as.numeric(na.omit(seq_along(res$res)[res$res>res$u]))
  pdf(paste0("sparks/", x$main_id[1], ".pdf"), width = 2.5, height = 1.2)
  par(mar = c(0.5,0.2,0.5,0.2), oma = c(0,0,1.0,0), cex = 0.8)
  par(xpd = NA)
  with(x, plot(seq_along(population_untransformed), population, type = "l", axes = FALSE,
      xlab = "", ylab = "", yaxs = "i", xaxs = "i"))
  if(length(bsw_l) > 0)
    points(bsw_l, x$population_untransformed[bsw_l], col = "red", pch = 20, cex = 3)
  if(length(bsw_u) > 0)
    points(bsw_u, x$population_untransformed[bsw_u], col = "blue", pch = 20, cex = 3)
  dev.off()
}

#make_spark(subset(heavy, main_id == 10139))

plyr::d_ply(heavy, "main_id", function(y) {
  make_spark(y)
})

#get_gomp_res(6528)
#get_gomp_res(20580)

plot(res)
abline(h = 0)
l <- qnorm(0.01, 0, sd = p$sigma_proc_50)
u <- qnorm(0.99, 0, sd = p$sigma_proc_50)
bsw <- seq_along(res)[res<l]

with(x, plot(seq_along(population), population, type = "l", axes = FALSE,
  xlab = "", ylab = "", yaxs = "i", xaxs = "i"))
points(bsw, pop[bsw], col = "red", pch = 20, cex = 2)

# output .csv to get started
h <- heavy %>% plyr::ddply("main_id", function(x) x[1,]) %>% arrange(desc(p10))
h <- h[,c("main_id", "p10", "common_name", "taxon_name", "exact_name", "ref", "notes", "data_notes")]


write.csv(h, file = "heavy.csv")
