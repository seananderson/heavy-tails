# plot order-level posterior predictions
# start by running 5.8...

library("rstan")
library("ggplot2")

if(!file.exists("betareg5.rds") {
  stan_beta5 <- stan_model("betareg5.stan")
  saveRDS(stan_beta5, "betareg5.rds")
} else {
  stan_beta5 <- readRDS("betareg5.rds")
}
d <- readRDS("beta-modelling-dat.rds")

m.stan.beta5 <- sampling(stan_beta5,
  data = list(
    N = nrow(d),
    n_order = max(d$order_id),
    n_class = max(d$class_id),
    n_sp = max(d$sp_id),
    x1 = as.numeric(d$log_dataset_length_scaled),
    order_id = d$order_id,
    class_id = d$class_id,
    sp_id = as.numeric(d$sp_id),
    y = d$p10),
  pars = c("b1", "mu_a",
    "sigma_a_class", "sigma_a_order", "sigma_a_sp", "phi",
    "a_class", "a_order"),
  iter = 5000, chains = 4, thin = 2)
saveRDS(m.stan.beta5, file = "beta-stan-samples-n-only.rds")
sink("beta-stan-stamples-n-only.txt")
print(m.stan.beta5)
sink()

m <- readRDS("beta-stan-samples-n-only.rds")

lu <- d[,c("order_id", "class_id", "taxonomic_class", "taxonomic_order")]
lu <- lu[!duplicated(lu), ]

a_class <- apply(extract(m, pars = "a_class")[[1]], 2, median)
a_class_df <- data.frame(a_class = a_class, class_id = 1:length(a_class))
a_class_df <- plyr::join(a_class_df, lu)

a_order <- apply(extract(m, pars = "a_order")[[1]], 2, median)
a_order_df <- data.frame(a_order = a_order, order_id = 1:length(a_order))
a_order_df <- plyr::join(a_order_df, lu)

a_df <- plyr::join(a_order_df, a_class_df)
a_df$mu_a <- median(extract(m, pars = "mu_a")[[1]])
a_df$a <- a_df$mu_a + a_df$a_order + a_df$a_class
a_df$invlogit_a <- plogis(a_df$a)

# get sample sizes for each order:
ns <- d %>% group_by(taxonomic_order) %>% dplyr::summarise(n_pops = n())
a_df <- plyr::join(a_df, ns)
a_df <- subset(a_df, taxonomic_order != "Unknown")
a_df <- subset(a_df, n_pops >= 5)
a_df <- droplevels(a_df)

a_df$sorted_order <- factor(a_df$taxonomic_order,
  levels(a_df$taxonomic_order)[rev(order(a_df$invlogit_a))])

###
mu_a <- extract(m, pars = "mu_a")[[1]]
a_class <- extract(m, pars = "a_class")[[1]]
a_order <- extract(m, pars = "a_order")[[1]]

# get order:
o <- plyr::ldply(unique(a_df$order_id), function(x) {
  i_order <- x
  i_class <- unique(lu[lu$order_id == x, "class_id"])
  a <- plogis(mu_a + a_order[, i_order] + a_class[,i_class])
  median(a)
})
o <- data.frame(order_id = unique(a_df$order_id), median_dens = o$V1)
o <- o[order(o$median_dens), ]

a_df <- plyr::join(a_df, o)
a_df <- a_df[order(a_df$median_dens), ]

pal <- RColorBrewer::brewer.pal(4, "Set3")

cols_df <-
  data.frame(col = c(RColorBrewer::brewer.pal(4, "Set3")[c(3, 4, 1, 2)],
    rep("#FFFFFF", 3)),
  taxonomic_class = c("Aves", "Mammalia", "Insecta", "Osteichthyes",
    "Chondrichtyhes", "Crustacea", "Gastropoda"),
  stringsAsFactors = FALSE)
lu <- plyr::join(lu, cols_df)

x <- rexp(1e7, 0.01)
x <- x[x > 2]
prior_p10 <- length(x[x < 10])/length(x)
