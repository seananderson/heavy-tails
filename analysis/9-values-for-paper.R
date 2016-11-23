# this file produces some values for inclusion of the text of the associated paper

library("dplyr")
source("5-shape-data.R")
p_inc <- readRDS("prob_inc_heavy_with_n.rds")

write_tex <- function(x, macro, ...) {
  out <- paste0("\\newcommand{\\", macro, "}{", x, "}")
  cat(out, file = zz)
  cat("\n", file = zz)
}
zz <- file("values.tex", "w") # open .tex file to write to throughout

# # what's the median and mean nu in the exponential prior?
set.seed(1)
x <- rexp(3e6, 0.01)
x <- x[x > 2]
write_tex(round(mean(x), 0) , "basePriorMean")
write_tex(round(median(x), 0) , "basePriorMedian")
write_tex(round(100*length(x[x < 10])/length(x), 1) , "basePriorProbHeavy")

# how much was imputed?
gpdd <- readRDS("gpdd-clean.rds")
perc_imputed <- gpdd %>% group_by(main_id) %>%
  summarise(has_imputed = ifelse(max(interpolated) == 1, TRUE, FALSE)) %>%
  summarise(percent_imputed = sum(has_imputed)/length(has_imputed)) %>%
  round(2) * 100
perc_imputed_pops <- as.character(perc_imputed$percent_imputed)

perc_imputed_points <- round(sum(gpdd$interpolated)/nrow(gpdd) * 100, 1)

fig2_pop_n <- filter(gpdd, taxonomic_class %in%
    c("Aves", "Mammalia", "Insecta", "Osteichthyes")) %>%
  summarise(n = length(unique(main_id))) %>% summarise(N = sum(n))
fig2_pop_n <- fig2_pop_n$N

NPops <- length(unique(gpdd$main_id))
NOrders <- length(unique(gpdd$taxonomic_order))
NClasses <- length(unique(gpdd$taxonomic_class))

cols_include <- c("taxonomic_class", "main_id", "p10", "nu_50", "type",
  "max_rhat", "min_neff")
gomp_hat_base$ type <- "base"
gomp_hat_ar1$type <- "ar1"
gomp_hat_logistic$type <- "logistic"
gomp_hat_obs_0.2$type <- "obs"

gtemp <- rbind(
  gomp_hat_base[,cols_include],
  gomp_hat_ar1[,cols_include],
  gomp_hat_logistic[,cols_include],
  gomp_hat_obs_0.2[,cols_include])

pheavy_class <- gtemp %>% filter(max_rhat < 1.05) %>%
  group_by(taxonomic_class, type) %>%
  summarise(n = n(), h =  length(which(p10 > 0.5)), p = round(100 * h / n)) %>%
  group_by(taxonomic_class) %>%
  summarise(min_p = min(p), max_p = max(p)) %>%
  filter(min_p > 0)

pheavy_type <- gtemp %>% filter(max_rhat < 1.05) %>%
  group_by(type) %>%
  summarise(n = n(), h =  length(which(p10 > 0.5)), p = round(100 * h / n))

pheavy_overall <- gtemp %>%
  group_by(type) %>%
  summarise(n = n(), h = length(which(p10 > 0.5)), p = round(100 * h / n))

## To obtain:
## Methods:
## - median time steps
## - range of time steps (methods 2nd para.)
step_stats <- gpdd %>% group_by(main_id) %>% summarise(time_steps = n()) %>%
  summarise(
    p25 = round(quantile(time_steps, probs = 0.25)),
    p50 = quantile(time_steps, probs = 0.50),
    p75 = round(quantile(time_steps, probs = 0.75)),
    max = max(time_steps),
    min = min(time_steps),
    mean = round(mean(time_steps), 1))
write_tex(step_stats$p50, "medianTimeSteps")
write_tex(step_stats$mean, "meanTimeSteps")
write_tex(step_stats$min, "minTimeSteps")
write_tex(step_stats$max, "maxTimeSteps")

## Results:
## - Pr(v < 10) > 0.5:
## - counts for birds, mammals, insects, fishes
fract_heavy <- gomp_hat_base %>%
  group_by(taxonomic_class) %>%
  summarise(n_all = n(), n_heavy =length(which(p10 > 0.5))) %>%
  filter(n_all > 15) %>%
  mutate(p_heavy = round(100 * n_heavy / n_all))
write_tex(subset(fract_heavy, taxonomic_class == "Aves")$n_all, "birdN")
write_tex(subset(fract_heavy, taxonomic_class == "Insecta")$n_all, "insectsN")
write_tex(subset(fract_heavy, taxonomic_class == "Mammalia")$n_all, "mammalsN")
write_tex(subset(fract_heavy, taxonomic_class == "Osteichthyes")$n_all, "fishN")
write_tex(subset(fract_heavy, taxonomic_class == "Aves")$n_heavy, "birdNH")
write_tex(subset(fract_heavy, taxonomic_class == "Insecta")$n_heavy, "insectsNH")
write_tex(subset(fract_heavy, taxonomic_class == "Mammalia")$n_heavy, "mammalsNH")
write_tex(subset(fract_heavy, taxonomic_class == "Osteichthyes")$n_heavy, "fishNH")
# percentages:
write_tex(subset(fract_heavy, taxonomic_class == "Aves")$p_heavy, "birdPH")
write_tex(subset(fract_heavy, taxonomic_class == "Insecta")$p_heavy, "insectsPH")
write_tex(subset(fract_heavy, taxonomic_class == "Mammalia")$p_heavy, "mammalsPH")
write_tex(subset(fract_heavy, taxonomic_class == "Osteichthyes")$p_heavy, "fishPH")

## - how many orders had at least one black swan population?
NOrdersHeavy <- gomp_hat_base %>% group_by(taxonomic_order) %>%
  summarise(n_heavy = length(which(p10 > 0.5))) %>%
  filter(n_heavy >= 1) %>%
  nrow
write_tex(NOrdersHeavy, "NOrdersHeavy")
write_tex(round(100*NOrdersHeavy/NOrders), "POrdersHeavy")

## - how many populations with nu < 10 to nu > 10 with obs. error?
base_heavy_main_ids_0.50 <- subset(gomp_hat_base, p10 > 0.50)[, "main_id"]
base50Obs50Switch <- filter(gomp_hat_obs_0.2,
  main_id %in% base_heavy_main_ids_0.50 & p10 <= 0.5) %>% nrow
write_tex(base50Obs50Switch, "baseFiftyObsFiftySwitch")
base_heavy_main_ids_0.75 <- subset(gomp_hat_base, p10 > 0.75)[, "main_id"]
base75Obs50Switch <- filter(gomp_hat_obs_0.2,
  main_id %in% base_heavy_main_ids_0.75 & p10 <= 0.5) %>% nrow

write_tex(base75Obs50Switch, "baseSeventyFiveObsFiftySwitch")
write_tex(length(base_heavy_main_ids_0.50), "totalHeavyFifty")
write_tex(length(base_heavy_main_ids_0.75), "totalHeavySeventyFive")
write_tex(round(100*base50Obs50Switch/length(base_heavy_main_ids_0.50)),
  "baseFiftyObsFiftySwitchPerc")
write_tex(round(100*base75Obs50Switch/length(base_heavy_main_ids_0.50)),
  "baseSeventyFiveObsFiftySwitchPerc")

# now with respect to median estimates:
base_heavy_main_ids_lt10 <- subset(gomp_hat_base, nu_50 <= 10)[, "main_id"]
baseNuTenObsTenSwitch <- filter(gomp_hat_obs_0.2,
  main_id %in% base_heavy_main_ids_lt10 & nu_50 > 10) %>% nrow
write_tex(baseNuTenObsTenSwitch, "baseNuTenObsTenSwitch")
write_tex(length(base_heavy_main_ids_lt10), "baseNuTen")

base_heavy_main_ids_lt5 <- subset(gomp_hat_base, nu_50 <= 5)[, "main_id"]
baseNuFiveObsTenSwitch <- filter(gomp_hat_obs_0.2,
  main_id %in% base_heavy_main_ids_lt5 & nu_50 > 10) %>% nrow
write_tex(baseNuFiveObsTenSwitch, "baseNuFiveObsTenSwitch")

## - how much more probable were heavy tails with 60 time steps. vs. 30 time
##   steps? (percentage and absolute)
pHeavyN30 <- p_inc[p_inc$N == 30, "p"]
pHeavyN60 <- p_inc[p_inc$N == 60, "p"]
pIncHeavyN30N60 <- pHeavyN60 / pHeavyN30
write_tex(sprintf("%.2f", round(pHeavyN30, 2)), "pHeavyNThirty")
write_tex(sprintf("%.2f", round(pHeavyN60, 2)), "pHeavyNSixty")
write_tex(sprintf("%.1f", round(pIncHeavyN30N60, 1)), "pIncHeavyNThirtyNSixty")

## Discussion:
## - Ricker / gompertz range of pops with black swans as percent
##
## Supplement:
## - "the model still categorized XX\% of cases as heavy tailed when v = 5"
##   (with observation error)
check_nu <- readRDS("check_nu.rds")

n_correct <- nrow(filter(check_nu, sigma_obs_true == 0.2,
  sigma_obs_assumed == 0.001, p_0.1 > 0.5, nu_true == 5))
n_tot <- max(check_nu$id)
obsErrorNuFivePerc <- round(100 * (n_correct / n_tot))
write_tex(obsErrorNuFivePerc, "obsErrorNuFivePerc")

## - how many autocorrelation residual Gompertz models didn't converge ("MCMC
##   chains for a small number..."
modelsNoConvergeAROne <- nrow(filter(gomp_hat_ar1,
  max_rhat > 1.05 | min_neff <  200))
write_tex(modelsNoConvergeAROne, "modelsNoConvergeAROne")

ar1_no_converge <- filter(gomp_hat_ar1,
  max_rhat > 1.05 | min_neff <  200) %>% select(main_id)

modelsNoConvergeAROneHeavyBase <-
  filter(gomp_hat_base, main_id %in% ar1_no_converge$main_id) %>%
  filter(p10 > 0.5) %>% nrow
write_tex(modelsNoConvergeAROneHeavyBase, "modelsNoConvergeAROneHeavyBase")

write_tex(perc_imputed_pops, "percImputedPops")
write_tex(perc_imputed_points, "percImputedPoints")
write_tex(fig2_pop_n, "nuCoefPopN")
#write_tex(total_assumed_log10, "totalAssumedLog")
for(cl in unique(pheavy_class$taxonomic_class)) {
  write_tex(paste0(
    pheavy_class[pheavy_class$taxonomic_class == cl, "min_p"],
    "--",
    pheavy_class[pheavy_class$taxonomic_class == cl, "max_p"]),
    paste0(cl, "RangePerc"))
}
if(!"Osteichthyes" %in% unique(pheavy_class$taxonomic_class)) {
  write_tex(0, "OsteichthyesRangePerc")
}

write_tex(min(pheavy_overall$p), "overallMinPerc")
write_tex(max(pheavy_overall$p), "overallMaxPerc")
write_tex(pheavy_overall$p[pheavy_overall$type == "base"], "overallBasePerc")
write_tex(NPops, "NPops")
write_tex(NOrders, "NOrders")
write_tex(NClasses, "NClasses")

stat_table <- read.csv("stat_table.csv", stringsAsFactors = FALSE)

interpPointsPerc <- round((sum(stat_table$n_zero_sub) + sum(stat_table$n_interpolated)) / sum(stat_table$n_points) * 100, 0)

write_tex(interpPointsPerc, "interpPointsPerc")

# count ups and downs:
heavy_res <- readRDS("heavy_residuals.rds")
heavy_res <- heavy_res %>% mutate(bs_l = res < l, bs_u = res > u)
n_bs_up <- sum(heavy_res$bs_u, na.rm = TRUE)
n_bs_down <- sum(heavy_res$bs_l, na.rm = TRUE)
write_tex(n_bs_up, "nBSUp")
write_tex(n_bs_down, "nBSDown")
write_tex(sprintf("%.1f", round(n_bs_down/n_bs_up, 1)), "ratioBSDownToUp")
write_tex(round(n_bs_down / (sum(n_bs_down, n_bs_up))*100), "percBSDown")

qq <- readRDS("skew-understimates.rds")
q <- quantile(1/qq$r, probs = c(0.25, 0.5, 0.75)) %>% round(1)
write_tex(paste0(q[[1]], "--", q[[3]]), "crashUnderRange")
write_tex(q[[2]], "crashUnderMedian")

sk <- readRDS("skew_samples.rds")
sk <- sk %>% group_by(main_id) %>%
  mutate(
    median_nu = median(nu),
    h1 = ifelse(median_nu <= 10, "1 heavy",
      ifelse(median_nu <= 70, "2 moderate", "3 normal"))) %>%
  ungroup()

normal <- filter(sk, h1 == "3 normal")$log_skew %>% sort
moderate <- filter(sk, h1 == "2 moderate")$log_skew %>% sort
heavy <- filter(sk, h1 == "1 heavy")$log_skew %>% sort
write_tex(round(mean(heavy < 0), 2) * 100, "probDensSkewedForHeavyPops")

overlap <- sk %>% filter(h1 == "3 normal") %>%
  group_by(main_id) %>%
    summarise(
      l = exp(quantile(log_skew, prob = 0.025)),
      u = exp(quantile(log_skew, prob = 0.975))) %>%
  mutate(overlap1 = l <= 1 & u >= 1)
write_tex(round(mean(overlap$overlap1), 2) * 100, "percNormPopsNotSkewed")

close(zz)
