# plot the simulation output from 3.3-test-gomp-models.R

library("dplyr")

check_nu <- readRDS("check_nu.rds")

load("nu_effective_seeds.rda")

effective_df <- data.frame(N = c(50, 50, 50, 100, 100, 100), nu_true = c(3, 5, 1e9, 3, 5, 1e9), nu_effective = c(median(nu_3_seeds_N50$nu_ests), median(nu_5_seeds_N50$nu_ests), 1e9, median(nu_3_seeds_N100$nu_ests), median(nu_5_seeds_N100$nu_ests), 1e9))

check_nu <- mutate(check_nu, inverse_nu = 1/nu, inverse_nu_l = 1/nu_l,
  inverse_nu_u = 1/nu_u, inverse_nu_true = 1/nu_true)

id.vars <- c("id", "nu_true", "N", "sigma_obs_true", "sigma_obs_assumed")
measure.vars <- c("inverse_nu", "b", "sigma_proc", "lambda")
check_nu_est <- reshape2::melt(check_nu, id.vars = id.vars,
  measure.vars = measure.vars)
check_nu_l <- reshape2::melt(check_nu, id.vars = id.vars,
  measure.vars = paste0(measure.vars, "_l"))
check_nu_u <- reshape2::melt(check_nu, id.vars = id.vars,
  measure.vars = paste0(measure.vars, "_u"))

check_nu_l <- plyr::rename(check_nu_l, c("value" = "l"))
check_nu_u <- plyr::rename(check_nu_u, c("value" = "u"))
check_nu_est <- plyr::rename(check_nu_est, c("value" = "est"))

check_nu_long <- data.frame(check_nu_l, est = check_nu_est[,"est"],
  u = check_nu_u[,"u"])
rm(check_nu_l, check_nu_u, check_nu_est)
check_nu_long$variable <- as.character(check_nu_long$variable)
check_nu_long$variable <- sub("([a-z]*)_l", "\\1",
  check_nu_long$variable)

## now bring in the true values for plotting horizontal lines:
measure.vars <- paste0(measure.vars, "_true")
check_nu_true <- reshape2::melt(check_nu, id.vars = id.vars,
  measure.vars = measure.vars)
check_nu_true <- plyr::rename(check_nu_true, c("variable" = "true_variable",
    "value" = "true_value"))
check_nu_long <- data.frame(check_nu_long,
  true_value = check_nu_true[, "true_value"])
check_nu_long$true_variable <- NULL # not needed

check_nu_long <- left_join(check_nu_long, effective_df)
check_nu_long$true_value[check_nu_long$variable == "inverse_nu"] <-
  1/check_nu_long$nu_effective[check_nu_long$variable == "inverse_nu"]

limits <- data.frame(variable = c("inverse_nu", "sigma_proc", "lambda", "b"), ylim_l = c(0, 0, 0, -0.25), ylim_u = c(0.5, 1.35, 12, 1), stringsAsFactors = FALSE)
check_nu_long <- plyr::join(check_nu_long, limits)

# and sort levels for plotting:
check_nu_long$variable <- factor(check_nu_long$variable,
  levels = c("inverse_nu", "sigma_proc", "lambda", "b"))

check_nu_long <- mutate(check_nu_long, scenario = paste(N, sigma_obs_true, sigma_obs_assumed, sep = "-"))
## long version of plots:
scenario_df <- data.frame(scenario = c("50-0.001-0.001", "100-0.001-0.001",
  "50-0.2-0.001", "50-0.2-0.2"), scen_short = c(2, 1, 3, 4),
  Scenario = c("2. N50, no obs. error", "1. N100, no obs. error",
    "3. N50, obs. error ignored", "4. N50, obs. error modelled"),
  stringsAsFactor = FALSE)
check_nu_long <- plyr::join(check_nu_long, scenario_df)

# and bring in p_0.1 for plotting:
measure.vars <- c("p_0.1")
check_nu_p10_long <- reshape2::melt(check_nu, id.vars = id.vars,
  measure.vars = measure.vars)
check_nu_p10_long <- plyr::rename(check_nu_p10_long, c("value" = "p10"))
check_nu_p10_long$variable <- NULL

check_nu_long <- plyr::join(check_nu_long, check_nu_p10_long)

library("ggplot2")
theme_set(theme_bw())

make_boxpanel <- function(dat, title) {
  ggplot(dat, aes(inverse_nu_true, p_0.1, group = inverse_nu_true)) + geom_point(position = position_jitter(width = 0.01), alpha = 0.6, colour = "black") + ylab("Probability nu < 10") + ylim(0, 1.01) + ggtitle(title) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.333, 0.5), labels = c("Infinity", "10", "5", "3", "2")) + xlab("nu true")
}

p1 <- make_boxpanel(filter(check_nu, sigma_obs_true == 0.001, N == 100), title = "N = 100, sigma_obs = 0.001, sigma_obs ignored")
p2 <- make_boxpanel(filter(check_nu, sigma_obs_true == 0.001, N == 50), title = "N = 50, sigma_obs = 0.001, sigma_obs ignored")
p3 <- make_boxpanel(filter(check_nu, sigma_obs_true == 0.2, sigma_obs_assumed == 0.001), title = "N = 100, sigma_obs = 0.2, sigma_obs ignored")
p4 <- make_boxpanel(filter(check_nu, sigma_obs_true == 0.2, sigma_obs_assumed == 0.2), title = "N = 100, sigma_obs = 0.2, sigma_obs assumed")
pdf("check-sim-box.pdf", width = 11, height = 6)
gridExtra::grid.arrange(p1, p2, p3, p4)
dev.off()

p <- ggplot(check_nu_long, aes(scen_short, est, group = id, colour = Scenario)) + geom_pointrange(aes(ymin = l, ymax = u), lwd = 0.1, position = position_dodge(width = 0.7)) + geom_hline(aes(yintercept = true_value), col = "black", lty = 2, lwd = 1) + facet_grid(variable~nu_true, scales = "free_y") + xlab("Scenario") + ylab("Parameter estimate")
ggsave("sim-gompertz.pdf", width = 11, height = 6)

p <- ggplot(check_nu_long, aes(scen_short, est, colour = Scenario)) + geom_boxplot() +  facet_grid(variable~nu_true, scales = "free_y") + xlab("Scenario") + ylab("Parameter estimate") + geom_hline(aes(yintercept = true_value), col = "black", lty = 2, lwd = 1)
ggsave("sim-gompertz-boxplots.pdf", width = 10, height = 6)

p <- filter(check_nu_long, variable == "inverse_nu") %>%
  ggplot(aes(scen_short, est, colour = Scenario, group = Scenario)) +
  geom_boxplot(alpha = 0.5, outlier.size = 0, colour = "grey50") +
  #geom_point(colour = "grey50") +
  geom_point(position = position_jitter(width = 0.15), alpha = 0.7) +
  ggtitle(expression(nu~true)) +
  facet_grid(~nu_true) +
  xlab("Scenario") + ylab(expression(Median~widehat(nu))) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.333, 0.5),
    labels = c("Infinity", "10", "5", "3", "2"))  +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) +
  geom_hline(aes(yintercept = true_value), col = "black", lty = 2, lwd = 1)
ggsave("sim-gompertz-median-dist.pdf", width = 9, height = 3)

p <- filter(check_nu_long, variable == "inverse_nu") %>%
  ggplot(aes(1/nu_true, p10, group = nu_true)) +
  geom_point(position = position_jitter(width = 0.01),
    alpha = 0.6, colour = "black") +
  facet_grid(~Scenario) +
  ylab("Probability nu < 10") +
  xlab(expression(nu~true)) +
  ylim(0, 1.01) +
  ggtitle("") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()) +
  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.333, 0.5),
    labels = c("Infinity", "10", "5", "3", "2"))

p <- filter(check_nu_long, variable == "inverse_nu") %>%
  ggplot(aes(scen_short, p10, group = Scenario, colour = Scenario)) +
  geom_boxplot(alpha = 0.5, outlier.size = 0, colour = "grey50") +
  geom_point(position = position_jitter(width = 0.15),
    alpha = 0.7) +
  facet_grid(~nu_true) +
  ylab(expression(Pr(nu < 10))) +
  xlab("Scenario") +
  ylim(0, 1.01) +
  ggtitle(expression(nu~true)) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank())
  #scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.333, 0.5),
    #labels = c("Infinity", "10", "5", "3", "2"))

ggsave("sim-gompertz-p10.pdf", width = 9, height = 3)
