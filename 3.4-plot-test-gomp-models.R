# plot the simulation output from 3.3-test-gomp-models.R

check_nu <- readRDS("check_nu.rds")

check_nu <- mutate(check_nu, inverse_nu = 1/nu, inverse_nu_l = 1/nu_l,
  inverse_nu_u = 1/nu_u, inverse_nu_true = 1/nu_true)

id.vars <- c("id", "nu_true", "N", "sigma_obs_true", "sigma_obs_assumed")
measure.vars <- c("inverse_nu", "b", "sigma_proc", "lambda", "phi")
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

limits <- data.frame(variable = c("inverse_nu", "sigma_proc", "lambda", "b", "phi"), ylim_l = c(0, 0, 0, -0.25, -1), ylim_u = c(0.5, 1.35, 12, 1, 1), stringsAsFactors = FALSE)
check_nu_long <- plyr::join(check_nu_long, limits)

# and sort levels for plotting:
check_nu_long$variable <- factor(check_nu_long$variable,
  levels = c("inverse_nu", "sigma_proc", "lambda", "b", "phi"))

check_nu_long <- mutate(check_nu_long, scenario = paste(N, sigma_obs_true, sigma_obs_assumed, sep = "-"))
## long version of plots:
scenario_df <- data.frame(scenario = c("50-0.001-0.001", "100-0.001-0.001",
  "50-0.3-0.001", "50-0.3-0.3"), scen_short = c(2, 1, 3, 4),
  Scenario = c("2. N50, no obs. error", "1. N100, no obs. error",
    "3. N50, obs. error", "4. N50, obs. error modelled"), stringsAsFactor = FALSE)
check_nu_long <- plyr::join(check_nu_long, scenario_df)

library(ggplot2)

# make_panel <- function(dat, title) {`
#   ggplot(dat, aes(id, est, colour = variable)) + geom_pointrange(aes(ymin = l, ymax = u)) + geom_hline(aes(yintercept = true_value), col = "black", lty = 2, lwd = 1) + facet_grid(variable~nu_true, scales = "free_y") + xlab("Iteration") + theme_bw() + ylab("Parameter estimate")  + ggtitle(title) + theme(panel.grid.major = element_blank(),
#     panel.grid.minor = element_blank(), panel.background = element_blank()) + geom_hline(aes(yintercept = c(ylim_l, ylim_u)), col = "white", lwd = 0) + scale_colour_brewer(type = "qual", palette= "Set1", guide = FALSE)
# }
#
#
# p1 <- make_panel(filter(check_nu_long, N == 100), title = "N = 100, sigma_obs = 0.001, sigma_obs ignored")
# p2 <- make_panel(filter(check_nu_long, N == 50, sigma_obs_true == 0.001),
#   title = "N = 50, sigma_obs = 0.001, sigma_obs ignored")
# p3 <- make_panel(filter(check_nu_long, N == 50, sigma_obs_true == 0.3,
#     sigma_obs_assumed == 0.001), title = "N = 50, sigma_obs = 0.2; sigma_obs ignored")
# p4 <- make_panel(filter(check_nu_long, N == 50, sigma_obs_true == 0.3,
#     sigma_obs_assumed == 0.3), title = "N = 50, sigma_obs = 0.2; sigma_obs assumed")
#
# pdf("sim-check.pdf", width = 14, height = 11)
# gridExtra::grid.arrange(p1, p2, p3, p4)
# dev.off()
#
# p1 <- ggplot(check_nu_long, aes(, p_0.1, group = nu_true)) + geom_boxplot() + geom_point(position = position_jitter(width = 0.01), alpha = 0.3) + ylab("Probability nu < 10") + ylim(0, 1)

make_boxpanel <- function(dat, title) {
  ggplot(dat, aes(inverse_nu_true, p_0.1, group = inverse_nu_true)) + geom_point(position = position_jitter(width = 0.01), alpha = 0.6, colour = "black") + ylab("Probability nu < 10") + ylim(0, 1.01) + ggtitle(title) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.333, 0.5), labels = c("Infinity", "10", "5", "3", "2")) + xlab("nu true")
}

p1 <- make_boxpanel(filter(check_nu, sigma_obs_true == 0.001, N == 100), title = "N = 100, sigma_obs = 0.001, sigma_obs ignored")
p2 <- make_boxpanel(filter(check_nu, sigma_obs_true == 0.001, N == 50), title = "N = 50, sigma_obs = 0.001, sigma_obs ignored")
p3 <- make_boxpanel(filter(check_nu, sigma_obs_true == 0.3, sigma_obs_assumed == 0.001), title = "N = 100, sigma_obs = 0.3, sigma_obs ignored")
p4 <- make_boxpanel(filter(check_nu, sigma_obs_true == 0.3, sigma_obs_assumed == 0.3), title = "N = 100, sigma_obs = 0.3, sigma_obs assumed")
pdf("check-sim-box.pdf", width = 11, height = 6)
gridExtra::grid.arrange(p1, p2, p3, p4)
dev.off()


p <- ggplot(check_nu_long, aes(scen_short, est, group = id, colour = Scenario)) + geom_pointrange(aes(ymin = l, ymax = u), position = position_dodge(width = 0.5)) + geom_hline(aes(yintercept = true_value), col = "black", lty = 2, lwd = 1) + facet_grid(variable~nu_true, scales = "free_y") + xlab("Scenario") + ylab("Parameter estimate")
ggsave("sim-gompertz.pdf", width = 10, height = 6)

# + theme_bw() + ylab("Parameter estimate")

# + theme(panel.grid.major = element_blank()

#   panel.grid.minor = element_blank(), panel.background = element_blank()) + geom_hline(aes(yintercept = c(ylim_l, ylim_u)), col = "white", lwd = 0) + scale_colour_brewer(type = "qual", palette= "Set1", guide = FALSE)


## need to adapt these plots for the extra parameters now:
# library(ggplot2)
# p1 <- ggplot(check_nu, aes(id, 1/nu)) + geom_pointrange(aes(ymin = 1/nu_l, ymax = 1/nu_u)) + geom_hline(aes(yintercept = 1/nu_true), col = "red") + ylim(0, 0.5) + facet_grid(~nu_true) + xlab("Iteration")
#
# p2 <- ggplot(check_nu, aes(id, sigma_proc)) + geom_pointrange(aes(ymin = sigma_proc_l, ymax = sigma_proc_u)) + geom_hline(yintercept = 0.5, col = "red") + facet_grid(~nu_true)+ xlab("Iteration")
#
# p3 <- ggplot(check_nu, aes(id, b)) + geom_pointrange(aes(ymin = b_l, ymax = b_u)) + geom_hline(yintercept = 0.75, col = "red") + facet_grid(~nu_true)+ xlab("Iteration")
#
# p4 <- ggplot(check_nu, aes(id, lambda)) + geom_pointrange(aes(ymin = lambda_l, ymax = lambda_u)) + geom_hline(yintercept = 0.75, col = "red") + facet_grid(~nu_true)+ xlab("Iteration")
#
# p5 <- ggplot(check_nu, aes(id, phi)) + geom_pointrange(aes(ymin = phi_l, ymax = phi_u)) + geom_hline(yintercept = 0.1, col = "red") + facet_grid(~nu_true)+ xlab("Iteration")
#
# pdf("check-sim.pdf", width = 12, height = 7)
# gridExtra::grid.arrange(p1, p2, p3, p4, p5)
# dev.off()
#
# p1 <- ggplot(check_nu, aes(1/nu_true, p_0.1, group = nu_true)) + geom_boxplot() + geom_point(position = position_jitter(width = 0.01), alpha = 0.3) + ylab("Probability nu < 10") + ylim(0, 1)
#
# p2 <- ggplot(check_nu, aes(1/nu_true, p_0.2, group = nu_true)) + geom_boxplot() + geom_point(position = position_jitter(width = 0.01), alpha = 0.3) + ylab("Probability nu < 20") + ylim(0, 1)
#
# pdf("check-sim-p-nu.pdf", width = 6, height = 7)
# gridExtra::grid.arrange(p1, p2)
# dev.off()


# take a basic t- dist
# and try and re-capture it with our prior at different nu values
# use a ton of data (2000 points, 500 points, 100 points, 50 points)
# show how probability of noticing goes down
