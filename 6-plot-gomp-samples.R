# this file makes a variety of plots of the model posteriors

library(ggplot2)

gomp_hat_ar1 <- filter(gomp_hat_ar1, taxonomic_class != "Crustacea")
#gomp_hat_ar1 <- filter(gomp_hat_base, !taxonomic_class %in% c("Crustacea", "Chondrichtyhes", "Gastropoda"))

p <- ggplot(ar1_vs_base, aes(1/nu_50, 1/nu_50_base, col = phi_50)) + geom_point(cex = 3.5, alpha = 0.8) + geom_segment(aes(x = 1/nu_50, xend = 1/nu_50, y = 1/nu_25_base, yend = 1/nu_75_base), lwd = 0.5, alpha = 0.2) + geom_segment(aes(x = 1/nu_25, xend = 1/nu_75, y = 1/nu_50_base, yend = 1/nu_50_base), lwd = 0.5, alpha = 0.2)+ facet_wrap(~taxonomic_class) + geom_abline(intercept = 0, slope = 1) + xlab(expression(nu~with~AR1)) + ylab(expression(nu~without~AR1)) + coord_fixed() + theme_bw() +  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.333, 0.5), labels = c("Infinity", "10", "5", "3", "2")) +  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.333, 0.5), labels = c("Infinity", "10", "5", "3", "2")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + labs(colour = expression(phi)) + scale_colour_gradient2(low = scales::muted("blue"), high = scales::muted("red"), mid = "grey")
ggsave("effect-of-ar1-on-nu-by-class.pdf", width = 7, height = 7)

p <- ggplot(ar1_vs_base, aes(1/nu_50, 1/nu_50_base, col = phi_50)) + geom_point(cex = 3.5, alpha = 0.8) + geom_segment(aes(x = 1/nu_50, xend = 1/nu_50, y = 1/nu_25_base, yend = 1/nu_75_base), lwd = 0.5, alpha = 0.2) + geom_segment(aes(x = 1/nu_25, xend = 1/nu_75, y = 1/nu_50_base, yend = 1/nu_50_base), lwd = 0.5, alpha = 0.2)+  geom_abline(intercept = 0, slope = 1) + xlab(expression(nu~with~AR1)) + ylab(expression(nu~without~AR1)) + coord_fixed() + theme_bw() +  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.333, 0.5), labels = c("Infinity", "10", "5", "3", "2")) +  scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.333, 0.5), labels = c("Infinity", "10", "5", "3", "2")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + labs(colour = expression(phi)) + scale_colour_gradient2(low = scales::muted("blue"), high = scales::muted("red"), mid = "grey")
ggsave("effect-of-ar1-on-nu.pdf", width = 7, height = 7)

# and the effect of assumed observation error?
p <- ggplot(ar1_no_obs_vs_obs, aes(1/nu_50, 1/nu_50_ar1_no_obs, colour = taxonomic_order)) + geom_point() + geom_abline(intercept = 0, slope = 1) + xlab(expression(1/nu~with~assumed~observation~error)) + ylab(expression(1/nu~without~observation~error)) + coord_fixed() + theme_bw()
ggsave("effect-of-obs-on-nu.pdf", width = 7, height = 7)

mytheme <- theme(panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(), panel.background = element_blank())

# no_fur <- filter(gomp_hat_ar1, fur == FALSE) %>%
#   arrange(taxonomic_class, nu_50) %>%
#   group_by(taxonomic_class) %>%
#   mutate(sort_id = seq_along(main_id),
#     sort_id_perc = 100 * (sort_id / max(sort_id)))

p <- ggplot(filter(gomp_hat_base, !taxonomic_class %in% c("Crustacea", "Gastropoda", "Chondrichtyhes")),
  aes(sort_id_perc, 1/nu_50)) +
  facet_wrap(~taxonomic_class) +
  geom_linerange(aes(x = sort_id_perc,
    ymax = 1/nu_95, ymin = 1/nu_5), alpha = 0.3, width = 0.02) +
  geom_linerange(aes(x = sort_id_perc,
    ymax = 1/nu_75, ymin = 1/nu_25), alpha = 0.8, width = 0.15) +
  geom_point(aes(x = sort_id_perc, y = 1/nu_50), col = "black") +
  theme_bw() +
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.333, 0.5),
    labels = c("Infinity", "10", "5", "3", "2")) +
  ylab("Estimated t degrees of freedom") + xlab("Percent of populations") +
 mytheme
 p

ggsave("gomp-nu-base.pdf", width = 6.5, height = 6)

# filter(gomp_hat_ar1, taxonomic_class == "Mammalia") %>%
#   ggplot(aes(sort_id_perc, 1/nu_50, colour = taxonomic_order)) +
#   facet_wrap(~taxonomic_order) +
#   geom_linerange(aes(x = sort_id_perc,
#     ymax = 1/nu_95, ymin = 1/nu_5), alpha = 0.2, width = 0.01) +
#   geom_linerange(aes(x = sort_id_perc,
#     ymax = 1/nu_75, ymin = 1/nu_25), alpha = 0.2, width = 0.02) +
#   geom_point(aes(x = sort_id_perc, y = 1/nu_50), col = "black") +
#   theme_bw() +
#   scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.333, 0.5),
#     labels = c("Infinity", "10", "5", "3", "2")) +
#   ylab("Estimated t degrees of freedom") + xlab("Percent of populations") +
#   mytheme


p <- ggplot(subset(gomp_hat_base, !taxonomic_class %in% c("Crustacea", "Gastropoda", "Chondrichtyhes")), aes(p10, taxonomic_order)) + geom_point(position = position_jitter(height = 0.18), alpha = 0.3, colour = "black") + ylab("Taxonomic order") + xlab(expression(Probability~nu<10~(heavy~tails))) + theme_bw() + facet_wrap(~taxonomic_class, scales = "free_y")
ggsave("gomp-base-p10-dot-order.pdf", width = 9, height = 5)

p <- ggplot(gomp_hat_ar1, aes(p10)) + geom_histogram(binwidth = 0.02) + ylab("Frequency") + xlab(expression(Probability~nu<10~(probability~of~heavy~tails))) + theme_bw() +
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  facet_wrap(~taxonomic_class, scales = "free_y")
ggsave("gomp-ar1-p10-hist.pdf", width = 6, height = 4)

inv_nu_y_labs <- scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.333, 0.5),
  labels = c("Infinity", "10", "5", "3", "2"))

p <- ggplot(gomp_hat_ar1_long, aes(value, 1/nu_50, colour = taxonomic_class)) + geom_point() + inv_nu_y_labs + ylab(expression(nu)) + facet_wrap(~variable, scales = "free_x") + xlab("Covariate value") + labs(colour = "Class")  + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_linerange(aes(x = value, ymax = 1/nu_25, ymin = 1/nu_75), alpha = 0.2, width = 0.05)

p <- ggplot(subset(gomp_hat_base_long, variable != "log_adult_body_mass_g"), aes(value, 1/nu_50, colour = taxonomic_class)) + geom_point() + inv_nu_y_labs + ylab(expression(nu)) + facet_wrap(~variable, scales = "free_x") + xlab("Covariate value") + labs(colour = "Class")  + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_linerange(aes(x = value, ymax = 1/nu_25, ymin = 1/nu_75), alpha = 0.2, width = 0.05)

ggsave("nu-covariates-base.pdf", width = 8, height = 6)

# which ones have a very high lambda_50? pink salmon:
#  ggplot(gpdd, subset(main_id %in% subset(gomp_hat_base, lambda_50 > 12)$main_id), aes(series_step, population_untransformed)) + geom_line() + facet_wrap(~label, scales = "free")

p <- ggplot(subset(gomp_hat_base_long_p, variable != "log_adult_body_mass_g"), aes(value, p10, colour = taxonomic_class)) + geom_point(alpha = 0.6) + ylab(expression(p(nu<10))) + facet_wrap(~variable, scales = "free_x") + xlab("Covariate value") + labs(colour = "Class")  + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

ggsave("nu-covariates-p-base.pdf", width = 8, height = 6)

p1 <- ggplot(filter(gomp_hat_ar1, taxonomic_class == "Mammalia"), aes(log10(adult_body_mass_g), p10, colour = taxonomic_order)) + geom_point() + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab("Probability nu < 10") + labs(colour = "Taxonomic order")

p2 <- ggplot(filter(gomp_hat_ar1, taxonomic_class == "Mammalia"), aes(log10(sexual_maturity_age_d/365), p10, colour = taxonomic_order)) + geom_point() + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ylab("Probability nu < 10") + labs(colour = "Taxonomic order")

pdf("p10-mammals-cross.pdf", width = 7, height = 10)
gridExtra::grid.arrange(p1, p2)
dev.off()

ticks <- c(2, 5, 20, 50, 100)

p <- ggplot(gomp_vs_logistic, aes(log(1/nu_50), log(1/logistic_nu_50), colour = taxonomic_class)) + geom_point(alpha = 0.8) +geom_segment(aes(x = log(1/nu_25), xend = log(1/nu_75), y = log(1/logistic_nu_50), yend = log(1/logistic_nu_50)), alpha = 0.2) + geom_segment(aes(x = log(1/nu_50), xend = log(1/nu_50), y = log(1/logistic_nu_25), yend = log(1/logistic_nu_75)), alpha = 0.2) +  geom_abline(intercept = 0, slope = 1) + scale_y_continuous(breaks = log(1/ticks), labels = ticks) + scale_x_continuous(breaks = log(1/ticks), labels = ticks) + xlab("Gompertz nu estimate") + ylab("Ricker-logistic nu estimate")

ticks <- c(2, 3, 5, 10, Inf)

p <- ggplot(gomp_vs_logistic, aes((1/nu_50), (1/logistic_nu_50), colour = taxonomic_class)) + geom_point(cex = 3 , alpha = 0.8) +geom_segment(aes(x = (1/nu_25), xend = (1/nu_75), y = (1/logistic_nu_50), yend = (1/logistic_nu_50)), alpha = 0.3) + geom_segment(aes(x = (1/nu_50), xend = (1/nu_50), y = (1/logistic_nu_25), yend = (1/logistic_nu_75)), alpha = 0.3) +  geom_abline(intercept = 0, slope = 1) + scale_y_continuous(breaks = (1/ticks), labels = ticks, limits = c(0, 0.5)) + scale_x_continuous(breaks = (1/ticks), labels = ticks, limits = c(0, 0.5)) + xlab("Gompertz nu estimate") + ylab("Ricker-logistic nu estimate") + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) + coord_fixed()

ggsave("gomp-vs-logistic-nu.pdf", width = 7, height = 5)

#ggplot(gomp_vs_logistic, aes(1/nu_50, 1/logistic_nu_50, colour = taxonomic_class)) + geom_point() +  geom_abline(intercept = 0, slope = 1) + scale_y_continuous(breaks = 1/ticks, labels = ticks)

#p <- ggplot(gomp_vs_logistic, aes(log(sigma_proc_50), log(logistic_sigma_proc_50), colour = taxonomic_class)) + geom_point() + geom_abline(intercept = 0, slope = 1)
#

# plot Brook et al. 2006 characteristics against nu estimates:

p <- ggplot(gomp_hat_base_long_brooks, aes(value, 1/nu_50, colour = taxonomic_class)) + geom_point() + inv_nu_y_labs + ylab(expression(nu)) + facet_wrap(~variable, scales = "free_x") + xlab("Covariate value") + labs(colour = "Class")  + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + geom_linerange(aes(x = value, ymax = 1/nu_25, ymin = 1/nu_75), alpha = 0.2, width = 0.05)
ggsave("brook-etal-correlates.pdf", width = 8, height = 6)
