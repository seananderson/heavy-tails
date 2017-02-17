# take the output from 3-test-t-sampling.R and plot it

sample_t <- readRDS("sample-t-sim-check.rds")

library("ggplot2")
theme_set(theme_bw())
library("dplyr")

sample_t <- filter(sample_t, nu_true > 1)
p <- ggplot(sample_t, aes(iter, 1/med_nu)) + geom_pointrange(aes(ymin = 1/l_nu, ymax = 1/u_nu), size = 0.2) +
  facet_grid(nu_true~N) +
  scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.333, 0.5),
    labels = c("Infinity", "10", "5", "3", "2")) +
  geom_hline(aes(yintercept = 1/nu_true), col = "red") +
  xlab("Iteration") + ylab(expression(t-distribution~degrees~of~freedom~(nu)))

  ggsave("t-dist-sampling-sim-prior-exp0point01.pdf", width = 8, height = 5)

p <- ggplot(sample_t, aes(iter, med_sigma)) + geom_pointrange(aes(ymin = l_sigma, ymax = u_sigma), size = 0.2) +
  facet_grid(nu_true~N) +
  geom_hline(aes(yintercept = 1), col = "red") +
  xlab("Iteration") + ylab("Scale parameter")

  ggsave("t-dist-sampling-sim-sigma-prior-exp0point01.pdf", width = 8, height = 5)

#sample_t <- readRDS("sample-t-sim-check-prior-exp-0.05.rds")
#p <- ggplot(sample_t, aes(iter, 1/med_nu)) + geom_pointrange(aes(ymin = 1/l_nu, ymax = 1/u_nu)) +
  #facet_grid(nu_true~N) +
  #scale_y_continuous(breaks = c(0, 0.1, 0.2, 0.333, 0.5),
    #labels = c("Infinity", "10", "5", "3", "2")) +
  #geom_hline(aes(yintercept = 1/nu_true), col = "red") +
  #xlab("Iteration") + ylab(expression(t-distribution~degrees~of~freedom~(nu)))
  #ggsave("t-dist-sampling-sim-prior-exp0.05.pdf", width = 8, height = 5)
