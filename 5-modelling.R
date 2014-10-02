# try some figure making

library(ggplot2)
library(mgcv)
newdata <- data.frame(log10_Len = seq(min(gomp_hat_base$log10_Len, na.rm = TRUE), max(gomp_hat_base$log10_Len, na.rm = TRUE), length.out = 100), log10_Lifesp = seq(min(gomp_hat_base$log10_Lifesp, na.rm = TRUE), max(gomp_hat_base$log10_Lifesp, na.rm = TRUE), length.out = 100))

m1 <- gam(1/nu_50 ~ s(log10_Len), data = gomp_hat_base)

ggplot(gomp_hat_base, aes(log10_Len, p10, colour = taxonomic_class)) + geom_point()

m2 <- gam(p10 ~ s(log10_Len), data = gomp_hat_base, family = betar(link = "logit"))
plot(m2)
newdata$p <- predict(m2, newdata = newdata, type = "response")
plot(newdata$log10_Len, newdata$p)

m3 <- gamm(p10 ~ s(log10_Len), data = gomp_hat_base, family = betar(link = "logit"), random = list(taxonomic_class=~1, taxonomic_order=~1), method = "REML")
plot(m3$gam)

ggplot(gomp_hat_base, aes(log10_Lifesp, p10, colour = taxonomic_class)) + geom_point()
m4 <- gamm(p10 ~ s(log10_Lifesp), data = gomp_hat_base, family = betar(link = "logit"), random = list(taxonomic_class=~1, taxonomic_order=~1), method = "REML")
plot(m4$gam)
newdata$log10_Lifesp_pred <- predict(m4$gam, newdata = newdata, type = "response")
with(newdata, plot(log10_Lifesp, log10_Lifesp_pred))

m5 <- glm(p10 ~ log10_Lifesp, data = gomp_hat_base, family = quasibinomial)
newdata$m5 <- predict(m5, newdata = newdata, type = "response")
plot(newdata$log10_Lifesp, newdata$m5)


library(gbm)
m6 <- gbm(p10~log10_Lifesp + as.factor(taxonomic_class), data = gomp_hat_base)

plot.gbm(m6, i.var = 1)

m7 <- gbm(p10~b_50 + as.factor(taxonomic_class), data = gomp_hat_base)
plot(m7)

#library(gamm4)
#m4 <- gamm4(p10 ~ s(log10_Len), data = gomp_hat_base, family = betar(link = "logit"), random = ~ (1|taxonomic_class/taxonomic_order))

gomp_hat_base$heavy <- ifelse(gomp_hat_base$p10 > 0.5, 1, 0)
with(gomp_hat_base, plot(b_50, jitter(heavy, 0.1)))

library(lme4)

gomp_hat_base$log10_Lifesp_scaled <- scale(gomp_hat_base$log10_Lifesp)
gomp_hat_base$sigma_proc_50_scaled <- scale(gomp_hat_base$sigma_proc_50)
gomp_hat_base$b_50_scaled <- scale(gomp_hat_base$b_50)
gomp_hat_base$dataset_length_scaled <- scale(gomp_hat_base$dataset_length)
gomp_hat_base$log10_Len_scaled <- scale(gomp_hat_base$log10_Len)

m8 <- glmer(heavy ~ log10_Lifesp_cent + (1 | taxonomic_class / taxonomic_order), family = binomial, data = gomp_hat_base)
summary(m8)
#newdata2 <- data.frame(
#newdata <- data.frame(log10_Len = seq(min(gomp_hat_base$log10_Len, na.rm = TRUE), max(gomp_hat_base$log10_Len, na.rm = TRUE), length.out = 100), log10_Lifesp = seq(min(gomp_hat_base$log10_Lifesp, na.rm = TRUE), max(gomp_hat_base$log10_Lifesp, na.rm = TRUE), length.out = 100))
#newdata$

newdata$m8 <- predict(m8, newdata = newdata, type = "response")

m9 <- glmer(heavy ~ b_50 + (1 | taxonomic_class / taxonomic_order), family = binomial, data = gomp_hat_base)
summary(m9)

m10 <- glmer(heavy ~ lambda_50 + (1 | taxonomic_class / taxonomic_order), family = binomial, data = gomp_hat_base)
summary(m10)

m11 <- glmer(heavy ~ log10_Len + (1 | taxonomic_class / taxonomic_order), family = binomial, data = gomp_hat_base)
summary(m11)
