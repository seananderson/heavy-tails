library(lme4)

#gomp_hat_base$log10_Lifesp_scaled <- scale(gomp_hat_base$log10_Lifesp)
gomp_hat_base$log_Lifesp_scaled <- scale(log(gomp_hat_base$Lifesp))
gomp_hat_base$log_sigma_proc_50_scaled <- scale(log(gomp_hat_base$sigma_proc_50))
gomp_hat_base$b_50_scaled <- scale(gomp_hat_base$b_50)
gomp_hat_base$lambda_50_scaled <- scale(gomp_hat_base$lambda_50)
gomp_hat_base$log_dataset_length_scaled <- scale(log(gomp_hat_base$dataset_length))
gomp_hat_base$log_Len_scaled <- scale(log(gomp_hat_base$Len))
gomp_hat_base$heavy <- ifelse(gomp_hat_base$p10 > 0.5, 1, 0)

# models <<- list()
# i <<- 0
# get_lmer_coefs <- function(param) {
#   i <<- i + 1
#   m1 <- glmer(heavy ~ gomp_hat_base[, param] + (1 | taxonomic_class / taxonomic_order), family = binomial, data = gomp_hat_base)
#   models[[i]] <- m1
#   est <- fixef(m1)[[2]]
#   ci.50 <- tryCatch(confint(m1, parm = 4, level = 0.5), warning = function(w) c(NA, NA))
#   ci.90 <- tryCatch(confint(m1, parm = 4, level = 0.90), warning = function(w) c(NA, NA))
#   ret <- data.frame(est = est, l.50 = ci.50[1], u.50 = ci.50[2], l.90 = ci.90[1],
#     u.90 = ci.90[2])
#   print(ret)
#   ret
# }
#
# out <- data.frame(predictor =
#   names(gomp_hat_base)[grep("scaled", names(gomp_hat_base))],
#   stringsAsFactors = FALSE)
#
# out <- plyr::ddply(out, "predictor", function(x) {
#   print(x$predictor)
#   get_lmer_coefs(x$predictor)
# })
#
# out$pos <- 1:6
#
# par(mar = c(4, 10, .5, .5), cex = 0.8)
# plot(1, 1, xlim = c(-1, 1), ylim = c(1,6), type = "n", axes = FALSE,
#   xlab = "Scaled coefficient", ylab = "")
# with(out, segments(l.90, pos, u.90, pos))
# with(out, segments(l.50, pos, u.50, pos, lwd = 2))
# with(out, points(est, pos, pch = 21, bg = "grey40", col = "grey10"))
# box()
# axis(1)
# axis(2, at = out$pos, labels = out$predictor, las = 2)
# abline(v = 0, lty = 2)

#m.all <- glmer(heavy ~ log_Lifesp_scaled + log_sigma_proc_50_scaled + log_dataset_length_scaled + b_50_scaled + lambda_50_scaled + (1 | taxonomic_class / taxonomic_order), family = binomial, data = gomp_hat_base)

library(glmmADMB)
temp.dat <- gomp_hat_base[,c("log_Lifesp_scaled", "log_sigma_proc_50_scaled", "log_dataset_length_scaled", "b_50_scaled", "lambda_50_scaled", "taxonomic_class", "taxonomic_order", "heavy", "p10")]
temp.dat <- na.omit(temp.dat)
temp.dat$taxonomic_class <- as.factor(temp.dat$taxonomic_class)
temp.dat$taxonomic_order <- as.factor(temp.dat$taxonomic_order)
m.admb.1 <- glmmadmb(heavy ~ log_Lifesp_scaled + log_sigma_proc_50_scaled + log_dataset_length_scaled + b_50_scaled + lambda_50_scaled + (1 | taxonomic_class / taxonomic_order), family = "binomial", data = temp.dat)
save(m.admb.1, file = "m.admb.1.rda")

m.admb.2 <- glmmadmb(p10 ~ log_Lifesp_scaled + log_sigma_proc_50_scaled + log_dataset_length_scaled + b_50_scaled + lambda_50_scaled + (1 | taxonomic_class / taxonomic_order), family = "beta", data = temp.dat)
save(m.admb.2, file = "m.admb.2.rds")

temp.dat2 <- gomp_hat_base[,c("log_Lifesp_scaled", "log_sigma_proc_50_scaled", "dataset_length", "log_dataset_length_scaled", "b_50_scaled", "lambda_50_scaled", "taxonomic_class", "taxonomic_order", "heavy", "p10")]
temp.dat2 <- na.omit(temp.dat2)
temp.dat2$taxonomic_class <- as.factor(temp.dat2$taxonomic_class)
temp.dat2$taxonomic_order <- as.factor(temp.dat2$taxonomic_order)

m.admb.3 <- glmmadmb(heavy ~ log_Lifesp_scaled + log_sigma_proc_50_scaled + dataset_length + b_50_scaled + lambda_50_scaled + (1 | taxonomic_class / taxonomic_order), family = "binomial", data = temp.dat2)

m.admb.4 <- glmmadmb(p10 ~ log_Lifesp_scaled + log_sigma_proc_50_scaled + dataset_length + b_50_scaled + lambda_50_scaled + (1 | taxonomic_class / taxonomic_order), family = "beta", data = temp.dat2)

newdat <- data.frame(log_Lifesp_scaled = 0, log_sigma_proc_50_scaled = 0, dataset_length = seq(20, 100), b_50_scaled = 0, lambda_50_scaled = 0)

p <- predict(m.admb.3, newdata = newdat, se = TRUE)
x <- newdat$dataset_length
p.beta <- predict(m.admb.4, newdata = newdat, se = TRUE)

pdf("length-vs-prob-tails.pdf", width = 7, height = 4)
par(mfrow = c(1, 2), las = 1)
par(mar = c(4, 1, 2, 1), cex = 0.8, oma = c(1, 4, 1, 1))

plot(x, plogis(p.beta$fit), type = "l", ylab = "", xlab = "Time-series length", ylim = c(0, 1.01), log = "x", yaxs = "i")
with(temp.dat2, points(dataset_length, p10, pch = 21, col = "#00000070", bg = "#00000030"))
polygon(c(x, rev(x)), c(plogis(p.beta$fit + qnorm(0.75) * p.beta$se.fit), rev(plogis(p.beta$fit - qnorm(0.75) * p.beta$se.fit))), border = NA, col = "#FF000020")
polygon(c(x, rev(x)), c(plogis(p.beta$fit + qnorm(0.975) * p.beta$se.fit), rev(plogis(p.beta$fit - qnorm(0.975) * p.beta$se.fit))), border = NA, col = "#FF000020")
par(xpd = NA)
mtext("Probability of\nobserving heavy tails", side = 2, line = 2.5, cex = 0.8, las = 0)
mtext("Beta model", side = 3, line = 1, cex = 0.8, las = 1)
par(xpd = FALSE)

plot(x, plogis(p$fit), type = "l", ylab = "", xlab = "Time-series length", ylim = c(0, 1), log = "x", yaxt = "n")
with(temp.dat2, points(dataset_length, jitter(heavy, 0.1), pch = 21, col = "#00000070", bg = "#00000030"))
polygon(c(x, rev(x)), c(plogis(p$fit + qnorm(0.75) * p$se.fit), rev(plogis(p$fit - qnorm(0.75) * p$se.fit))), border = NA, col = "#FF000020")
polygon(c(x, rev(x)), c(plogis(p$fit + qnorm(0.975) * p$se.fit), rev(plogis(p$fit - qnorm(0.975) * p$se.fit))), border = NA, col = "#FF000020")
mtext("Binomial model", side = 3, line = 1, cex = 0.8, las = 1)
dev.off()





#est.lmer.bin <- fixef(m.all)
est.admb.bin <- fixef(m.admb.1)
est.admb.beta <- fixef(m.admb.2)
#ci.lmer.bin <- confint(m.all)
ci.admb.bin.95 <- confint(m.admb.1)
ci.admb.beta.95 <- confint(m.admb.2)
ci.admb.bin.50 <- confint(m.admb.1, level = 0.5)
ci.admb.beta.50 <- confint(m.admb.2, level = 0.5)

ord <- rev(order(est.admb.beta[-1]))

make_coef_plot <- function(cis1, cis2, fe) {
  xlim <- range(rbind(cis1[-1, ], cis2[-1, ]))

  cis1 <- cis1[-1, ][ord, ]
  cis2 <- cis2[-1, ][ord, ]
  fe <- fe[-1][ord]

  plot(1, 1, xlim = xlim, ylim = c(1,5), type = "n", axes = FALSE,
    xlab = "", ylab = "")
  pos <- 1:5
  for(i in 1:5) {
    segments(cis1[i, 1], pos[i], cis1[i, 2], pos[i], lwd = 2.5)
    segments(cis2[i, 1], pos[i], cis2[i, 2], pos[i], lwd = 0.6)
    points(fe[i], pos[i], pch = 21, bg = "grey40", col = "grey10")
  }
  box()
  axis(1)
  abline(v = 0, lty = 2)
}

pdf("admb-coefs.pdf", width = 6.5, height = 4)
par(mar = c(4, .5, .5, .5), cex = 0.8, mfrow = c(1, 2), oma = c(.5, 11.5, 0, .5))
par(tck = -0.02, mgp = c(2, 0.5, 0), col.axis = "grey25", col = "grey25")
make_coef_plot(ci.admb.beta.50, ci.admb.beta.95, est.admb.beta)
mtext("Beta model\ncoefficient", line = 3, side = 1)
axis(2, at = 1:5, labels = names(est.admb.bin)[-1][ord], las = 2)
make_coef_plot(ci.admb.bin.50, ci.admb.bin.95, est.admb.bin)
mtext("Binomial model\ncoefficient", side = 1, line = 3)
dev.off()

#m.all2 <- glmer(heavy ~ log10_Len_scaled + log_sigma_proc_50_scaled + log_dataset_length_scaled + b_50_scaled + lambda_50_scaled + (1 | taxonomic_class / taxonomic_order), family = binomial, data = gomp_hat_base)
#
# test for overdispersion, function from http://glmm.wikidot.com/faq
#
overdisp_fun <- function(model) {
  ## number of variance parameters in
  ##   an n-by-n variance-covariance matrix
  vpars <- function(m) {
    nrow(m)*(nrow(m)+1)/2
  }
  model.df <- sum(sapply(VarCorr(model),vpars))+length(fixef(model))
  rdf <- nrow(model.frame(model))-model.df
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
overdisp_fun(m.admb.1)
overdisp_fun(m.admb.2)
#overdisp_fun(m.all)

#est <- fixef(m.all)
#ci.95 <- confint(m.all, level = 0.95)

#ci.50 <- confint(m.all, level = 0.50)

# get_lmer_coefs("sigma_proc_50_scaled")
# get_lmer_coefs("b_50_scaled")
# get_lmer_coefs("lambda_50_scaled")
# get_lmer_coefs("dataset_length_scaled")
# get_lmer_coefs("log10_Len_scaled")
