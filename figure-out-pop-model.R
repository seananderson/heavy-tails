# play with gompertz and logistic growth model parameterizations
# figure out what's the same and what's different and how they
# model growth

logistic1 <- function(Nt, r, K, theta = 1) {
  r * (1 - (Nt/K)^theta)
}

#logistic2 <- function(Nt, r, K) {
  #r * Nt * (1 - Nt/K)
#}

# both return log(abundance_{t + 1})
logistic3 <- function(Nt, r, K) {
  log(Nt) + r * (1 - Nt/K)
}

gompertz1 <- function(Nt, lambda, b) {
  lambda + b * log(Nt)
}

gompertz2 <- function(Nt, lambda, b) {
  log(Nt) + lambda + b * log(Nt)
}

N <- 50
logx <- log(seq(1, 20))
r <- 0.8
K <- 10

par(mfrow = c(2, 1))
plot(exp(logx), log(exp(logistic3(exp(logx), 0.9, 5.29))/exp(logx)))
plot(logx, log(exp(gompertz1(exp(logx), 0.5, 0.6))/exp(logx)))
#plot(logx, log(exp(gompertz2(exp(logx), 0.5, 0.6-1))/exp(logx)))

#exp(1.0/0.6)

par(mfrow = c(2, 1))
plot(exp(logx), log(exp(logistic3(exp(logx), 0.5, 5.29))/exp(logx)), main = "logistic")
plot(exp(logx), log(exp(gompertz1(exp(logx), 0.5, 0.6))/exp(logx)), main = "gompertz")
