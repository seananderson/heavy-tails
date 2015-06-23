# consider increase in probability of extinction varlues of nu


N <- 100
nu <- 9999
sigma_proc <- 0.7
lambda <- 1
b <- 0.7
y1 <- 3
burnin <- 1:20

y <- vector(length = N, mode = "numeric")
y[1] <- 3
proc_error <- metRology::rt.scaled(N, df = nu, mean = 0, sd = sigma_proc)
for(i in 2:N) {
  y[i] <- lambda + b * y[i-1] + proc_error[i-1]
}
y <- y[-burnin]
plot(y, type = "o", ylim = c(0, 10))
#print(mean(y))
abline(h = 3.3, col = "red")
abline(h = 3.3 * 0.7, col = "blue")
