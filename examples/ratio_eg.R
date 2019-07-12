n <- 5000
Ti <- runif(n)
Ci <- rnorm(n, 0.8, 0.5)
Yi <- pmin(Ti, Ci)
delta <- 1*(Ti <= Ci)
# lambda <- sort(Yi)
lambda <- seq(min(Yi), max(Yi), length.out = 1000)
censor_ratio <- mean(1-delta)

sum1 <- sapply(lambda, function (x) {sum(1*(Ci >= x) * (1 - delta))})
sum2 <- sapply(lambda, function (x) {sum(1*(Ci >= x) * delta)})
plot(lambda, (sum1/sum2)^0.5)
abline(a = (censor_ratio/(1-censor_ratio))^0.5, b = 0, col = 'red')