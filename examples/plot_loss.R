library(survival)

rho <- function(u, tau) {
  u*(tau - (u<0))
}

loss <- function(Y, tau, q) {
  input <- Y - q
  return (mean(rho(input, tau)))
}

oracle_loss <- function(Y, C, tau, q) {
  input <- Y - pmin(C, q)
  return (mean(rho(input, tau)))
}

attach(mtcars)
op <- par(mfrow = c(2,2),
          oma = c(5,4,0,0) + 0.1,
          mar = c(1,1,1,1) + 0.1)

for (n in c(100, 500, 1000, 5000)) {
  Latent <- runif(n)
  C <- rnorm(n, mean = 0.8, sd = 0.2)
  Y <- pmin(Latent, C)
  delta <- 1*(Latent <= C)
  sum(delta)/length(delta)
  tau <- 0.5
  q <- seq(0, 1, length.out = 1000)
  
  # survival estimate of C
  C.km <- survfit(Surv(C, 1-delta) ~ 1, type = 'kaplan-meier')
  C.surv <- stepfun(C.km$time, c(1, C.km$surv))
  
  plot(q, sapply(q, function (x) {loss(Y, tau, x)}), type='l', ylab = 'loss', col = 1)
  lines(q, sapply(q, function (x) {loss(Y, 1 - (1-tau)*C.surv(x), x)}), col = 2)
  
  min1 <- min(sapply(q, function (x) {loss(Y, tau, x)}))
  argmin1 <- q[which(sapply(q, function (x) {loss(Y, tau, x)}) == min1)]
  min2 <- min(sapply(q, function (x) {loss(Y, 1 - (1-tau)*C.surv(x), x)}))
  argmin2 <- q[which(sapply(q, function (x) {loss(Y, 1 - (1-tau)*C.surv(x), x)}) == min2)]
  
  abline(v = tau, col = 4)
  abline(v = argmin1, col = 1, lty = 2)
  abline(v = argmin2, col = 2, lty = 2)
  # Add a legend
  # legend(0.35, 0.25, legend=c("truth", "our method"), lty = 1, col = c(1,2)) 
}
title(xlab = "q",
      ylab = "sample quantile loss",
      outer = TRUE, line = 3)
par(op)