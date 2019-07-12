library(survival)

est <- function(Y, tau, q) {
  input <- 1*(Y > q)
  return (mean((1 - tau) - input))
}

new_est <- function(Y, tau, q, G) {
  input <- 1*(Y > q)
  return (mean((1 - tau)*G(q) - input))
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
  
  plot(q, sapply(q, function (x) {est(Y, tau, x)}), type='l', ylab = 'est', col = 1)
  lines(q, sapply(q, function (x) {new_est(Y, tau, x, C.surv)}), col = 2)
  
  min1 <- min(sapply(q, function (x) {est(Y, tau, x)})^2)
  root1 <- q[which(sapply(q, function (x) {est(Y, tau, x)})^2 == min1)]
  
  min2 <- min(sapply(q, function (x) {new_est(Y, tau, x, C.surv)})^2)
  root2 <- q[which(sapply(q, function (x) {new_est(Y, tau, x, C.surv)})^2 == min2)]
  
  abline(h = 0, lty = 3)
  abline(v = tau, col = 4)
  abline(v = root1, col = 1, lty = 2)
  abline(v = root2, col = 2, lty = 2)
  # Add a legend
  # legend(0.35, 0.25, legend=c("truth", "our method"), lty = 1, col = c(1,2)) 
}
title(xlab = "q",
      ylab = "U(q)",
      outer = TRUE, line = 3)
par(op)