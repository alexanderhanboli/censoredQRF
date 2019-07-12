library(quantreg)
library(survival)

rho <- function(u, tau) {
  u*(tau - (u<0))
}

# loss <- function(Y, delta, tau, lambda) {
#   Y.max <- max(Y)
#   Y.min <- min(Y)
#   n <- length(Y)
#   censoring_prob <- sum(1-delta)/n
#   C.M <- Y.max + 100
#   C.m <- Y.min - 100
#   loss1 <- (1-delta)*rho(Y - pmin(Y, lambda) ,tau)
#   loss2 <- delta*(censoring_prob*rho(Y - pmin(C.m, lambda), tau) + (1-censoring_prob)*rho(Y - pmin(C.M, lambda), tau))
#   loss <- loss1 + loss2
#   return (sum(loss))
# }

# loss <- function(Y, delta, tau, lambda) {
#   Y.max <- max(Y)
#   Y.min <- min(Y)
#   n <- length(Y)
#   censoring_prob <- sum(1-delta)/n
#   kappa <- 1 - pnorm(tau, mean = 1-censoring_prob, sd = sqrt(censoring_prob*(1-censoring_prob)/1))
#   loss1 <- (1-delta)*rho(Y - pmin(Y, lambda), tau)
#   loss2 <- kappa*delta*rho(Y - lambda, tau)
#   loss <- loss1 + loss2
#   return (sum(loss))
# }

loss <- function(Y, delta, tau, lambda, C.surv) {
  Y.max <- max(Y[delta==1])
  Y.min <- min(Y[delta==1])
  n <- length(Y)
  kappa <- C.surv(lambda)
  # loss1 <- (1-delta)*1*(Y > lambda)*tau
  # loss2 <- delta*1*(Y > lambda)*tau + delta*1*(Y < lambda)*(tau-1)*kappa
  # loss <- loss1 + loss2
  loss <- 1*(Y > lambda) - kappa*(1-tau)
  # print(kappa)
  return ((sum(loss))^2)
}

latent_loss <- function(Y, tau, lambda) {
  input <- Y - lambda
  return (sum(rho(input, tau)))
}

oracle_loss <- function(Y, C, tau, lambda) {
  input <- Y - pmin(C, lambda)
  return (sum(rho(input, tau)))
}

n <- 2000
Latent <- runif(n)
#b <- rbinom(n, 1, 0.8)
#C <- b*rnorm(n, mean = 0.1, sd = 5) + (1-b)*rnorm(n, mean = 1.9, sd = 5)
#C <- 5*rbinom(n, 1, 0.2) + runif(n)
# C <- 0.8
C <- rnorm(n, mean = 0.3)
Y <- pmin(Latent, C)
delta <- 1*(Latent <= C)
lambda <- seq(min(Latent), max(Latent), length.out = 100)
C.km <- survfit(Surv(Y, 1-delta) ~ 1, type = 'kaplan-meier')
C.surv <- stepfun(C.km$time, c(1, C.km$surv))

findmin <- function(lambda, Y, Latent, delta, tau, C.surv, option = 'approx') {
  if (option == 'approx') {
    l <- sapply(lambda, function (x) {loss(Y, delta, tau, x, C.surv)})
    return (lambda[which(l == min(l))[1]])
  } else if (option == 'truth') {
    l <- sapply(lambda, function (x) {latent_loss(Latent, tau, x)})
    return (lambda[which(l == min(l))[1]])
  } else if (option == 'oracle') {
    l <- sapply(lambda, function (x) {oracle_loss(Y, C, tau, x)})
    return (lambda[which(l == min(l))[1]])
  } else {
    l <- sapply(lambda, function (x) {latent_loss(Y, tau, x)})
    return (lambda[which(l == min(l))[1]])
  }
}

sum(delta)/length(delta)
tau_list <- seq(0,1,length.out = 101)
plot(tau_list, sapply(tau_list, function (x) {findmin(lambda, Y, Latent, delta, x, C.surv, option = 'approx')}), 
     pch=17, xlab = 'tau', ylab = 'quantile', ylim = c(0,1), cex = 0.5)
lines(tau_list, sapply(tau_list, function (x) {findmin(lambda, Y, Latent, delta, x, C.surv, option = 'biased')}), 
       col = 3, type = 'l')
lines(tau_list, sapply(tau_list, function (x) {findmin(lambda, Y, Latent, delta, x, C.surv, option = 'truth')}), 
       col = 2, type = 'l')
# lines(tau_list, sapply(tau_list, function (x) {findmin(lambda, Y, Latent, delta, x, C.surv, option = 'oracle')}), 
#        col = 4, type = 'l')
abline(a = 0, b = 1, lty = 2)

# KM
# fit.km <- survfit(Surv(Y, delta) ~ 1, type = 'kaplan-meier')
# lines(tau_list, quantile(fit.km, probs = tau_list, conf.int=FALSE), col = 6)
