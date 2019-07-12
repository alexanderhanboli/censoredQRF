rho <- function(u, tau) {
  u*(tau - (u<0))
}

loss <- function(Y, delta, tau, lambda) {
  input <- Y - (1-delta)*pmin(Y, lambda) - delta*lambda
  return (sum(rho(input, tau)))
}

latent_loss <- function(Y, tau, lambda) {
  input <- Y - lambda
  return (sum(rho(input, tau)))
}

oracle_loss <- function(Y, C, tau, lambda) {
  input <- Y - pmin(C, lambda)
  return (sum(rho(input, tau)))
}

n <- 1000
Latent <- runif(n)
# C <- rnorm(n, mean = 0.8, sd = 0.2)
C <- 1*rbinom(n, 1, 0.4) + runif(n)
Y <- pmin(Latent, C)
delta <- 1*(Latent <= C)
tau <- 0.6
lambda <- seq(0, 1, length.out = 100)

findmin <- function(lambda, Y, Latent, delta, tau, option = 'approx') {
  if (option == 'approx') {
    l <- sapply(lambda, function (x) {loss(Y, delta, tau, x)})
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
tau_list <- seq(0,1,length.out = 11)
plot(tau_list, 
     sapply(tau_list, function (x) {findmin(lambda, Y, Latent, delta, x, option = 'approx')}), 
     pch=17)
points(tau_list, 
       sapply(tau_list, function (x) {findmin(lambda, Y, Latent, delta, x, option = 'biased')}), 
       col = 3)
points(tau_list, 
       sapply(tau_list, function (x) {findmin(lambda, Y, Latent, delta, x, option = 'truth')}), 
       col = 2)
points(tau_list, 
       sapply(tau_list, function (x) {findmin(lambda, Y, Latent, delta, x, option = 'oracle')}), 
       col = 4)
abline(a = 0, b = 1, lty = 2)
