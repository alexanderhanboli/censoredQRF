n <- 2000
Latent <- runif(n)
b <- rbinom(n, 1, 0.5)
# C <- b*rnorm(n, mean = 0.1, sd = .5) + (1-b)*rnorm(n, mean = 1.9, sd = .5)
C <- 1*rbinom(n, 1, 0.5) + rnorm(n)
# C <- 0.8
Y <- pmin(Latent, C)
delta <- 1*(Latent <= C)

hist(C, breaks = 30)

obs_C <- Y[delta == 0]
hist(obs_C, breaks = 30)
