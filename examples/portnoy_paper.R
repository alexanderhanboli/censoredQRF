setwd("~/Projects/generalizedForest/examples")

source("metrics.R")
source("help_functions.R")
library(generalizedForest)
library(ggplot2)
library(grf)
library(randomForestSRC)
library(survival)

# Construct data
gen_data <- function(model = 1, censor = 1, n = 500, taus = c(0.3, 0.5, 0.7)) {
  X <- sort(runif(n, 0, 2))
  U <- rnorm(n)
  V <- rnorm(n)
  beta0 <- 5
  beta1 <- 1
  sigma0 <- 0.39
  sigma1 <- 0.09
  sigma2 <- 0.30
  # Latent variables
  if (model == 1) {
    yLatent <- beta0 + beta1 * X + sigma0 * U
    quantiles <- sapply(taus, function(x) {sigma0*qnorm(x) + beta0 + beta1 * X})
  } else {
    yLatent <- beta0 + beta1 * X + (sigma1 + sigma2 * X^2) * U
    quantiles <- sapply(taus, function(x) {(sigma1 + sigma2 * X^2)*qnorm(x) + beta0 + beta1 * X})
  }
  # Censoring variables
  if (censor == 1) {
    c <- 6.5
  } else {
    c <- 5.5 + 0.75 * X + sigma2 * V
  }
  # Response variables
  y <- pmin(yLatent, c)
  status <- 1*(yLatent <= c)
  return (
    list(
      'data' = cbind.data.frame(X, y, status),
      'oracle_data' = cbind.data.frame(X, yLatent),
      'quantiles' = quantiles,
      'censor' = c
      )
  )
}

taus <- c(0.4, 0.5, 0.6)
obj.train <- gen_data(model = 2, censor = 2, n = 2000, taus)
obj.test <- gen_data(model = 2, censor = 2, n = 500, taus)
data.train <- obj.train$data
quantiles.train <- obj.train$quantiles
censor.train <- obj.train$censor
data.test <- obj.test$data
quantiles.test <- obj.test$quantiles

# plot data
plot(data.test$X[data.test$status==1], data.test$y[data.test$status==1], cex = 0.2, 
     ylim = c(4.5, 7.0), xlab = 'X', ylab = 'y')
points(data.test$X[data.test$status==0], data.test$y[data.test$status==0], type = 'p', col = 'green', cex = 0.3)
lines(data.test$X, quantiles.test[,1])
lines(data.test$X, quantiles.test[,2])
lines(data.test$X, quantiles.test[,3])

# build generalizedForest model
source('crf_function.R')
source('crf_km.R')
tau <- 0.5
nodesize <- 200
Yc <- crf.km(y ~ X, ntree = 2000, nodesize = nodesize, data_train = data.train, data_test = data.test, 
          yname = 'y', iname = 'status', tau = tau)$predicted

lines(data.test$X, Yc, col = 'red')

# Yc_bias <- crf(y ~ X, ntree = 2000, nodesize = nodesize, data_train = data.train, data_test = data.test, 
#           yname = 'y', iname = 'status', tau = tau, fixed_censoring = FALSE, debias = FALSE)$predicted
# 
# lines(data.test$X, Yc_bias, col = 'blue')

#oracle model
#generalized random forest (Stefan's)
grf.latent <- quantile_forest(obj.train$oracle_data[,1,drop=FALSE],
                              obj.train$oracle_data$y,
                              quantiles = tau, num.trees = 2000, min.node.size = nodesize)
Ygrf.latent <- predict(grf.latent, data.test$X, quantiles = tau)
lines(data.test$X, Ygrf.latent, col = 'cyan')
# 
# grf <- quantile_forest(data.train[,1,drop=FALSE], 
#                               data.train$y, 
#                               quantiles = 0.6, num.trees = 2000, min.node.size = 200)
# Ygrf <- predict(grf, data.test$X, quantiles = 0.6)
# lines(data.test$X, Ygrf, col = 6)
