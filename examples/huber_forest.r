source("examples/metrics.R")
source("examples/help_functions.R")
library(generalizedForest)
library(xgboost)
library(ggplot2)

# Load in the data
n <- 2000
n_test <- 1000
Xtrain <- sort(runif(n = n, min = -5, max = 5))
Xtest <- sort(runif(n = n_test, min = -5, max = 5))
Ytrain <- Xtrain^2 + rnorm(n)
Ytest <- Xtest^2 + rnorm(n_test)
data_train <- cbind.data.frame(Xtrain, Ytrain)
data_test <- cbind.data.frame(Xtest, Ytest)
colnames(data_train) <- c('x', 'y')
colnames(data_test) <- c('x', 'y')

# Add noise (t-distribution)
noise_idx <- sample(c(FALSE, TRUE), size = n, replace = TRUE, prob = c(0.8, 0.2))
data_train$y[noise_idx] <- data_train$y[noise_idx] + 2*rt(n = sum(noise_idx), df = 2)

# plot truth
par(mfrow=c(1,1))
plot(data_train$x, data_train$y)

# build generalizedForest model
grf <- generalizedForest(y ~ ., data = data_train, ntree = 800, nodesize = 5)

# get proximity matrix
proxSup <- getProximityMtx(grf, newdata = data_test)
proxMtx <- proxSup$proximityMtx
dim(proxMtx) # test by train dimension

# fixed point method
iters <- 50
delta <- 0.2
Yrf <- proxSup$predicted
YpH <- runif(n=length(Yrf), min = min(Yrf), max = max(Yrf))
error <- rep(NA, iters)

for (t in 1:iters) {
  Yold <- YpH
  for (i in 1:NROW(data_test)) {
    new_weight <- proxMtx[i, ]/sqrt(1 + ((Yold[i]-grf$data$y)/delta)^2)
    YpH[i] <- t(new_weight)%*%grf$data$y / sum(new_weight)
  }
  error[t] <- sqrt(sum((YpH - Yold)^2))
  gc()
}

plot(error, type = "l")

print(error[iters])

plot(Xtest, YpH, type = 'l', col='red')
lines(Xtest, Yrf, col='blue')
