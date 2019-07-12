source("examples/metrics.R")
source("examples/help_functions.R")
library(generalizedForest)
library(ggplot2)

# Load in the data
n <- 2000
n_test <- 1000

Xtrain <- sort(runif(n = n, min = -5, max = 5))
Ttrain <- Xtrain^2 + rnorm(n)
ctrain <- 2*Xtrain^2*rbinom(n = n, size = 1, prob = 0.5)
Ytrain <- pmax(Ttrain, ctrain)
data_train <- cbind.data.frame(Xtrain, Ytrain)

Xtest <- sort(runif(n = n_test, min = -5, max = 5))
Ytest <- Xtest^2 + rnorm(n_test)
data_test <- cbind.data.frame(Xtest, Ytest)

colnames(data_train) <- c('x', 'y')
colnames(data_test) <- c('x', 'y')

# plot truth
par(mfrow=c(1,1))
plot(data_train$x, data_train$y)

# build generalizedForest model
grf <- generalizedForest(y ~ ., data = data_train, ntree = 800, nodesize = 10)

# get proximity matrix
proxSup <- getProximityMtx(grf, newdata = data_test)
proxMtx <- proxSup$proximityMtx
dim(proxMtx) # test by train dimension

# get quantiles
idx <- order(Ytrain)
orderedMtx <- proxMtx[ ,idx]
orderedCensor <- ctrain[idx]

Yrf <- proxSup$predicted
Yc <- rep(NA, length(Yrf))
Ytrain_sorted <- sort(Ytrain)
tau <- 0.5
for (r in 1:NROW(data_test)) {
  Yc[r] <- Ytrain_sorted[1]
  min_diff <- 10000
  for (lambda in Ytrain_sorted) {
    tmp <- abs(((tau - (Ytrain_sorted < lambda))*orderedMtx[r, ])%*%(orderedCensor < lambda))
    if (tmp > 0 && tmp < min_diff) {
      Yc[r] <- lambda
      min_diff <- tmp
    }
  }
}

par(mfrow=c(1,1))
plot(data_train$x, data_train$y)
lines(Xtest, Yc, col='red')
lines(Xtest, Yrf, col='blue')
