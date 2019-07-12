source("examples/metrics.R")
source("examples/help_functions.R")
library(generalizedForest)
library(grf)

# Load in the data
n <- 2000
n_test <- 2000

# training data
Xtrain <- sort(runif(n = n, min = 0, max = 2))
sigma <- 0.4
Ttrain <- exp(Xtrain + rnorm(n, mean = 0, sd = sigma))
ctrain <- 30*(rbinom(n = n, size = 1, prob = 0.9) - 0.2) # missing value
Ytrain <- pmin(Ttrain, ctrain)
censorInd <- (Ttrain < ctrain)
data_train <- cbind.data.frame(Xtrain, Ytrain, censorInd)
# plot training data
plot(Xtrain, Ytrain, cex = 0.2)
points(Xtrain[!censorInd], Ytrain[!censorInd], type = 'p', col = 'red', cex = 0.3)
points(Xtrain[!censorInd], Ttrain[!censorInd], type = 'p', col = 'green', cex = 0.3)

# test data
Xtest <- sort(runif(n = n_test, min = 0, max = 2))
Ytest <- exp(Xtest + rnorm(n_test, mean = 0, sd = sigma))
data_test <- cbind.data.frame(Xtest, Ytest, rep(TRUE, n_test))

colnames(data_train) <- c('x', 'y', 'ind')
colnames(data_test) <- c('x', 'y', 'ind')

# build generalizedForest model
grf <- generalizedForest(y ~ x, data = data_train, ntree = 1000, nodesize = 200)

# get proximity matrix
proxSup <- getProximityMtx(grf, newdata = data_test)
proxMtx <- proxSup$proximityMtx
dim(proxMtx) # test by train dimension

# get quantiles
tau <- 0.2

# censor forest
Yrf <- proxSup$predicted # RF prediction
Yc <- rep(NA, length(Yrf))
Ytrain_sorted <- sort(Ytrain)
candidateY <- seq(min(Ytrain), max(Ytrain), length.out = 2000)
right_censoring <- TRUE

for (r in 1:NROW(data_test)) {
  Yc[r] <- candidateY[1]
  min_loss <- 10000
  for (lambda in candidateY) {
    if (right_censoring) {
      # oracle
      # loss <- proxMtx[r, ]%*%quantile_loss(Ytrain - pmin(ctrain, lambda), tau)
      # practice
      loss <- proxMtx[r, ]%*%quantile_loss(Ytrain - (1-censorInd)*pmin(Ytrain, lambda) - censorInd*lambda, tau)
    } else {
      #
    }
    if (loss < min_loss) {
      Yc[r] <- lambda
      min_loss <- loss
    }
  }
}

# generalized random forest (Stefan's)
grf_qf <- quantile_forest(data_train[,1,drop=FALSE], Ytrain, quantiles = tau, num.trees = 1000, min.node.size = 200)
Ygrf <- predict(grf_qf, Xtest, quantiles = tau)

# comparison
plot(Xtest, Ytest, cex = 0.1)
quantiles <- exp(Xtest + qnorm(tau, 0, sigma))
lines(Xtest, quantiles, col = 'red', cex = 10) # true tau quantile
lines(Xtest, Yc, col='blue') # my method
lines(Xtest, Ygrf, col = 'green') # Stefan's forest