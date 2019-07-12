source("metrics.R")
source("help_functions.R")
library(generalizedForest)
library(ggplot2)
library(grf)
library(randomForestSRC)

# Load in the data
n <- 2000
n_test <- 1000
p <- 10

# training data
Xtrain <- matrix(rnorm(n = n*p), nrow = n, ncol = p)
beta <- c(0.3,0.2,0.5,0.2,0.7,1.0,0,0,0,0)
sigma <- 1.0
Ttrain <- exp(Xtrain%*%beta + rnorm(n, mean = 0, sd = sigma))
ctrain <- rexp(n = n, rate = 0.08)
Ytrain <- pmin(Ttrain, ctrain)
censorInd <- (Ttrain <= ctrain)
data_train <- cbind.data.frame(Xtrain, Ytrain, censorInd)
# plot training data
plot(Xtrain[,1], Ytrain, cex = 0.2)
points(Xtrain[!censorInd,1], Ytrain[!censorInd], type = 'p', col = 'red', cex = 0.3)
points(Xtrain[!censorInd,1], Ttrain[!censorInd], type = 'p', col = 'green', cex = 0.3)

# test data
Xtest <- matrix(rnorm(n = n_test*p), nrow = n_test, ncol = p)
Ytest <- exp(Xtest%*%beta + rnorm(n_test, mean = 0, sd = sigma))
data_test <- cbind.data.frame(Xtest, Ytest, rep(TRUE, n_test))

xnam <- paste0('x', 1:p)
colnames(data_train) <- c(xnam, 'y', 'ind')
colnames(data_test) <- c(xnam, 'y', 'ind')

# build generalizedForest model
fmla <- as.formula(paste("y ~ ", paste(xnam, collapse= "+")))
grf <- generalizedForest(fmla, 
                         data = data_train, ntree = 1000, nodesize = 100)

# get proximity matrix
proxSup <- getProximityMtx(grf, newdata = data_test)
proxMtx <- proxSup$proximityMtx
dim(proxMtx) # test by train dimension

# get quantiles
tau <- 0.5

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

# quantile forest
# Yqf <- rep(NA, length(Yrf))
# orderedProx <- proxMtx[ ,order(Ytrain)]
# for (r in 1:NROW(data_test)) {
#   rhs <- 0
#   j <- 1
#   for (lambda in Ytrain_sorted) {
#     rhs <- rhs + orderedProx[r, j]
#     j <- j + 1
#     if (rhs >= tau) {
#       break
#     }
#   }
#   Yqf[r] <- lambda
# }

# generalized random forest (Stefan's)
grf_qf_latent <- quantile_forest(data_train[,1:p,drop=FALSE], Ttrain, quantiles = tau, num.trees = 1000, min.node.size = 100)
Ygrf_latent <- predict(grf_qf_latent, Xtest, quantiles = tau)

grf_qf <- quantile_forest(data_train[,1:p,drop=FALSE], Ytrain, quantiles = tau, num.trees = 1000, min.node.size = 100)
Ygrf <- predict(grf_qf, Xtest, quantiles = tau)

# survival forest
# surv_rf <- rfsrc(Surv(y, ind) ~ ., data = data_train, ntree = 1000, nodesize = 100, splitrule = 'logrank')
# Ysurv <- predict(object = surv_rf, newdata = data_test)

# results
metrics(Ytest, Yc)
metrics(Ytest, Ygrf)
metrics(Ytest, Ygrf_latent)
