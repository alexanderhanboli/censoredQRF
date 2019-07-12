setwd("~/Projects/generalizedForest/examples")

source("metrics.R")
source("help_functions.R")
source('crf_km.R')
library(generalizedForest)
library(ggplot2)
library(grf)
library(quantregForest)
library(randomForestSRC)
library(survival)

attach(mtcars)
op <- par(mar=c(4,4,1,1)+0.1, oma = c(0,0,0,0) + 0.1, pty="s")

# Load in the data
n <- 300
n_test <- 300

# training data
Xtrain <- sort(runif(n = n, min = 0, max = 2*pi))
sigma <- 0.3
Ttrain <- sin(Xtrain) + rnorm(n, mean = 0, sd = sigma)
ctrain <- sin(Xtrain) + rexp(n = n, rate = 0.2) - 1.5
Ytrain <- pmin(Ttrain, ctrain)
censorInd <- 1*(Ttrain <= ctrain)
print(mean(censorInd))
data_train <- cbind.data.frame(Xtrain, Ytrain, censorInd)
# plot training data
plot(Xtrain, Ytrain, cex = 0.2)
points(Xtrain[!censorInd], Ytrain[!censorInd], type = 'p', col = 'red', cex = 0.3)
points(Xtrain[!censorInd], Ttrain[!censorInd], type = 'p', col = 'green', cex = 0.3)

# test data
Xtest <- sort(runif(n = n_test, min = 0, max = 2*pi))
Ytest <- sin(Xtest) + rnorm(n_test, mean = 0, sd = sigma)
data_test <- cbind.data.frame(Xtest, Ytest, rep(1, n_test))

colnames(data_train) <- c('x', 'y', 'ind')
colnames(data_test) <- c('x', 'y', 'ind')

# build an AFT model
# aft <- survreg(Surv(y, ind) ~ ., data_train, dist='weibull', scale=1)
# Yaft <- predict(aft, data_test)

# build generalizedForest model
# get quantiles
tau <- 0.7
nodesize <- 5
ntree <- 1000
Yc <- crf.km(y ~ x, ntree = ntree, nodesize = 4*nodesize, data_train = data_train, data_test = data_test, 
             yname = 'y', iname = 'ind', 
             tau = tau)

# # generalized random forest (Stefan's)
# # latent T ~ x
grf_qf_latent <- quantile_forest(data_train[,1,drop=FALSE], Ttrain, quantiles = tau, num.trees = ntree, min.node.size = nodesize)
Ygrf_latent <- predict(grf_qf_latent, Xtest, quantiles = tau)
# # censored Y ~ x
grf_qf <- quantile_forest(data_train[,1,drop=FALSE], Ytrain, quantiles = tau, num.trees = ntree, min.node.size = nodesize)
Ygrf <- predict(grf_qf, Xtest, quantiles = tau)

# # survival forest
# surv_rf <- rfsrc(Surv(y, ind) ~ x, data = data_train, ntree = 1000, nodesize = 100)
# Ysurv <- predict(surv_rf, newdata = data_test)$predicted

# Meinshasen
#qrf_latent <- quantregForest(x=data_train[,1,drop=FALSE], y=Ttrain, nodesize=3*nodesize, ntree=ntree)
#Yqrf_latent <- predict(qrf_latent, data_test[,1,drop=FALSE], what = tau)

# comparison
plot(Xtest, Ytest, cex = 0.04, xlab = 'x', ylab = 'y')
quantiles <- sin(Xtest) + qnorm(tau, 0, sigma)
# lines(Xtest, Yaft)
lines(Xtest, quantiles, col = 'black', cex = 2)
lines(Xtest, Yc$predicted, col='red', lty = 5, cex = 1)
#lines(Xtest, Yc_nodebias$predicted, col='blue', lty = 3, cex = 1)
lines(Xtest, Ygrf, col = 'blue', lty = 2, cex = 1)
#lines(Xtest, Yqrf_latent, col = 'yellow')
lines(Xtest, Ygrf_latent, col = 'black', type = 'b', pch = 18, lty = 1, cex = .5)

# Add a legend
legend(4, 1, legend=c("true quantile", "cRF", "gRF", "gRF-oracle"), 
       lty=c(1, 5, 2, 1), cex=0.8, pch = c(-1,-1,-1, 18), col = c('black', 'red', 'blue', 'black'))

# results
metrics(Ytest, Yc$predicted, quantiles, tau)
metrics(Ytest, Ygrf, quantiles, tau)
metrics(Ytest, Ygrf_latent, quantiles, tau)
# 
# randomForestSRC:::cindex(data_test$y, data_test$ind, Yc)
# 1-randomForestSRC:::cindex(data_test$y, data_test$ind, Ysurv)
# randomForestSRC:::cindex(data_test$y, data_test$ind, Ygrf)
# randomForestSRC:::cindex(data_test$y, data_test$ind, Ygrf_latent)