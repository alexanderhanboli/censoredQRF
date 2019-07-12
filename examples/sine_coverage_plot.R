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
nodesize <- 10
ntree <- 1000
Yc_l <- crf.km(y ~ x, ntree = ntree, nodesize = 3*nodesize, data_train = data_train, data_test = data_test, 
             yname = 'y', iname = 'ind', 
             tau = 0.025)$predicted
Yc_u <- crf.km(y ~ x, ntree = ntree, nodesize = 3*nodesize, data_train = data_train, data_test = data_test, 
               yname = 'y', iname = 'ind', 
               tau = 0.975)$predicted

# # generalized random forest (Stefan's)
# # latent T ~ x
grf_qf_latent <- quantile_forest(data_train[,1,drop=FALSE], Ttrain, 
                                 quantiles = c(0.025,0.1,0.5,0.9,0.975), 
                                 num.trees = ntree, min.node.size = nodesize)
Ygrf_latent_l <- predict(grf_qf_latent, Xtest, quantiles = 0.025)
Ygrf_latent_u <- predict(grf_qf_latent, Xtest, quantiles = 0.975)
# # censored Y ~ x
grf_qf <- quantile_forest(data_train[,1,drop=FALSE], Ytrain, 
                          quantiles = c(0.025,0.1,0.5,0.9,0.975), 
                          num.trees = ntree, min.node.size = nodesize)
Ygrf_l <- predict(grf_qf, Xtest, quantiles = 0.025)
Ygrf_u <- predict(grf_qf, Xtest, quantiles = 0.975)

# Meinshasen
qrf_latent <- quantregForest(x=data_train[,1,drop=FALSE], y=Ttrain, nodesize=3*nodesize, ntree=ntree)
Yqrf_latent_l <- predict(qrf_latent, data_test[,1,drop=FALSE], what = 0.025)
Yqrf_latent_u <- predict(qrf_latent, data_test[,1,drop=FALSE], what = 0.975)

qrf <- quantregForest(x=data_train[,1,drop=FALSE], y=Ytrain, nodesize=3*nodesize, ntree=ntree)
Yqrf_l <- predict(qrf, data_test[,1,drop=FALSE], what = 0.025)
Yqrf_u <- predict(qrf, data_test[,1,drop=FALSE], what = 0.975)

# comparison
data_test$crf_l <- Yc_l
data_test$crf_u <- Yc_u
data_test$grf_l <- Ygrf_l
data_test$grf_u <- Ygrf_u
data_test$grf_oracle_l <- Ygrf_latent_l
data_test$grf_oracle_u <- Ygrf_latent_u
data_test$qrf_l <- Yqrf_l
data_test$qrf_u <- Yqrf_u
data_test$qrf_oracle_l <- Yqrf_latent_l
data_test$qrf_oracle_u <- Yqrf_latent_u
data_test$ci_l <- sin(Xtest) + qnorm(0.025, 0, sigma)
data_test$ci_u <- sin(Xtest) + qnorm(0.975, 0, sigma)

ggplot(data = data_test, aes(x=x, y=y)) +
  geom_point() +
  geom_line(aes(y = crf_l)) + 
  geom_line(aes(y = crf_u)) + 
  geom_line(aes(y = ci_l), colour='red') + 
  geom_line(aes(y = ci_u), colour='red') + 
  geom_ribbon(aes(ymin = crf_l, ymax = crf_u), fill = "grey", alpha = .2)
ggsave("crf_sine_ci_plot.pdf", width = 5, height = 5, path = "./results/")

ggplot(data = data_test, aes(x=x, y=y)) +
  geom_point() +
  geom_line(aes(y = grf_l)) + 
  geom_line(aes(y = grf_u)) + 
  geom_line(aes(y = ci_l), colour='red') + 
  geom_line(aes(y = ci_u), colour='red') + 
  geom_ribbon(aes(ymin = crf_l, ymax = crf_u), fill = "grey", alpha = .2)
ggsave("grf_sine_ci_plot.pdf", width = 5, height = 5, path = "./results/")

ggplot(data = data_test, aes(x=x, y=y)) +
  geom_point() +
  geom_line(aes(y = grf_oracle_l)) + 
  geom_line(aes(y = grf_oracle_u)) + 
  geom_line(aes(y = ci_l), colour='red') + 
  geom_line(aes(y = ci_u), colour='red') + 
  geom_ribbon(aes(ymin = crf_l, ymax = crf_u), fill = "grey", alpha = .2)
ggsave("grf_oracle_sine_ci_plot.pdf", width = 5, height = 5, path = "./results/")

ggplot(data = data_test, aes(x=x, y=y)) +
  geom_point() +
  geom_line(aes(y = qrf_l)) + 
  geom_line(aes(y = qrf_u)) + 
  geom_line(aes(y = ci_l), colour='red') + 
  geom_line(aes(y = ci_u), colour='red') + 
  geom_ribbon(aes(ymin = crf_l, ymax = crf_u), fill = "grey", alpha = .2)
ggsave("qrf_sine_ci_plot.pdf", width = 5, height = 5, path = "./results/")

ggplot(data = data_test, aes(x=x, y=y)) +
  geom_point() +
  geom_line(aes(y = qrf_oracle_l)) + 
  geom_line(aes(y = qrf_oracle_u)) + 
  geom_line(aes(y = ci_l), colour='red') + 
  geom_line(aes(y = ci_u), colour='red') + 
  geom_ribbon(aes(ymin = crf_l, ymax = crf_u), fill = "grey", alpha = .2)
ggsave("qrf_oracle_sine_ci_plot.pdf", width = 5, height = 5, path = "./results/")