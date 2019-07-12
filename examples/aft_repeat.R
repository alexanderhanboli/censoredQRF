setwd("~/Projects/generalizedForest/examples")

source("metrics.R")
source("help_functions.R")
source("crf_km.R")
library(generalizedForest)
library(quantregForest)
library(ggplot2)
library(grf)
library(randomForestSRC)
library(survival)

# Load in the data
n <- 300
n_test <- 300

one_run <- function(n, n_test, tau) {
  # training data
  Xtrain <- sort(runif(n = n, min = 0, max = 2))
  sigma <- 0.3
  Ttrain <- exp(Xtrain + rnorm(n, mean = 0, sd = sigma))
  ctrain <- rexp(n = n, rate = 0.08)
  Ytrain <- pmin(Ttrain, ctrain)
  censorInd <- 1*(Ttrain <= ctrain)
  print(mean(censorInd))
  data_train <- cbind.data.frame(Xtrain, Ytrain, censorInd)
  # test data
  Xtest <- sort(runif(n = n_test, min = 0, max = 2))
  quantile_test <- exp(Xtest + qnorm(tau, 0, sigma))
  Ytest <- exp(Xtest + rnorm(n_test, mean = 0, sd = sigma))
  data_test <- cbind.data.frame(Xtest, Ytest, rep(1, n_test))
  # column names
  colnames(data_train) <- c('x', 'y', 'ind')
  colnames(data_test) <- c('x', 'y', 'ind')
  
  nodesize <- 5
  ntree <- 1000
  Yc <- crf.km(y ~ x, ntree = ntree, nodesize = 8*nodesize, data_train = data_train, data_test = data_test, 
               yname = 'y', iname = 'ind', 
               tau = tau)
  
  # generalized random forest (Stefan's)
  # latent T ~ x
  grf_qf_latent <- quantile_forest(data_train[,1,drop=FALSE], Ttrain, quantiles = tau, num.trees = ntree, min.node.size = nodesize)
  Ygrf_latent <- predict(grf_qf_latent, Xtest, quantiles = tau)
  # censored Y ~ x
  grf_qf <- quantile_forest(data_train[,1,drop=FALSE], Ytrain, quantiles = tau, num.trees = ntree, min.node.size = nodesize)
  Ygrf <- predict(grf_qf, Xtest, quantiles = tau)
  
  # quantile random forest (Meinshasen)
  qrf_latent <- quantregForest(x=data_train[,1,drop=FALSE], y=Ttrain, nodesize=8*nodesize, ntree=ntree)
  Yqrf_latent <- predict(qrf_latent, data_test[,1,drop=FALSE], what = tau)
  
  qrf <- quantregForest(x=data_train[,1,drop=FALSE], y=Ytrain, nodesize=8*nodesize, ntree=ntree)
  Yqrf <- predict(qrf, data_test[,1,drop=FALSE], what = tau)
  
  # survival forest
  surv_rf <- rfsrc(Surv(y, ind) ~ x, data = data_train, ntree = ntree, nodesize = nodesize)
  Ysurv <- predict(surv_rf, newdata = data_test)$predicted
  
  # results
  return(
    list(
      'crf'=metrics(data_test$y, Yc$predicted, quantile_test, tau),
      'qrf'=metrics(data_test$y, Yqrf, quantile_test, tau),
      'qrf_latent'=metrics(data_test$y, Yqrf_latent, quantile_test, tau),
      'grf'=metrics(data_test$y, Ygrf, quantile_test, tau),
      'grf_latent'=metrics(data_test$y, Ygrf_latent, quantile_test, tau),
      
      'crf_c'=randomForestSRC:::cindex(data_test$y, data_test$ind, Yc$predicted),
      'qrf_c'=randomForestSRC:::cindex(data_test$y, data_test$ind, Yqrf),
      'qrf_latent_c'=randomForestSRC:::cindex(data_test$y, data_test$ind, Yqrf_latent),
      'rsf_c'=1-randomForestSRC:::cindex(data_test$y, data_test$ind, Ysurv),
      'grf_c'=randomForestSRC:::cindex(data_test$y, data_test$ind, Ygrf),
      'grf_latent_c'=randomForestSRC:::cindex(data_test$y, data_test$ind, Ygrf_latent)
    )
  )
}

B = 40
tau <- 0.1
mse_result <- list('crf'=rep(NA,B), 'qrf'=rep(NA,B), 'grf'=rep(NA,B), 'qrf_oracle'=rep(NA,B), 'grf_oracle'=rep(NA,B))
mad_result <- list('crf'=rep(NA,B), 'qrf'=rep(NA,B), 'grf'=rep(NA,B), 'qrf_oracle'=rep(NA,B), 'grf_oracle'=rep(NA,B))
quantile_result <- list('crf'=rep(NA,B), 'qrf'=rep(NA,B), 'grf'=rep(NA,B), 'qrf_oracle'=rep(NA,B), 'grf_oracle'=rep(NA,B))
cindex_result <- list('crf'=rep(NA,B), 'rsf'=rep(NA,B), 'qrf'=rep(NA,B), 'grf'=rep(NA,B), 'qrf_oracle'=rep(NA,B), 'grf_oracle'=rep(NA,B))
for (t in 1:B) {
  print(t)
  
  tmp <- one_run(n, n_test, tau)
  
  mse_result$crf[t] <- tmp$crf['mse']
  mse_result$qrf[t] <- tmp$qrf['mse']
  mse_result$qrf_oracle[t] <- tmp$qrf_latent['mse']
  mse_result$grf[t] <- tmp$grf['mse']
  mse_result$grf_oracle[t] <- tmp$grf_latent['mse']
  
  mad_result$crf[t] <- tmp$crf['mad']
  mad_result$qrf[t] <- tmp$qrf['mad']
  mad_result$qrf_oracle[t] <- tmp$qrf_latent['mad']
  mad_result$grf[t] <- tmp$grf['mad']
  mad_result$grf_oracle[t] <- tmp$grf_latent['mad']
  
  quantile_result$crf[t] <- tmp$crf['quantile_loss']
  quantile_result$qrf[t] <- tmp$qrf['quantile_loss']
  quantile_result$qrf_oracle[t] <- tmp$qrf_latent['quantile_loss']
  quantile_result$grf[t] <- tmp$grf['quantile_loss']
  quantile_result$grf_oracle[t] <- tmp$grf_latent['quantile_loss']
  
  cindex_result$crf[t] <- tmp$crf_c
  cindex_result$rsf[t] <- tmp$rsf_c
  cindex_result$qrf[t] <- tmp$qrf_c
  cindex_result$qrf_oracle[t] <- tmp$qrf_latent_c
  cindex_result$grf[t] <- tmp$grf_c
  cindex_result$grf_oracle[t] <- tmp$grf_latent_c
}

#sapply(mse_result, mean)
#sapply(mse_result, sd)

#sapply(quantile_result, mean)
#sapply(quantile_result, sd)

#sapply(cindex_result, mean)
#sapply(cindex_result, sd)

# boxplot
require(reshape2)
dd <- as.data.frame(quantile_result)
ggplot(data = melt(dd), aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=variable)) + labs(x = "Model", y = paste("Quantile loss, tau =", tau), fill = "Model")
ggsave(paste0("quantile_loss_result_", 10*tau, ".pdf"), width = 5, height = 5, path = "./results/")

dd <- as.data.frame(mse_result)
ggplot(data = melt(dd), aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=variable)) + labs(x = "Model", y = paste("MSE, tau =", tau), fill = "Model")
ggsave(paste0("mse_result_", 10*tau, ".pdf"), width = 5, height = 5, path = "./results/")

dd <- as.data.frame(mad_result)
ggplot(data = melt(dd), aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=variable)) + labs(x = "Model", y = paste("MAD, tau =", tau), fill = "Model")
ggsave(paste0("mad_result_", 10*tau, ".pdf"), width = 5, height = 5, path = "./results/")

dd <- as.data.frame(cindex_result)
ggplot(data = melt(dd), aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=variable)) + labs(x = "Model", y = paste("C-index, tau =", tau), fill = "Model")
ggsave(paste0("cindex_result_", 10*tau, ".pdf"), width = 5, height = 5, path = "./results/")
