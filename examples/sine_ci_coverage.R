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
n <- 1000
n_test <- 300

one_run <- function(n, n_test, nodesize) {
  # training data
  Xtrain <- sort(runif(n = n, min = 0, max = 2*pi))
  sigma <- 0.3
  Ttrain <- sin(Xtrain) + rnorm(n, mean = 0, sd = sigma) + 2.5
  ctrain <- sin(Xtrain) + rexp(n = n, rate = 0.2) + 1
  Ytrain <- pmin(Ttrain, ctrain)
  censorInd <- 1*(Ttrain <= ctrain)
  print(mean(censorInd))
  data_train <- cbind.data.frame(Xtrain, Ytrain, censorInd)
  # test data
  Xtest <- sort(runif(n = n_test, min = 0, max = 2*pi))
  Ytest <- sin(Xtest) + rnorm(n_test, mean = 0, sd = sigma) + 2.5
  data_test <- cbind.data.frame(Xtest, Ytest, rep(1, n_test))
  # column names
  colnames(data_train) <- c('x', 'y', 'ind')
  colnames(data_test) <- c('x', 'y', 'ind')
  
  ntree <- 1000
  Yc_l <- crf.km(y ~ x, ntree = ntree, nodesize = nodesize, data_train = data_train, data_test = data_test, 
               yname = 'y', iname = 'ind', 
               tau = 0.025)$predicted
  Yc_u <- crf.km(y ~ x, ntree = ntree, nodesize = nodesize, data_train = data_train, data_test = data_test, 
                 yname = 'y', iname = 'ind', 
                 tau = 0.975)$predicted
  Yc_coverage <- mean((Yc_l <= Ytest) & (Ytest <= Yc_u))
  
  # generalized random forest (Stefan's)
  # latent T ~ x
  grf_qf_latent <- quantile_forest(data_train[,1,drop=FALSE], Ttrain, 
                                   quantiles = c(0.025,0.1,0.5,0.9,0.975), 
                                   num.trees = ntree, min.node.size = nodesize)
  Ygrf_latent_l <- predict(grf_qf_latent, Xtest, quantiles = 0.025)
  Ygrf_latent_u <- predict(grf_qf_latent, Xtest, quantiles = 0.975)
  Ygrf_latent_coverage <- mean((Ygrf_latent_l <= Ytest) & (Ytest <= Ygrf_latent_u))
  
  # censored Y ~ x
  grf_qf <- quantile_forest(data_train[,1,drop=FALSE], Ytrain, 
                            quantiles = c(0.025,0.1,0.5,0.9,0.975), 
                            num.trees = ntree, min.node.size = nodesize)
  Ygrf_l <- predict(grf_qf, Xtest, quantiles = 0.025)
  Ygrf_u <- predict(grf_qf, Xtest, quantiles = 0.975)
  Ygrf_coverage <- mean((Ygrf_l <= Ytest) & (Ytest <= Ygrf_u))
  
  # quantile random forest (Meinshasen)
  qrf_latent <- quantregForest(x=data_train[,1,drop=FALSE], y=Ttrain, nodesize=nodesize, ntree=ntree)
  Yqrf_latent_l <- predict(qrf_latent, data_test[,1,drop=FALSE], what = 0.025)
  Yqrf_latent_u <- predict(qrf_latent, data_test[,1,drop=FALSE], what = 0.975)
  Yqrf_latent_coverage <- mean((Yqrf_latent_l <= Ytest) & (Ytest <= Yqrf_latent_u))
  
  qrf <- quantregForest(x=data_train[,1,drop=FALSE], y=Ytrain, nodesize=nodesize, ntree=ntree)
  Yqrf_l <- predict(qrf, data_test[,1,drop=FALSE], what = 0.025)
  Yqrf_u <- predict(qrf, data_test[,1,drop=FALSE], what = 0.975)
  Yqrf_coverage <- mean((Yqrf_l <= Ytest) & (Ytest <= Yqrf_u))
  
  # results
  return(
    list(
      'crf'=Yc_coverage,
      'qrf'=Yqrf_coverage,
      'grf'=Ygrf_coverage,
      'qrf_latent'=Yqrf_latent_coverage,
      'grf_latent'=Ygrf_latent_coverage
    )
  )
}

B <- 300
nodesize <- 0
ntree <- 1000
result <- list('nodesize'=rep(NA,B), 'crf'=rep(NA,B), 'qrf'=rep(NA,B), 'grf'=rep(NA,B), 
               'qrf_oracle'=rep(NA,B), 'grf_oracle'=rep(NA,B))
for (t in 1:B) {
  print(t)
  
  if (t%%20 == 1) {
    nodesize <- nodesize + 5
  }
  
  tmp <- one_run(n, n_test, nodesize)
  
  result$nodesize[t] <- nodesize
  result$crf[t] <- tmp$crf
  result$qrf[t] <- tmp$qrf
  result$grf[t] <- tmp$grf
  result$qrf_oracle[t] <- tmp$qrf_latent
  result$grf_oracle[t] <- tmp$grf_latent
}

# Plot
require(reshape2)
dd <- melt(as.data.frame(result), id.vars = 'nodesize')
dd.agg <- aggregate(value ~ nodesize + variable, dd, function(x) c(mean = mean(x), sd = sd(x)))
dd.agg$mean <- dd.agg[-1][[2]][,1]
dd.agg$sd <- dd.agg[-1][[2]][,2]
dd.agg$value <- NULL
ggplot(data = dd.agg, aes(x=nodesize, y=mean, colour=variable)) + 
  # geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1) +
  geom_line() +
  geom_point() +
  labs(x = "Nodesize", y = "Coverage", fill = "Nodesize")
ggsave(paste0("sine_coverage_nodesize.pdf"), width = 5, height = 5, path = "./results/")
