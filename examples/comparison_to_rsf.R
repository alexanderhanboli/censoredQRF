setwd("~/Projects/generalizedForest/examples")
source("metrics.R")
source("help_functions.R")
source("crf_km.R")
library(generalizedForest)
library(quantregForest)
library(ggplot2)
library(grf)
library(randomForestSRC)

data("veteran", package = "randomForestSRC")
v.rsf <- rfsrc(Surv(time, status) ~ ., data = veteran[c(1:100),], ntree = 100)
plot.survival(v.rsf)
pred <- predict(v.rsf, newdata = veteran[c(101:137),])
pred$survival
dim(pred$survival)

find_median <- function(surv, max_value){
  bin_len <- dim(surv)[2]
  b_size <- dim(surv)[1]
  bin <- seq(0, max_value, length.out = bin_len)
  quantiles <- rep(0, b_size)
  for (i in 1:b_size){
    check <- surv[i,] > 0.5
    quantile <- bin[sum(check)]
    quantiles[i] <- quantile
  }
  
  return(quantiles)
}