setwd("~/Projects/generalizedForest/examples")

source("metrics.R")
source("help_functions.R")
source("crf_km.R")
library(generalizedForest)
library(quantregForest)
library(ggplot2)
library(randomForestSRC)
library(survival)
library(prodlim)
library(quantreg)
library(grf)
library(mlbench)
library(alr3)

## data, formula specifications
dataset <- 'BostonHousing'
if (dataset == 'engel') {
  data(engel, package = 'quantreg')
  plot(engel, log = "xy",
       main = "'engel' data (log - log scale)")
  plot(log10(foodexp) ~ log10(income), data = engel,
       main = "'engel' data (log10 - transformed)")
  data.na <- na.omit(engel) ##remove NA's
  data.na$income <- log10(data.na$income)
  data.na$foodexp <- log10(data.na$foodexp)
  n <- nrow(data.na)
  c <- rexp(n = n, rate = 0.08)
  data.na$y <- pmin(data.na$foodexp, c)
  data.na$status <- 1*(data.na$foodexp <= c)
  xnam <- 'income'
  ynam <- 'y'
  tnam <- 'foodexp'
  inam <- 'status'
  print(mean(data.na[[inam]]))
} else if (dataset == 'Mammales') {
  data(Mammals)
  attach(Mammals)
  x <- log(weight)
  y <- log(speed)
  plot(x,y, xlab="Weight in log(Kg)", ylab="Speed in log(Km/hour)",type="n")
  points(x[hoppers],y[hoppers],pch = "h", col="red")
  points(x[specials],y[specials],pch = "s", col="blue")
  others <- (!hoppers & !specials)
  points(x[others],y[others], col="black",cex = .75)
  
  data.na <- na.omit(Mammals) ##remove NA's
  data.na$speed <- log(data.na$speed)
  data.na$weight <- log(data.na$weight)
  n <- nrow(data.na)
  c <- rexp(n = n, rate = 0.08)
  data.na$y <- pmin(data.na$speed, c)
  data.na$status <- 1*(data.na$speed <= c)
  xnam <- c('weight', 'hoppers', 'specials')
  ynam <- 'y'
  tnam <- 'speed'
  inam <- 'status'
  print(mean(data.na[[inam]]))
} else if (dataset == 'airquality') {
  data(airquality)
  airquality <- airquality[ !apply(is.na(airquality), 1,any), ]
  
  data.na <- na.omit(airquality) ##remove NA's
  data.na$Ozone <- log(data.na$Ozone)
  n <- nrow(data.na)
  c <- rexp(n = n, rate = 0.08)
  data.na$y <- pmin(data.na$Ozone, c)
  data.na$status <- 1*(data.na$Ozone <= c)
  ynam <- 'y'
  tnam <- 'Ozone'
  inam <- 'status'
  xnam <- colnames(data.na)[-which(colnames(data.na) %in% c(ynam, tnam, inam))]
  print(mean(data.na[[inam]]))
} else if (dataset == 'BostonHousing') {
  data(BostonHousing, package = 'mlbench')
  BostonHousing <- BostonHousing[ !apply(is.na(BostonHousing), 1,any), ]
  data.na <- na.omit(BostonHousing) ##remove NA's
  data.na <- as.data.frame(sapply(data.na, as.numeric))
  n <- nrow(data.na)
  c <- as.double(rexp(n = n, rate = 0.01))
  data.na$y <- pmin(data.na$medv, c)
  data.na$status <- 1*(data.na$medv <= c)
  #data.na$medv = log(data.na$medv)
  #data.na$y = log(data.na$y)
  ynam <- 'y'
  tnam <- 'medv'
  inam <- 'status'
  xnam <- colnames(data.na)[-which(colnames(data.na) %in% c(ynam, tnam, inam))]
  print(mean(data.na[[inam]]))
} else if (dataset == 'Ozone') {
  data(Ozone)
  data.na <- na.omit(Ozone) ##remove NA's
  data.na <- as.data.frame(sapply(data.na, as.numeric))
  n <- nrow(data.na)
  c <- as.double(rexp(n = n, rate = 0.10))
  data.na$y <- pmin(data.na$V4, c)
  data.na$status <- 1*(data.na$V4 <= c)
  data.na$V4 = log(data.na$V4)
  data.na$y = log(data.na$y)
  ynam <- 'y'
  tnam <- 'V4'
  inam <- 'status'
  xnam <- colnames(data.na)[-which(colnames(data.na) %in% c(ynam, tnam, inam))]
  print(mean(data.na[[inam]]))
} else if (dataset == "BigMac") {
  data("BigMac2003")
  data.na <- na.omit(BigMac2003) ##remove NA's
  data.na <- as.data.frame(sapply(data.na, as.numeric))
  n <- nrow(data.na)
  c <- as.double(rexp(n = n, rate = 0.008))
  data.na$y <- pmin(data.na$BigMac, c)
  data.na$status <- 1*(data.na$BigMac <= c)
  ynam <- 'y'
  tnam <- 'BigMac'
  inam <- 'status'
  xnam <- colnames(data.na)[-which(colnames(data.na) %in% c(ynam, tnam, inam))]
  print(mean(data.na[[inam]]))
} else if (dataset == "Abalone") {
  Abalone <- read.csv("~/Downloads/Abalone.csv", row.names=NULL)
  data.na <- na.omit(Abalone) ##remove NA's
  data.na <- as.data.frame(sapply(data.na, as.numeric))
  data.na <- data.na[sample(nrow(data.na), size = 1000, replace = FALSE), ]
  n <- nrow(data.na)
  c <- as.double(rexp(n = n, rate = 0.1))
  data.na$rings <- log(data.na$rings)
  data.na$y <- pmin(data.na$rings, c)
  data.na$status <- 1*(data.na$rings <= c)
  ynam <- 'y'
  tnam <- 'rings'
  inam <- 'status'
  xnam <- colnames(data.na)[-which(colnames(data.na) %in% c(ynam, tnam, inam))]
  print(mean(data.na[[inam]]))
}

one_run <- function(b, data, xnam, ynam, inam, tau, nodesize, ntree) {
  cat("bootstrap step:", b, "\n")
  
  train <- sample(1:nrow(data), nrow(data), replace = TRUE)
  data_train <- data[train, ]
  data_test <- data[-train, ][data[-train, inam] == 1, ]
  print(mean(data_train[[inam]])) 
  
  # build generalizedForest model
  fmla <- as.formula(paste(ynam, "~", paste(xnam, collapse= "+")))
  
  Yc <- crf.km(fmla, ntree = ntree, nodesize = 2*nodesize, data_train = data_train, data_test = data_test, 
               yname = ynam, iname = inam, tau = tau)$predicted
  
  # generalized random forest (Stefan's)
  grf_qf <- quantile_forest(as.matrix(data_train[,xnam,drop=FALSE]), data_train[,ynam], quantiles = c(tau), 
                            num.trees = ntree, min.node.size = nodesize)
  Ygrf <- predict(grf_qf, as.matrix(data_test[,xnam,drop=FALSE]), quantiles = tau)
  
  grf_qf_latent <- quantile_forest(as.matrix(data_train[,xnam,drop=FALSE]), data_train[,tnam], quantiles = tau, 
                            num.trees = ntree, min.node.size = nodesize)
  Ygrf_latent <- predict(grf_qf_latent, as.matrix(data_test[,xnam,drop=FALSE]), quantiles = tau)
  
  # quantile random forest (Meinshasen)
  qrf <- quantregForest(x=data_train[,xnam,drop=FALSE], y=data_train[,ynam], nodesize=2*nodesize, ntree=ntree)
  Yqrf <- predict(qrf, data_test[,xnam,drop=FALSE], what = tau)
  
  qrf_latent <- quantregForest(x=data_train[,xnam,drop=FALSE], y=data_train[,tnam], nodesize=2*nodesize, ntree=ntree)
  Yqrf_latent <- predict(qrf_latent, data_test[,xnam,drop=FALSE], what = tau)
  
  # survival forest
  surv_rf <- rfsrc(Surv(y, status) ~ ., data = data_train, ntree = ntree, nodesize = nodesize)
  Ysurv <- predict(object = surv_rf, newdata = data_test)$predicted
  
  # results
  return(
    list(
      'crf'=metrics(data_test[[tnam]], Yc, data_test[[tnam]], tau),
      'qrf'=metrics(data_test[[tnam]], Yqrf, data_test[[tnam]], tau),
      'grf'=metrics(data_test[[tnam]], Ygrf, data_test[[tnam]], tau),
      
      'qrf_latent'=metrics(data_test[[tnam]], Yqrf_latent, data_test[[tnam]], tau),
      'grf_latent'=metrics(data_test[[tnam]], Ygrf_latent, data_test[[tnam]], tau),
      
      'crf_c'=randomForestSRC:::cindex(data_test[[tnam]], data_test[[inam]], Yc),
      'qrf_c'=randomForestSRC:::cindex(data_test[[tnam]], data_test[[inam]], Yqrf),
      'rsf_c'=1-randomForestSRC:::cindex(data_test[[tnam]], data_test[[inam]], Ysurv),
      'grf_c'=randomForestSRC:::cindex(data_test[[tnam]], data_test[[inam]], Ygrf),
      
      'qrf_c_latent'=randomForestSRC:::cindex(data_test[[tnam]], data_test[[inam]], Yqrf_latent),
      'grf_c_latent'=randomForestSRC:::cindex(data_test[[tnam]], data_test[[inam]], Ygrf_latent)
    )
  )
}

B <- 40
nodesize <- 10
ntree <- 1000
tau <- 0.9
quantile_result <- list('crf'=rep(NA,B), 'qrf'=rep(NA,B), 'grf'=rep(NA,B), 
                        'qrf_oracle'=rep(NA,B), 'grf_oracle'=rep(NA,B))
cindex_result <- list('crf'=rep(NA,B), 'qrf'=rep(NA,B), 'rsf'=rep(NA,B), 'grf'=rep(NA,B),
                      'qrf_oracle'=rep(NA,B), 'grf_oracle'=rep(NA,B))
for (t in 1:B) {
  tmp <- one_run(t, data.na, xnam, ynam, inam, tau, nodesize, ntree)
  
  quantile_result$crf[t] <- tmp$crf['quantile_loss']
  quantile_result$qrf[t] <- tmp$qrf['quantile_loss']
  quantile_result$grf[t] <- tmp$grf['quantile_loss']
  
  quantile_result$qrf_oracle[t] <- tmp$qrf_latent['quantile_loss']
  quantile_result$grf_oracle[t] <- tmp$grf_latent['quantile_loss']
  
  cindex_result$crf[t] <- tmp$crf_c
  cindex_result$qrf[t] <- tmp$qrf_c
  cindex_result$rsf[t] <- tmp$rsf_c
  cindex_result$grf[t] <- tmp$grf_c
  
  cindex_result$qrf_oracle[t] <- tmp$qrf_c_latent
  cindex_result$grf_oracle[t] <- tmp$grf_c_latent
} 

# boxplot
require(reshape2)

dd <- as.data.frame(cindex_result)
ggplot(data = melt(dd), aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=variable)) + labs(x = "Model", y = paste("C-index, tau =", tau), fill = "Model")
ggsave(paste0(dataset, "_cindex_result_", 10*tau, ".pdf"), width = 5, height = 5, path = "./results/")

dd <- as.data.frame(quantile_result)
ggplot(data = melt(dd), aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=variable)) + labs(x = "Model", y = paste("Quantile loss, tau =", tau), fill = "Model")
ggsave(paste0(dataset, "_quantile_loss_result_", 10*tau, ".pdf"), width = 5, height = 5, path = "./results/")