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
dataset <- 'Abalone'
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
  c <- as.double(rexp(n = n, rate = 0.02))
  data.na$y <- pmin(data.na$V4, c)
  data.na$status <- 1*(data.na$V4 <= c)
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
  data.na <- data.na[sample(nrow(data.na), size = 500, replace = FALSE), ]
  n <- nrow(data.na)
  data.na$rings <- log(data.na$rings)
  c <- as.double(rexp(n = n, rate = 0.10))
  data.na$y <- pmin(data.na$rings, c)
  data.na$status <- 1*(data.na$rings <= c)
  ynam <- 'y'
  tnam <- 'rings'
  inam <- 'status'
  xnam <- colnames(data.na)[-which(colnames(data.na) %in% c(ynam, tnam, inam))]
  print(mean(data.na[[inam]]))
}

one_run <- function(b, data, xnam, ynam, inam, nodesize, ntree) {
  cat("bootstrap step:", b, "\n")
  
  train <- sample(1:nrow(data), nrow(data), replace = TRUE)
  data_train <- data[train, ]
  data_test <- data[-train, ][data[-train, inam] == 1, ]
  print(mean(data_train[[inam]])) 
  
  # build generalizedForest model
  fmla <- as.formula(paste(ynam, "~", paste(xnam, collapse= "+")))
  
  Yc_l <- crf.km(fmla, ntree = ntree, nodesize = nodesize, data_train = data_train, data_test = data_test, 
               yname = ynam, iname = inam, tau = 0.025)$predicted
  Yc_u <- crf.km(fmla, ntree = ntree, nodesize = nodesize, data_train = data_train, data_test = data_test, 
                 yname = ynam, iname = inam, tau = 0.975)$predicted
  Yc_coverage <- mean((Yc_l <= data_test[[tnam]]) & (data_test[[tnam]] <= Yc_u))
  
  # generalized random forest (Stefan's)
  grf_qf <- quantile_forest(as.matrix(data_train[,xnam,drop=FALSE]), data_train[,ynam], quantiles = c(0.025,0.1,0.5,0.9,0.975), 
                            num.trees = ntree, min.node.size = nodesize)
  Ygrf_l <- predict(grf_qf, as.matrix(data_test[,xnam,drop=FALSE]), quantiles = 0.025)
  Ygrf_u <- predict(grf_qf, as.matrix(data_test[,xnam,drop=FALSE]), quantiles = 0.975)
  Ygrf_coverage <- mean((Ygrf_l <= data_test[[tnam]]) & (data_test[[tnam]] <= Ygrf_u))
  
  grf_qf_latent <- quantile_forest(as.matrix(data_train[,xnam,drop=FALSE]), data_train[,tnam], quantiles = c(0.025,0.1,0.5,0.9,0.975), 
                                   num.trees = ntree, min.node.size = nodesize)
  Ygrf_latent_l <- predict(grf_qf_latent, as.matrix(data_test[,xnam,drop=FALSE]), quantiles = 0.025)
  Ygrf_latent_u <- predict(grf_qf_latent, as.matrix(data_test[,xnam,drop=FALSE]), quantiles = 0.975)
  Ygrf_latent_coverage <- mean((Ygrf_latent_l <= data_test[[tnam]]) & (data_test[[tnam]] <= Ygrf_latent_u))
  
  # quantile random forest (Meinshasen)
  qrf <- quantregForest(x=data_train[,xnam,drop=FALSE], y=data_train[,ynam], nodesize=nodesize, ntree=ntree)
  Yqrf_l <- predict(qrf, data_test[,xnam,drop=FALSE], what = 0.025)
  Yqrf_u <- predict(qrf, data_test[,xnam,drop=FALSE], what = 0.975)
  Yqrf_coverage <- mean((Yqrf_l <= data_test[[tnam]]) & (data_test[[tnam]] <= Yqrf_u))
  
  qrf_latent <- quantregForest(x=data_train[,xnam,drop=FALSE], y=data_train[,tnam], nodesize=nodesize, ntree=ntree)
  Yqrf_latent_l <- predict(qrf_latent, data_test[,xnam,drop=FALSE], what = 0.025)
  Yqrf_latent_u <- predict(qrf_latent, data_test[,xnam,drop=FALSE], what = 0.975)
  Yqrf_latent_coverage <- mean((Yqrf_latent_l <= data_test[[tnam]]) & (data_test[[tnam]] <= Yqrf_latent_u))
  
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

B <- 20*6
nodesize <- 0
ntree <- 1000
result <- list('nodesize'=rep(NA,B), 'crf'=rep(NA,B), 'qrf'=rep(NA,B), 'grf'=rep(NA,B), 
                        'qrf_oracle'=rep(NA,B), 'grf_oracle'=rep(NA,B))
for (t in 1:B) {
  if (t%%20 == 1) {
    nodesize <- nodesize + 15
  }
  
  tmp <- one_run(t, data.na, xnam, ynam, inam, nodesize, ntree)
  
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
ggsave(paste0(dataset, "_coverage_nodesize.pdf"), width = 5, height = 5, path = "./results/")
