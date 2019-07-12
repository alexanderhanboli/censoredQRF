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

## data, formula specifications
dataset <- 'engel'
if (dataset == 'pbc') {
  data(pbc, package = "randomForestSRC")
  data.na <- na.omit(pbc) ##remove NA's
  xnam <- colnames(data.na)[-which(colnames(data.na) %in% c('status', 'days'))]
  yid <- which(colnames(data.na) == 'days')
  colnames(data.na)[yid] <- 'y'
  ynam <- 'y'
  inam <- 'status'
  print(mean(data.na[[inam]]))
  #surv.f <- as.formula(Surv(y, status) ~ .)
} else if (dataset == 'veteran') {
  data(veteran, package = "randomForestSRC")
  data.na <- na.omit(veteran) ##remove NA's
  xnam <- colnames(data.na)[-which(colnames(data.na) %in% c('status', 'time'))]
  yid <- which(colnames(data.na) == 'time')
  colnames(data.na)[yid] <- 'y'
  ynam <- 'y'
  inam <- 'status'
  print(mean(data.na[[inam]]))
  #surv.f <- as.formula(Surv(y, status) ~ .)
} else if (dataset == 'wihs') {
  data(wihs, package = "randomForestSRC")
  data.na <- na.omit(wihs) ##remove NA's
  xnam <- colnames(data.na)[-which(colnames(data.na) %in% c('status', 'time'))]
  yid <- which(colnames(data.na) == 'time')
  colnames(data.na)[yid] <- 'y'
  ynam <- 'y'
  inam <- 'status'
  data.na$status[data.na$status == 2] <- 1
  print(mean(data.na[[inam]]))
  #surv.f <- as.formula(Surv(y, status) ~ .)
} else if (dataset == 'vdv') {
  data(vdv, package = "randomForestSRC")
  data.na <- na.omit(vdv) ##remove NA's
  xnam <- colnames(data.na)[-which(colnames(data.na) %in% c('Censoring', 'Time'))]
  yid <- which(colnames(data.na) == 'Time')
  colnames(data.na)[yid] <- 'y'
  ynam <- 'y'
  inam <- 'Censoring'
  print(mean(data.na[[inam]]))
  #surv.f <- as.formula(Surv(y, status) ~ .)
} else if (dataset == 'uis') {
  data(uis, package = "quantreg")
  data.na <- na.omit(uis) ##remove NA's
  data.na$ID <- NULL
  data.na$Y <- NULL
  xnam <- colnames(data.na)[-which(colnames(data.na) %in% c('CENSOR', 'TIME'))]
  yid <- which(colnames(data.na) == 'TIME')
  colnames(data.na)[yid] <- 'y'
  data.na$y <- log(data.na$y)
  iid <- which(colnames(data.na) == 'CENSOR')
  colnames(data.na)[iid] <- 'status'
  ynam <- 'y'
  inam <- 'status'
  print(mean(data.na[[inam]]))
  #surv.f <- as.formula(Surv(y, status) ~ .)
} else if (dataset == 'engel') {
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
  inam <- 'status'
  print(mean(data.na[[inam]]))
}

one_run <- function(b, data, xnam, ynam, inam, tau, nodesize, ntree) {
  cat("bootstrap step:", b, "\n")
  
  train <- sample(1:nrow(data), nrow(data), replace = TRUE)
  data_train <- data[train, ]
  data_test <- data[-train, ][data[-train, inam] == 1, ]
  
  # build generalizedForest model
  fmla <- as.formula(paste(ynam, "~", paste(xnam, collapse= "+")))
  
  Yc <- crf.km(fmla, ntree = ntree, nodesize = 10*nodesize, data_train = data_train, data_test = data_test, 
               yname = ynam, iname = inam, tau = tau)$predicted
  
  # generalized random forest (Stefan's)
  grf_qf <- quantile_forest(data_train[,xnam,drop=FALSE], data_train[,ynam], quantiles = tau, 
                            num.trees = ntree, min.node.size = nodesize)
  Ygrf <- predict(grf_qf, data_test[,xnam,drop=FALSE], quantiles = tau)
  
  # quantile random forest (Meinshasen)
  qrf <- quantregForest(x=data_train[,xnam,drop=FALSE], y=data_train[,ynam], nodesize=5*nodesize, ntree=ntree)
  Yqrf <- predict(qrf, data_test[,xnam,drop=FALSE], what = tau)
  
  # survival forest
  surv_rf <- rfsrc(Surv(y, status) ~ ., data = data_train, ntree = ntree, nodesize = nodesize)
  Ysurv <- predict(object = surv_rf, newdata = data_test)$predicted
  
  # results
  return(
    list(
      'crf'=metrics(data_test[[ynam]], Yc, data_test[[ynam]], tau),
      'qrf'=metrics(data_test[[ynam]], Yqrf, data_test[[ynam]], tau),
      'grf'=metrics(data_test[[ynam]], Ygrf, data_test[[ynam]], tau),
      
      'crf_c'=randomForestSRC:::cindex(data_test[[ynam]], data_test[[inam]], Yc),
      'qrf_c'=randomForestSRC:::cindex(data_test[[ynam]], data_test[[inam]], Yqrf),
      'rsf_c'=1-randomForestSRC:::cindex(data_test[[ynam]], data_test[[inam]], Ysurv),
      'grf_c'=randomForestSRC:::cindex(data_test[[ynam]], data_test[[inam]], Ygrf)
    )
  )
}

B <- 20
nodesize <- 5
ntree <- 1000
tau <- 0.5
quantile_result <- list('crf'=rep(NA,B), 'qrf'=rep(NA,B), 'grf'=rep(NA,B))
cindex_result <- list('crf'=rep(NA,B), 'qrf'=rep(NA,B), 'rsf'=rep(NA,B), 'grf'=rep(NA,B))
for (t in 1:B) {
  tmp <- one_run(t, data.na, xnam, ynam, inam, tau, nodesize, ntree)
  
  quantile_result$crf[t] <- tmp$crf['quantile_loss']
  quantile_result$qrf[t] <- tmp$qrf['quantile_loss']
  quantile_result$grf[t] <- tmp$grf['quantile_loss']
  
  cindex_result$crf[t] <- tmp$crf_c
  cindex_result$qrf[t] <- tmp$qrf_c
  cindex_result$rsf[t] <- tmp$rsf_c
  cindex_result$grf[t] <- tmp$grf_c
} 

# boxplot
require(reshape2)

dd <- as.data.frame(cindex_result)
ggplot(data = melt(dd), aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=variable)) + labs(x = "Model", y = paste("C-index, tau =", tau), fill = "Model")
ggsave(paste0("pbc_cindex_result_", 10*tau, ".pdf"), width = 5, height = 5, path = "./results/")

dd <- as.data.frame(quantile_result)
ggplot(data = melt(dd), aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=variable)) + labs(x = "Model", y = paste("Quantile loss, tau =", tau), fill = "Model")
ggsave(paste0("pbc_quantile_loss_result_", 10*tau, ".pdf"), width = 5, height = 5, path = "./results/")