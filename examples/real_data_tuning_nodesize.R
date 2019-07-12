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
dataset <- 'uis'
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
}

one_run <- function(b, data, xnam, ynam, inam, tau, nodesize, ntree) {
  cat("bootstrap step:", b, "\n")
  
  train <- sample(1:nrow(data), nrow(data), replace = TRUE)
  data_train <- data[train, ]
  data_test <- data[-train, ][data[-train, inam] == 1, ]
  
  # build generalizedForest model
  fmla <- as.formula(paste(ynam, "~", paste(xnam, collapse= "+")))
  
  Yc1 <- crf.km(fmla, ntree = ntree, nodesize = 50*nodesize, data_train = data_train, data_test = data_test, 
               yname = ynam, iname = inam, tau = tau)$predicted
  
  Yc2 <- crf.km(fmla, ntree = ntree, nodesize = 40*nodesize, data_train = data_train, data_test = data_test, 
                yname = ynam, iname = inam, tau = tau)$predicted
  
  Yc3 <- crf.km(fmla, ntree = ntree, nodesize = 30*nodesize, data_train = data_train, data_test = data_test, 
                yname = ynam, iname = inam, tau = tau)$predicted
  
  Yc4 <- crf.km(fmla, ntree = ntree, nodesize = 20*nodesize, data_train = data_train, data_test = data_test, 
                yname = ynam, iname = inam, tau = tau)$predicted
  
  
  # results
  return(
    list(
      'crf1'=metrics(data_test[[ynam]], Yc1, data_test[[ynam]], tau),
      'crf2'=metrics(data_test[[ynam]], Yc2, data_test[[ynam]], tau),
      'crf3'=metrics(data_test[[ynam]], Yc3, data_test[[ynam]], tau),
      'crf4'=metrics(data_test[[ynam]], Yc4, data_test[[ynam]], tau)
    )
  )
}

B <- 20
nodesize <- 5
ntree <- 1000
tau <- 0.5
quantile_result <- list('crf1'=rep(NA,B), 'crf2'=rep(NA,B), 'crf3'=rep(NA,B), 'crf4'=rep(NA,B))
for (t in 1:B) {
  tmp <- one_run(t, data.na, xnam, ynam, inam, tau, nodesize, ntree)
  
  quantile_result$crf1[t] <- tmp$crf1['quantile_loss']
  quantile_result$crf2[t] <- tmp$crf2['quantile_loss']
  quantile_result$crf3[t] <- tmp$crf3['quantile_loss']
  quantile_result$crf4[t] <- tmp$crf4['quantile_loss']
} 

# boxplot
require(reshape2)

dd <- as.data.frame(quantile_result)
ggplot(data = melt(dd), aes(x=variable, y=value)) + 
  geom_boxplot(aes(fill=variable)) + labs(x = "Model", y = paste("Quantile loss, tau =", tau), fill = "Model")
ggsave(paste0("pbc_quantile_loss_result_", 10*tau, ".pdf"), width = 5, height = 5, path = "./results/")