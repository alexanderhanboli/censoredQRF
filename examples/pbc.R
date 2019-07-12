##------------------------------------------------------------
## Compare RF-SRC to Cox regression
## Illustrates C-index and Brier score measures of performance
## assumes "pec" and "survival" libraries are loaded
##------------------------------------------------------------
setwd("~/Projects/generalizedForest/examples")

source("metrics.R")
source("help_functions.R")
source("crf_km.R")
library(generalizedForest)
library(quantregForest)
library(ggplot2)
library(randomForestSRC)
if (library("survival", logical.return = TRUE)
    & library("prodlim", logical.return = TRUE)
    & library("grf", logical.return = TRUE))
{
  ## data, formula specifications
  data(pbc, package = "randomForestSRC")
  pbc.na <- na.omit(pbc) ##remove NA's
  surv.f <- as.formula(Surv(days, status) ~ .)
  xnam <- colnames(pbc.na)[-which(colnames(pbc.na) %in% c('status', 'days'))]
  ## compute out-of-bag C-index for cox regression and compare to rfsrc
  B <- 20
  nodesize <- 20
  ntree <- 2000
  
  cat("out-of-bag RSF Analysis ...", "\n")
  rfsrc.err <- sapply(1:B, function(b) {
    if (b%%5 == 0) cat("rfsrc bootstrap:", b, "\n")
    train <- sample(1:nrow(pbc.na), 0.8*nrow(pbc.na), replace = FALSE)
    rfsrc.obj <- rfsrc(surv.f, pbc.na[train, ], nodesize=nodesize, ntree = ntree)
    if (is.list(rfsrc.obj)) {
      1-randomForestSRC:::cindex(pbc.na$days[-train],
                                   pbc.na$status[-train],
                                   predict(rfsrc.obj, pbc.na[-train, ])$predicted)
    } else NA
  })
  
  cat("out-of-bag Cox Analysis ...", "\n")
  cox.err <- sapply(1:B, function(b) {
    if (b%%5 == 0) cat("cox bootstrap:", b, "\n")
    train <- sample(1:nrow(pbc.na), 0.8*nrow(pbc.na), replace = FALSE)
    cox.obj <- tryCatch({coxph(surv.f, pbc.na[train, ])}, error=function(ex){NULL})
    cox.pred <- tryCatch({survfit(cox.obj, newdata = pbc.na[-train, ])}, error=function(ex){NULL})
    cox.median <- summary(cox.pred)$table[,"median"]
    if (is.list(cox.obj)) {
      c(
        metrics(pbc.na[-train, 'days'], cox.median),
        'C' = 1-randomForestSRC:::cindex(pbc.na$days[-train],
                                   pbc.na$status[-train],
                                   predict(cox.obj, pbc.na[-train, ]))
      )
    } else NA
  })
  
  cat("out-of-bag wagerGRF Analysis ...", "\n")
  grf.err <- sapply(1:B, function(b) {
    if (b%%5 == 0) cat("grf bootstrap:", b, "\n")
    train <- sample(1:nrow(pbc.na), 0.8*nrow(pbc.na), replace = FALSE)
    grf.obj <- quantile_forest(X = pbc.na[train,xnam], Y = pbc.na[train,'days'], 
                               quantiles = 0.5, num.trees = ntree, min.node.size = nodesize)
    grf.pred <- predict(grf.obj, pbc.na[-train, xnam], quantiles=0.5)
    if (is.list(grf.obj)) {
      c(
        metrics(pbc.na[-train, 'days'], grf.pred),
        'C' = randomForestSRC:::cindex(pbc.na$days[-train],
                                 pbc.na$status[-train],
                                 grf.pred)
      )
    } else NA
  })
  
  cat("out-of-bag crf Analysis ...", "\n")
  crf.err <- sapply(1:B, function(b) {
    if (b%%5 == 0) cat("crf bootstrap:", b, "\n")
    train <- sample(1:nrow(pbc.na), 0.8*nrow(pbc.na), replace = FALSE)
    fmla <- as.formula(paste("days ~ ", paste(xnam, collapse= "+")))
    pred <- crf(fmla, ntree=ntree, nodesize=nodesize, pbc.na[train, ], 
                pbc.na[-train, ], 'days', 'status', 0.5)$predicted
    c(
      metrics(pbc.na[-train, 'days'], pred),
      'C' = randomForestSRC:::cindex(pbc.na$days[-train],
                               pbc.na$status[-train],
                               pred)
    )
  })
  
  cat("\n\tOOB C-index\n\n")
  cat("\tRSF : ", mean(rfsrc.err, na.rm = TRUE), "+-", sd(rfsrc.err, na.rm = TRUE), "\n")
  cat("\tCox regression : ", mean(cox.err['C',], na.rm = TRUE), "+-", sd(cox.err['C',], na.rm = TRUE), "\n")
  cat("\tGRF : ", mean(grf.err['C',], na.rm = TRUE), "+-", sd(grf.err['C',], na.rm = TRUE), "\n")
  cat("\tCRF : ", mean(crf.err['C',], na.rm = TRUE), "+-", sd(crf.err['C',], na.rm = TRUE), "\n")
  
  cat("\n\tOOB MSE\n\n")
  cat("\tCox regression : ", mean(cox.err['mse',], na.rm = TRUE), "+-", sd(cox.err['mse',], na.rm = TRUE), "\n")
  cat("\tGRF : ", mean(grf.err['mse',], na.rm = TRUE), "+-", sd(grf.err['mse',], na.rm = TRUE), "\n")
  cat("\tCRF : ", mean(crf.err['mse',], na.rm = TRUE), "+-", sd(crf.err['mse',], na.rm = TRUE), "\n")
  
  cat("\n\tOOB MAD\n\n")
  cat("\tCox regression : ", mean(cox.err['mad',], na.rm = TRUE), "+-", sd(cox.err['mad',], na.rm = TRUE), "\n")
  cat("\tGRF : ", mean(grf.err['mad',], na.rm = TRUE), "+-", sd(grf.err['mad',], na.rm = TRUE), "\n")
  cat("\tCRF : ", mean(crf.err['mad',], na.rm = TRUE), "+-", sd(crf.err['mad',], na.rm = TRUE), "\n")
  
  cat("\n\tOOB Bias\n\n")
  cat("\tCox regression : ", mean(cox.err['bias',], na.rm = TRUE), "+-", sd(cox.err['bias',], na.rm = TRUE), "\n")
  cat("\tGRF : ", mean(grf.err['bias',], na.rm = TRUE), "+-", sd(grf.err['bias',], na.rm = TRUE), "\n")
  cat("\tCRF : ", mean(crf.err['bias',], na.rm = TRUE), "+-", sd(crf.err['bias',], na.rm = TRUE), "\n")
}