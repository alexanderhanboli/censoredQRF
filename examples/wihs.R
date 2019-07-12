##------------------------------------------------------------
## Compare RF-SRC to Cox regression
## Illustrates C-index and Brier score measures of performance
## assumes "pec" and "survival" libraries are loaded
##------------------------------------------------------------
source("crf_function.R")
source("metrics.R")
if (library("survival", logical.return = TRUE)
    & library("prodlim", logical.return = TRUE)
    & library("grf", logical.return = TRUE))
{
  ## data, formula specifications
  data(wihs, package = "randomForestSRC")
  data.na <- na.omit(veteran) ##remove NA's
  surv.f <- as.formula(Surv(time, status) ~ .)
  xnam <- colnames(data.na)[-which(colnames(data.na) %in% c('status', 'time'))]
  ## compute out-of-bag C-index for cox regression and compare to rfsrc
  B <- 50
  nodesize <- 10
  ntree <- 2000
  
  cat("out-of-bag RSF Analysis ...", "\n")
  rfsrc.err <- sapply(1:B, function(b) {
    if (b%%5 == 0) cat("rfsrc bootstrap:", b, "\n")
    train <- sample(1:nrow(data.na), 0.8*nrow(data.na), replace = FALSE)
    rfsrc.obj <- rfsrc(surv.f, data.na[train, ], nodesize=nodesize, ntree = ntree)
    if (is.list(rfsrc.obj)) {
      1-randomForestSRC:::cindex(data.na$time[-train],
                                 data.na$status[-train],
                                 predict(rfsrc.obj, data.na[-train, ])$predicted)
    } else NA
  })
  
  cat("out-of-bag Cox Analysis ...", "\n")
  cox.err <- sapply(1:B, function(b) {
    if (b%%5 == 0) cat("cox bootstrap:", b, "\n")
    train <- sample(1:nrow(data.na), 0.8*nrow(data.na), replace = FALSE)
    cox.obj <- tryCatch({coxph(surv.f, data.na[train, ])}, error=function(ex){NULL})
    cox.pred <- tryCatch({survfit(cox.obj, newdata = data.na[-train, ])}, error=function(ex){NULL})
    cox.median <- tryCatch({summary(cox.pred)$table[,"median"]}, error=function(ex){NULL})
    if (is.list(cox.obj)) {
      c(
        metrics(data.na[-train, 'time'], cox.median),
        'C' = 1-randomForestSRC:::cindex(data.na$time[-train],
                                         data.na$status[-train],
                                         predict(cox.obj, data.na[-train, ]))
      )
    } else NA
  })
  
  cat("out-of-bag wagerGRF Analysis ...", "\n")
  grf.err <- sapply(1:B, function(b) {
    if (b%%5 == 0) cat("grf bootstrap:", b, "\n")
    train <- sample(1:nrow(data.na), 0.8*nrow(data.na), replace = FALSE)
    grf.obj <- quantile_forest(X = data.na[train,xnam], Y = data.na[train,'time'], 
                               quantiles = 0.5, num.trees = ntree, min.node.size = nodesize)
    grf.pred <- predict(grf.obj, data.na[-train, xnam], quantiles=0.5)
    if (is.list(grf.obj)) {
      c(
        metrics(data.na[-train, 'time'], grf.pred),
        'C' = randomForestSRC:::cindex(data.na$time[-train],
                                       data.na$status[-train],
                                       grf.pred)
      )
    } else NA
  })
  
  cat("out-of-bag crf Analysis ...", "\n")
  crf.err <- sapply(1:B, function(b) {
    if (b%%5 == 0) cat("crf bootstrap:", b, "\n")
    train <- sample(1:nrow(data.na), 0.8*nrow(data.na), replace = FALSE)
    fmla <- as.formula(paste("time ~ ", paste(xnam, collapse= "+")))
    pred <- crf(fmla, ntree=ntree, nodesize=nodesize, data.na[train, ], 
                data.na[-train, ], 'time', 'status', 0.5)$predicted
    c(
      metrics(data.na[-train, 'time'], pred),
      'C' = randomForestSRC:::cindex(data.na$time[-train],
                                     data.na$status[-train],
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