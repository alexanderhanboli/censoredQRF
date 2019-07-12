source("metrics.R")
source("help_functions.R")
library(generalizedForest)
library(randomForestSRC)

crf <- function(fmla, ntree, nodesize, data_train, data_test, yname, iname, 
                tau, fixed_censoring = FALSE, fixed_c = NULL, debias = FALSE) {
  # build generalizedForest model
  grf <- generalizedForest(fmla, data = data_train, ntree = ntree, nodesize=nodesize, proximity = TRUE)
  # get proximity matrix
  proxSup <- getProximityMtx(grf, newdata = data_test)
  proxMtx <- proxSup$proximityMtx
  proxMtxTrain <- proxSup$proximityMtxTrain
  # censor forest
  Yrf <- proxSup$predicted # RF predictions
  Yc <- rep(NA, length(Yrf)) # to store new predictions
  Ytrain <- data_train[[yname]]
  ctrain <- fixed_c
  censorInd <- data_train[[iname]]
  Ytrain_sorted <- sort(Ytrain)
  candidateY <- seq(min(Ytrain), max(Ytrain), length.out = 2000)
  # candidateY <- Ytrain_sorted
  right_censoring <- TRUE
  # find minimum
  for (r in 1:NROW(data_test)) {
    Yc[r] <- candidateY[1]
    min_loss <- 100000
    for (lambda in candidateY) {
      if (right_censoring) {
        if (fixed_censoring) {
          loss <- proxMtx[r, ]%*%quantile_loss(Ytrain - pmin(ctrain, lambda), tau)
        } else if (debias == FALSE) {
          loss1 <- (1-censorInd)*quantile_loss(Ytrain - pmin(Ytrain, lambda), tau)
          loss2 <- censorInd*quantile_loss(Ytrain - lambda, tau)
          loss <- proxMtx[r, ]%*%(loss1 + loss2)
        } else {
          censoring_prob <- as.vector(proxMtx[r, ]%*%(1-censorInd))
          kappa <- 1 - pnorm(tau, mean = 1-censoring_prob, sd = sqrt(censoring_prob*(1-censoring_prob)/1))
          loss1 <- (1-censorInd)*quantile_loss(Ytrain - pmin(Ytrain, lambda), tau)
          loss2 <- kappa*censorInd*quantile_loss(Ytrain - lambda, tau)
          loss <- proxMtx[r, ]%*%(loss1 + loss2)
        }
      } else {
        #
      }
      if (loss < min_loss) {
        Yc[r] <- lambda
        min_loss <- loss
      }
    }
  }
  return(
    list(
      'predicted'=Yc,
      'c'=c
    )
  )
}