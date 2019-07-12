source("metrics.R")
source("help_functions.R")
source("Csurv.R")
library(generalizedForest)
#library(randomForestSRC)
library(survival)

G.compare <- function(fmla, ntree, nodesize, data_train, data_test, yname, iname, nnb = 1) {
  # build generalizedForest model
  grf <- generalizedForest(fmla, data = data_train, ntree = ntree, nodesize=nodesize, proximity = TRUE)
  # get proximity matrix
  proxSup <- getProximityMtx(grf, newdata = data_test)
  proxMtx <- proxSup$proximityMtx
  # censor forest
  n <- nrow(data_test)
  ntrain <- nrow(data_train)
  Ytrain <- data_train[[yname]]
  censorInd <- data_train[[iname]]
  
  G.beran <- matrix(nrow = NROW(data_test), ncol = 500)
  G.nnb <- matrix(nrow = NROW(data_test), ncol = 500)

  for (r in 1:NROW(data_test)) {
    # C survival estimate
    candidateY <- seq(0, 8, length.out = 500)
    
    # Beran
    boot.idx <- proxMtx[r, ] > 0
    denom <- sapply(Ytrain, function(x) {1*(Ytrain >= x)%*%proxMtx[r, ]})
    denom[denom == 0] <- 1
    base <- (1 - proxMtx[r, ]/denom)^(1 - censorInd)
    
    for (i in 1:length(candidateY)) {
      G.beran[r,i] <- C.surv(candidateY[i], Ytrain, base)
    }
    
    # KM nnb
    boot.idx <- order(proxMtx[r, ], decreasing = TRUE)[1:nnb]
    C.boot <- Ytrain[boot.idx]
    i.boot <- 1 - censorInd[boot.idx]
    C.km <- survfit(Surv(C.boot, i.boot) ~ 1, type = 'kaplan-meier')
    C.surv.nnb <- stepfun(C.km$time, c(1, C.km$surv))
    
    for (i in 1:length(candidateY)) {
      G.nnb[r,i] <- C.surv.nnb(candidateY[i])
    }
  }
    
  return(
    list(
      'x'=candidateY,
      'beran'=G.beran,
      'nnb'=G.nnb
    )
  )
}

# MAIN
n_test <- 4
sigma <- 0.3
rate <- 0.2

# test data (fixed)
Xtest <- c(0.2, 0.4, 0.6, 0.8)*2
Ytest <- exp(Xtest + rnorm(n_test, mean = 0, sd = sigma))
quantile_test <- exp(Xtest + qnorm(tau, 0, sigma))
data_test <- cbind.data.frame(Xtest, Ytest, rep(1, n_test))

n <- 5000
# training data
Xtrain <- sort(runif(n = n, min = 0, max = 2))
Ttrain <- exp(Xtrain + rnorm(n, mean = 0, sd = sigma))
ctrain <- rexp(n = n, rate = rate)
Ytrain <- pmin(Ttrain, ctrain)
censorInd <- 1*(Ttrain <= ctrain)
print(mean(censorInd))
data_train <- cbind.data.frame(Xtrain, Ytrain, censorInd)

# column names
colnames(data_train) <- c('x', 'y', 'ind')
colnames(data_test) <- c('x', 'y', 'ind')

# survival functions
g.comp <- G.compare(y~x, ntree = 1000, nodesize = n/10, 
                    data_train = data_train, data_test = data_test, 
                    yname = 'y', iname = 'ind', nnb = n/10)

# true
G.true <- pexp(g.comp$x, rate = rate, lower.tail = FALSE)

# plot
require(reshape2)

result <- as.data.frame(list('q'=g.comp$x, 'G.true'=G.true, 'G.beran'=g.comp$beran[1,], 'G.nnb'=g.comp$nnb[1,]))
dd <- melt(result, id.vars = 'q')
ggplot(data = dd, aes(x=q, y=value, colour=variable)) + 
  geom_line(size=1) +
  labs(x = "q", y = "Survival function value", fill = "q")
ggsave(paste0("G_comparison_high_contamination_", 1, "samplesize_", n, ".pdf"), width = 5, height = 5, path = "./results/")

for (i in 2:n_test){
  result <- as.data.frame(list('q'=g.comp$x, 'G.true'=G.true, 'G.beran'=g.comp$beran[i,], 'G.nnb'=g.comp$nnb[i,]))
  dd <- melt(result, id.vars = 'q')
  ggplot(data = dd, aes(x=q, y=value, colour=variable)) + 
    geom_line(size=1) +
    labs(x = "q", y = "Survival function value", fill = "q")
  ggsave(paste0("G_comparison_high_contamination_", i, "samplesize_", n, ".pdf"), width = 5, height = 5, path = "./results/")
}