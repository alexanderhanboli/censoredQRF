"generalizedForest" <- function(formula, data = NULL, ntree = 100, nodesize = 5,
  proximity = 0, ...) {

    cl <- match.call()
    cl[[1]] <- as.name("generalizedForest")
    grf <- randomForest(formula, data=data, ntree=ntree, nodesize=nodesize,
      proximity=proximity, ...)
    grf[["call"]] = cl
    grf[["data"]] = data

    class(grf) <- c("generalizedForest", "randomForest")

    return(grf)
}
