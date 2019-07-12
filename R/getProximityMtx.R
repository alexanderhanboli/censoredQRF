getProximityMtx <- function(object, newdata=NULL, train_approx = FALSE, ...) {

  if (!inherits(object, "generalizedForest"))
      stop("object is not of class generalizedForest")

  if (is.null(newdata)) newdata <- object[["data"]]
  n <- NROW(newdata)
  train_data <- object[["data"]]
  response <- object$y
  train_nodes <- attr(predict(object,newdata=train_data,nodes=TRUE),"nodes")
  object$nodes = train_nodes

  rf_knn <- predict(object, newdata, proximity=2, nodes=TRUE)
  proximityMtxTrain <- NULL
  if (train_approx) {
      rf_knn_train <- predict(object, train_data, proximity=2, nodes=TRUE)
      proximityMtxTrain <- rf_knn_train$proximity
  }
  pred <- rf_knn$predicted
  proximityMtx <- rf_knn$proximity
  colnames(proximityMtx) <- row.names(train_data)
  row.names(proximityMtx) <- row.names(newdata)
  nodes <- attr(rf_knn, "nodes")

  return(list(
    'proximityMtxTrain' = proximityMtxTrain,
    'proximityMtx' = proximityMtx,
    'nodes' = nodes,
    'predicted' = pred
  ))
}
