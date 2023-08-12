loss.grp.lasso <- function(Y.i, yhat.i){
  yhat.i[yhat.i < 0.00001] <- 0.00001
  yhat.i[yhat.i > 0.99999] <- 0.99999
  
  if (is.matrix(yhat.i) == T) {
    val <- matrix(NA, nrow = nrow(yhat.i), ncol = ncol(yhat.i))
    if (sum(Y.i == 1)) {  # if all 1 or all zero, then we only need calculate one of the following "if"
      val[Y.i == 1,] <- - 2 * log(yhat.i[Y.i == 1, , drop = FALSE])
    }
    if (sum(Y.i == 0)){
      val[Y.i == 0,] <- -2 * log(1 - yhat.i[Y.i == 0, , drop = FALSE])
    }
  } else {
    val <- double(length(Y.i))  # a zero vector
    if (sum(Y.i == 1)) {  # if all 1 or all zero, then we only need calculate one of the following "if"
      val[Y.i == 1] <- - 2 * log(yhat.i[Y.i == 1])
    }
    if (sum(Y.i == 0)){
      val[Y.i == 0] <- -2 * log(1 - yhat.i[Y.i == 0])
    }
  }
  return(val)
}


loss.Disc.Surv <- function(Y.i, yhat.i){
  yhat.i[yhat.i < 0.00001] <- 0.00001
  yhat.i[yhat.i > 0.99999] <- 0.99999
  
  if (is.matrix(yhat.i) == T) {
    val <- matrix(NA, nrow = nrow(yhat.i), ncol = ncol(yhat.i))
    if (sum(Y.i == 1)) {  # if all 1 or all zero, then we only need calculate one of the following "if"
      val[Y.i == 1,] <- - 2 * log(yhat.i[Y.i == 1, , drop = FALSE])
    }
    if (sum(Y.i == 0)){
      val[Y.i == 0,] <- -2 * log(1 - yhat.i[Y.i == 0, , drop = FALSE])
    }
  } else {
    val <- double(length(Y.i))  # a zero vector
    if (sum(Y.i == 1)) {  # if all 1 or all zero, then we only need calculate one of the following "if"
      val[Y.i == 1] <- - 2 * log(yhat.i[Y.i == 1])
    }
    if (sum(Y.i == 0)){
      val[Y.i == 0] <- -2 * log(1 - yhat.i[Y.i == 0])
    }
  }
  return(val)
}

loss.strat_cox <- function(delta.obs, y.hat, ID, total=TRUE){
  y.hat <- as.matrix(y.hat)
  revCumsum.strat <- function(temp.ID){
    temp.y.hat <- y.hat[which(ID == temp.ID), , drop=FALSE] # "time" still keep increasing order
    temp.rsk <- apply(temp.y.hat, 2, function(x) rev(cumsum(rev(exp(x)))))
    return(temp.rsk)
  }
  rsk <- do.call(rbind, lapply(1:nrow(unique(ID)), function(i) revCumsum.strat(i)))
  
  if (total == TRUE) {
    return(-2 * (crossprod(delta.obs, y.hat) - crossprod(delta.obs, log(rsk))))
  } else {
    return(-2 * (y.hat[delta.obs == 1, , drop = FALSE] - log(rsk)[delta.obs == 1, , drop = FALSE]))
  }
}


