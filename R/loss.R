loss.grp.lasso <- function(Y.i, yhat.i){
  yhat.i[yhat.i < 0.00001] <- 0.00001
  yhat.i[yhat.i > 0.99999] <- 0.99999
  
  if (is.matrix(yhat.i) == TRUE) {
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
  
  if (is.matrix(yhat.i) == TRUE) {
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

loss.strat_cox <- function(delta.obs, y.hat, ID, weight = NULL, total = TRUE){
  y.hat <- as.matrix(y.hat)
  n <- nrow(y.hat)
  if (is.null(weight)) weight <- rep(1.0, n)
  weight <- as.vector(weight)

  revCumsum.strat <- function(temp.ID){
    idx <- which(ID == temp.ID)
    temp.y.hat <- y.hat[idx, , drop = FALSE]  # "time" still keeps increasing order
    temp.w    <- weight[idx]
    # weighted risk set: rev(cumsum(rev(w * exp(eta)))) for each lambda
    temp.rsk <- apply(temp.y.hat, 2, function(x) rev(cumsum(rev(temp.w * exp(x)))))
    return(temp.rsk)
  }
  rsk <- do.call(rbind, lapply(1:nrow(unique(ID)), function(i) revCumsum.strat(i)))

  if (total == TRUE) {
    # -2 * sum_i  w_i * delta_i * (eta_i - log(rsk_i))
    return(-2 * crossprod(weight * delta.obs, y.hat - log(rsk)))
  } else {
    # per-event weighted contribution (one row per event)
    idx1 <- which(delta.obs == 1)
    return(-2 * weight[idx1] * (y.hat[idx1, , drop = FALSE] - log(rsk)[idx1, , drop = FALSE]))
  }
}


