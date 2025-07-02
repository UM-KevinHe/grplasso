#' Extract coefficients of a DiscSurv object
#'
#' Return the model coefficients of a \code{DiscSurv} object
#'
#' @param fit a \code{DiscSurv} object.
#'
#' @param lambda values of the regularization parameter lambda at which coefficients are requested. For values of lambda not in the sequence of fitted models, linear interpolation is used.
#'
#' @param which indices of the penalty parameter lambda at which predictions are required. By default, all indices are returned. If lambda is specified, this will override which.
#'
#' @param drop whether to keep coefficient names
#'
#' @param ...
#'
#' @exportS3Method coef DiscSurv
#'
#' @examples
#' data(DiscTime)
#' data <- DiscTime$data
#' Event.char <- DiscTime$Event.char
#' Z.char <- DiscTime$Z.char
#' Time.char <- DiscTime$Time.char
#' fit <- DiscSurv(data, Event.char, Z.char, Time.char)
#' coef(fit, lambda = fit$lambda)$alpha[, 1:10]
#' coef(fit, lambda = fit$lambda)$beta[, 1:10]



coef.DiscSurv <- function(fit, lambda, which=1:length(fit$lambda), drop = TRUE, ...) {
  if (!missing(lambda)) {
    if (any(lambda > max(fit$lambda) | lambda < min(fit$lambda))){
      stop('lambda must lie within the range of the fitted coefficient path', call.=FALSE)
    }
    ind <- approx(fit$lambda, seq(fit$lambda), lambda)$y
    l <- floor(ind)
    r <- ceiling(ind)
    w <- ind %% 1
    # linearly interpolate between two lambda
    beta <- (1 - w) * fit$beta[, l, drop = FALSE] + w * fit$beta[, r, drop = FALSE]
    alpha <- (1 - w) * fit$alpha[, l, drop = FALSE] + w * fit$alpha[, r, drop = FALSE]
    colnames(beta) <- round(lambda, 4)
    colnames(alpha) <- round(lambda, 4)
  } else {  #specify lambda value as index
    beta <- fit$beta[, which, drop = FALSE]
    alpha <- fit$alpha[, which, drop = FALSE]
  }
  if (drop == TRUE){
    beta <- drop(beta)
    alpha <- drop(alpha)
    coef <- list(alpha = alpha, 
                 beta = beta)
  } else {
    coef <- list(alpha = alpha, 
                 beta = beta)
  }
  return(coef)
}



