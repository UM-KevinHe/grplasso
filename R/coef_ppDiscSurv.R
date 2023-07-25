#' Extract coefficients of a ppDiscSurv object
#'
#' Return the model coefficients of a \code{ppDiscSurv} object
#'
#' @param fit a \code{ppDiscSurv} object.
#'
#' @param lambda values of the regularization parameter lambda at which coefficients are requested. For values of lambda not in the sequence of fitted models, linear interpolation is used.
#'
#' @param which indices of the penalty parameter lambda at which predictions are required. By default, all indices are returned. If lambda is specified, this will override which.
#'
#' @param drop whether to keep coefficient names
#'
#' @param ...
#'
#' @exportS3Method coef ppDiscSurv
#'
#' @examples
#' data(Surv_Data)
#' data <- Surv_Data$data
#' Event.char <- Surv_Data$Event.char
#' prov.char <- Surv_Data$prov.char
#' Z.char <- Surv_Data$Z.char
#' Time.char <- Surv_Data$Time.char
#' fit <- pp.DiscSurv(data, Event.char, prov.char, Z.char, Time.char)
#' coef(fit, lambda = fit$lambda)$beta[, 1:10]
#' coef(fit, lambda = fit$lambda)$gamma[, 1:10]


coef.ppDiscSurv <- function(fit, lambda, which=1:length(fit$lambda), drop = TRUE, ...) {
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
    gamma <- (1 - w) * fit$gamma[, l, drop = FALSE] + w * fit$gamma[, r, drop = FALSE]
    alpha <- (1 - w) * fit$alpha[, l, drop = FALSE] + w * fit$alpha[, r, drop = FALSE]
    colnames(beta) <- round(lambda, 4)
    colnames(gamma) <- round(lambda, 4)
    colnames(alpha) <- round(lambda, 4)
  } else {  #specify lambda value as index
    beta <- fit$beta[, which, drop = FALSE]
    gamma <- fit$gamma[, which, drop = FALSE]
    alpha <- fit$alpha[, which, drop = FALSE]
  }
  if (drop == TRUE){
    beta <- drop(beta)
    gamma <- drop(gamma)
    alpha <- drop(alpha)
    coef <- list(alpha = alpha,
                 gamma = gamma, 
                 beta = beta)
  } else {
    coef <- list(alpha = alpha,
                 gamma = gamma, 
                 beta = beta)
  }
  return(coef)
}



