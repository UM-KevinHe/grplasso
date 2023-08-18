#' Extract coefficients of a strat_cox object
#'
#' Return the model coefficients of a \code{strat_cox} object
#'
#' @param fit a \code{strat_cox} object.
#'
#' @param lambda values of the regularization parameter lambda at which coefficients are requested. For values of lambda not in the sequence of fitted models, linear interpolation is used.
#'
#' @param which indices of the penalty parameter lambda at which predictions are required. By default, all indices are returned. If lambda is specified, this will override which.
#'
#' @param drop whether to keep coefficient names
#'
#' @param ...
#'
#' @exportS3Method coef strat_cox
#'
#' @examples
#' data(Cox_Data)
#' data <- Cox_Data$data
#' Event.char <- Cox_Data$Event.char
#' prov.char <- Cox_Data$prov.char
#' Z.char <- Cox_Data$Z.char
#' Time.char <- Cox_Data$Time.char
#' fit <- Strat.cox(data, Event.char, Z.char, Time.char, prov.char, group = c(1, 2, 2, 3, 3))
#' coef(fit, lambda = fit$lambda)[, 1:5]

coef.strat_cox <- function(fit, lambda, which=1:length(fit$lambda), drop = TRUE, ...) {
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
    colnames(beta) <- round(lambda, 4)
  } else {  #specify lambda value as index
    beta <- fit$beta[, which, drop = FALSE]
  }
  if (drop == TRUE){
    beta <- drop(beta)
  }
  return(beta)
}