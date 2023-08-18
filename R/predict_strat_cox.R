#' Predictions of a strat_cox object
#'
#' Return the model predictions of a \code{strat_cox} object
#'
#' @param fit a \code{strat_cox} object.
#'
#' @param data an `dataframe` or `list` object that contains the variables for prediction.
#'
#' @param Z.char names of covariates in `data` as vector of character strings.
#'
#' @param lambda values of the regularization parameter lambda at which predictions are requested. For values of lambda not in the sequence of fitted models, linear interpolation is used.
#'
#' @param which indices of the penalty parameter lambda at which predictions are required. By default, all indices are returned. If lambda is specified, this will override which.
#'
#' @param type type of prediction: 
#'   * `link`: linear predictors
#'   * `response`: risk (i.e., `exp(link)`)
#'   * `vars`: the indices for the non-zero coefficients
#'   * `nvars`: the number of non-zero coefficients
#'   * `groups`: the indices for the non-zero groups
#'   * `ngroups`: the number of non-zero coefficients
#'
#' @param ...
#'
#' @exportS3Method predict strat_cox
#'
#' @examples
#' data(Cox_Data)
#' data <- Cox_Data$data
#' Event.char <- Cox_Data$Event.char
#' prov.char <- Cox_Data$prov.char
#' Z.char <- Cox_Data$Z.char
#' Time.char <- Cox_Data$Time.char
#' fit <- Strat.cox(data, Event.char, Z.char, Time.char, prov.char, group = c(1, 2, 2, 3, 3))
#' predict(fit, data, Z.char, lambda = fit$lambda, type = "response")[1:5, 1:5]

predict.strat_cox <- function(fit, data, Z.char, lambda, which = 1:length(fit$lambda),
                               type = c("link", "response", "vars", "nvars", "groups", "ngroups"), ...){
  beta <- coef.strat_cox(fit, lambda = lambda, which = which, drop = FALSE)

  if (type == "vars"){
    return(drop(apply(beta != 0, 2, FUN = which)))
  }
  
  if (type == "groups") {
    if (ncol(beta) == 1) {
      return(unique(fit$group[beta != 0]))
    } else {
      return(drop(apply(beta != 0, 2, function(x) unique(fit$group[x]))))
    }
  }
  
  if (type == "nvars") {
    v <- as.list(apply(beta != 0, 2, FUN = which))
    nvars <- sapply(v, length)
    return(nvars)
  }
  
  if (type == "ngroups") {
    g <- as.data.frame(apply(beta != 0, 2, function(x) unique(fit$group[x])))
    ngroups <- sapply(g, length)
    if (length(ngroups) == 1){
      names(ngroups) <- colnames(beta)
    }
    return(ngroups)
  }
  
  # predict response
  if (missing(data) | is.null(data)) {
    stop("Must supply data for predictions", call. = FALSE)
  }
  
  eta <- as.matrix(data[, Z.char, drop = F]) %*% beta
  if (type == "link") {
    return(eta)
  }
  
  if (type == "response") {
    return(exp(eta))
  }
}