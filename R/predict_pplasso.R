#' Predictions of a ppLasso or gr_ppLasso object
#'
#' Return the model predictions of a \code{ppLasso} or \code{gr_ppLasso} object
#'
#' @param fit a \code{ppLasso} or \code{gr_ppLasso} object.
#'
#' @param data an `dataframe` or `list` object that contains the variables for prediction.
#'
#' @param Z.char names of covariates in `data` as vector of character strings.
#'
#' @param prov.char name of provider IDs variable in `data` as a character string.
#'
#' @param lambda values of the regularization parameter lambda at which predictions are requested. For values of lambda not in the sequence of fitted models, linear interpolation is used.
#'
#' @param which indices of the penalty parameter lambda at which predictions are required. By default, all indices are returned. If lambda is specified, this will override which.
#'
#' @param type type of prediction: 
#'   * `link`: linear predictors
#'   * `response`: fitted values (i.e., `exp(link)/(1 + exp(link))`)
#'   * `class`: the binomial outcome with the largest probability
#'   * `vars`: the indices for the non-zero coefficients
#'   * `nvars`: the number of non-zero coefficients
#'
#' @param ...
#'
#' @exportS3Method predict ppLasso
#'
#' @examples
#' data(GLM_Data)
#' data <- GLM_Data$data
#' Y.char <- GLM_Data$Y.char
#' prov.char <- GLM_Data$prov.char
#' Z.char <- GLM_Data$Z.char
#' fit <- pp.lasso(data, Y.char, Z.char, prov.char)
#' predict(fit, data, Z.char, prov.char, lambda = fit$lambda, type = "response")[1:10, 1:5]
#' predict(fit, data, Z.char, prov.char, lambda = 0.001, type = "class")[1:10]
#' predict(fit, data, Z.char, prov.char, lambda = 0.04, type = "vars")
#'

predict.ppLasso <- function(fit, data, Z.char, prov.char, lambda, which = 1:length(fit$lambda),
                            type = c("link", "response", "class", "vars", "nvars"),  ...){
  beta <- coef.ppLasso(fit, lambda = lambda, which = which, drop = FALSE)$beta
  gamma <- coef.ppLasso(fit, lambda = lambda, which = which, drop = FALSE)$gamma

  if (type == "vars"){
    return(drop(apply(beta != 0, 2, FUN = which)))
  }

  if (type == "nvars") {
    v <- as.list(apply(beta != 0, 2, FUN = which))
    nvars <- sapply(v, length)
    return(nvars)
  }

  # predict response
  if (missing(data) | is.null(data)) {
    stop("Must supply data for predictions", call. = FALSE)
  }

  Prov.id <- data[, prov.char, drop = F]
  obs.prov.effect <- as.matrix(apply(Prov.id, 1, function(x) gamma[x, ]))

  if (ncol(obs.prov.effect) == 1){
    eta <- as.matrix(data[, Z.char, drop = F]) %*% beta + obs.prov.effect
  } else {
    eta <- as.matrix(data[, Z.char, drop = F]) %*% beta + t(obs.prov.effect)
  }
  
  if (type == "link") {
    return(eta)
  }
  
  if (type == "response") {
    pred.prob <- plogis(eta)
    return(pred.prob)
  }

  if (type == "class") {
    return(1 * (eta > 0))
  }
}


#' @rdname predict.ppLasso
#'
#' @param fit a \code{ppLasso} or \code{gr_ppLasso}.
#'
#' @param data an `dataframe` or `list` object that contains the variables for prediction.
#'
#' @param Z.char names of covariates in `data` as vector of character strings.
#'
#' @param prov.char name of provider IDs variable in `data` as a character string.
#'
#' @param lambda values of the regularization parameter lambda at which predictions are requested. For values of lambda not in the sequence of fitted models, linear interpolation is used.
#'
#' @param which indices of the penalty parameter lambda at which predictions are required. By default, all indices are returned. If lambda is specified, this will override which.
#'
#' @param type type of prediction: 
#'   * `link`: linear predictors
#'   * `response`: fitted values (i.e., `exp(link)/(1 + exp(link))`)
#'   * `class`: the binomial outcome with the largest probability
#'   * `vars`: the indices for the non-zero coefficients
#'   * `nvars`: the number of non-zero coefficients
#'   * `groups`: the indices for the non-zero groups
#'   * `ngroups`: the number of non-zero coefficients
#'   * `beta.norm`: L2 norm of the coefficients in each group
#'
#'
#' @param ...
#'
#' @exportS3Method predict gr_ppLasso
#'
#' @examples
#' data(GLM_Data)
#' data <- GLM_Data$data
#' Y.char <- GLM_Data$Y.char
#' prov.char <- GLM_Data$prov.char
#' Z.char <- GLM_Data$Z.char
#' group <- GLM_Data$group
#' fit <- grp.lasso(data, Y.char, Z.char, prov.char, group = group)
#' predict(fit, data, Z.char, prov.char, lambda = fit$lambda, type = "response")[1:10, 1:5]
#' predict(fit, data, Z.char, prov.char, lambda = 0.001, type = "class")[1:10]
#' predict(fit, data, Z.char, prov.char, lambda = 0.04, type = "vars")
#' predict(fit, data, Z.char, prov.char, lambda = 0.04, type = "groups")

predict.gr_ppLasso <- function(fit, data, Z.char, prov.char, lambda, which = 1:length(fit$lambda),
                               type = c("link", "response", "class", "vars", "groups", "nvars", "ngroups", "beta.norm"),  ...){
  beta <- coef.gr_ppLasso(fit, lambda = lambda, which = which, drop = FALSE)$beta
  gamma <- coef.gr_ppLasso(fit, lambda = lambda, which = which, drop = FALSE)$gamma

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

  if (type == "beta.norm"){
    return(drop(apply(beta, 2, function(x) tapply(x, fit$group, function(x){sqrt(sum(x^2))}))))
  }

  # predict response
  if (missing(data) | is.null(data)) {
    stop("Must supply data for predictions", call. = FALSE)
  }

  Prov.id <- data[, prov.char, drop = F]
  obs.prov.effect <- as.matrix(apply(Prov.id, 1, function(x) gamma[x, ]))

  if (ncol(obs.prov.effect) == 1){
    eta <- as.matrix(data[, Z.char, drop = F]) %*% beta + obs.prov.effect
  } else {
    eta <- as.matrix(data[, Z.char, drop = F]) %*% beta + t(obs.prov.effect)
  }
  
  if (type == "link") {
    return(eta)
  }

  if (type == "response") {
    pred.prob <- plogis(eta)
    return(pred.prob)
  }

  if (type == "class") {
    return(1 * (eta > 0))
  }
}

