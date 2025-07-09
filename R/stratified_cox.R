#' Fit a penalized stratified cox model
#'
#' Main function for fitting a penalized stratified cox model. 
#'
#' @param data an `dataframe` or `list` object that contains the variables in the model.
#'
#' @param Event.char name of the event indicator in `data` as a character string. Event indicator should be a 
#' binary variable with 1 indicating that the event has occurred and 0 indicating (right) censoring.
#'
#' @param prov.char name of stratum indicator in `data` as a character string.
#' If "prov.char" is not specified, all observations are are considered to be from the same stratum.
#'
#' @param Z.char names of covariates in `data` as vector of character strings.
#'
#' @param Time.char name of the follow up time in `data` as a character string.
#'
#' @param group a vector describing the grouping of the coefficients. If there are coefficients to be included in the model without being penalized, assign them to group 0 (or "0").
#'
#' @param group.multiplier A vector of values representing multiplicative factors by which each covariate's penalty is to be multiplied. Default is a vector of 1's.
#'
#' @param standardize logical flag for x variable standardization, prior to fitting the model sequence. The coefficients are always returned on the original scale. Default is `standardize=TRUE`.
#'
#' @param lambda a user supplied lambda sequence. Typical usage is to have the program compute its own lambda sequence based on `nlambda` and `lambda.min.ratio`.
#'
#' @param nlambda the number of lambda values. Default is 100.
#'
#' @param lambda.min.ratio the fraction of the smallest value for lambda with `lambda.max` (smallest lambda for which all coefficients are zero) on log scale. Default is 1e-03.
#'
#' @param lambda.early.stop whether the program stop before running the entire sequence of lambda. Early stop based on the ratio of deviance for models under two successive lambda. Default is `FALSE`.
#'
#' @param nvar.max number of maximum selected variables. Default is the number of all covariates.
#' 
#' @param group.max number of maximum selected groups. Default is the number of all groups.
#'
#' @param stop.loss.ratio if `lambda.early.stop = TRUE`, the ratio of loss for early stopping. Default is 1e-3.
#
#' @param tol convergence threshold. For each lambda, the program will stop if the maximum change of covariate coefficient is smaller than `tol`. Default is 1e-4.
#'
#' @param max.each.iter maximum number of iterations for each lambda. Default is 1e4.
#'
#' @param max.total.iter maximum number of iterations for entire path. Default is `max.each.iter` * `nlambda`.
#'
#' @param actSet whether to use the active method for variable selection. Default is TRUE.
#'
#' @param actIter if `actSet = TRUE`, the maximum number of iterations for a new updated active set. Default is `max.each.iter` (i.e. we will update the current active set until convergence ).
#'
#' @param actGroupNum if `actSet = TRUE`, the maximum number of variables that can be selected into the new active set for each time when the active set is updated. Default is number of groups.
#'
#' @param actSetRemove if `actSet = TRUE`, whether we remove the zero coefficients from the current active set. Default is FALSE.
#'
#' @param returnX whether return the standardized design matrix. Default is FALSE.
#'
#' @param trace.lambda whether display the progress for fitting the entire path. Default is FALSE.
#' 
#' @param ... extra arguments to be passed to function.
#'
#'
#' @return An object with S3 class \code{strat_cox}.
#'
#' \item{beta}{the fitted matrix of covariate coefficients.
#' The number of rows is equal to the number of coefficients,
#' and the number of columns is equal to nlambda.}
#' 
#' \item{group}{a vector describing the grouping of the coefficients.}
#'
#' \item{lambda}{the sequence of `lambda` values in the path.}
#' 
#' \item{loss}{the likelihood of the fitted model at each value of `lambda`.}
#' 
#' \item{linear.predictors}{the linear predictors of the fitted model at each value of `lambda`.}
#'
#' \item{df}{the estimates of effective number of selected variables all the points along the regularization path.}
#'
#' \item{iter}{the number of iterations until convergence at each value of `lambda`.}
#'
#' @export
#'
#' @seealso \code{\link{coef}}, \code{\link{plot}} function.
#'
#' @examples
#' data(ContTime)
#' data <- ContTime$data
#' Event.char <- ContTime$Event.char
#' prov.char <- ContTime$prov.char
#' Z.char <- ContTime$Z.char
#' Time.char <- ContTime$Time.char
#' fit <- Strat.cox(data, Event.char, Z.char, Time.char, prov.char, group = c(1, 2, 2, 3, 3))
#' fit$beta[, 1:5]
#'
#' @importFrom Rcpp evalCpp
#'
#' @details
#' The model is fit by Newton method and coordinate descent method.
#'
#' @references
#' K. He, J. Kalbfleisch, Y. Li, and et al. (2013) Evaluating hospital readmission rates in dialysis facilities; adjusting for hospital effects.
#' \emph{Lifetime Data Analysis}, \strong{19}: 490-512.
#' \cr


Strat.cox <- function(data, Event.char, Z.char, Time.char, prov.char, group = 1:length(Z.char), group.multiplier,
                      standardize = T, lambda, nlambda = 100, lambda.min.ratio = 1e-3, lambda.early.stop = FALSE,
                      nvar.max = p, group.max = length(unique(group)), stop.loss.ratio = 1e-3, tol = 1e-4, 
                      max.each.iter = 1e4, max.total.iter = (max.each.iter * nlambda), actSet = TRUE, 
                      actIter = max.each.iter, actGroupNum = sum(unique(group) != 0), actSetRemove = F,
                      returnX = FALSE, trace.lambda = FALSE,...){
  
  if (missing(prov.char)){ #single intercept
    warning("Provider information not provided. All data is assumed to originate from a single provider!", call. = FALSE)
    ID <- matrix(1, nrow = nrow(data))
    colnames(ID) <- "intercept"
    data <- data[order(data[, Time.char]), ]
  } else {
    # re-order original data based on observed time (stratified by provider)
    data <- data[order(data[, prov.char], data[, Time.char]), ]
    #recode prov.ID as {1, 2, 3, .....}
    unique.prov <- unique(data[, prov.char])
    prov.ref <- cbind(1:length(unique.prov), unique.prov)
    colnames(prov.ref) <- c("New.ID", prov.char)
    ID <- as.matrix(data[, prov.char])  #stratum indicator
    colnames(ID) <- prov.char
    ID <- merge(ID, prov.ref, by = prov.char)[, 2, drop = F] # now, ID has been recoded as 1, 2, 3, ....
    colnames(ID) <- prov.char
  }

  n.each_prov <- table(ID)

  initial.group <- group
  if (standardize == T){
    std.Z <- newZG.Std.grplasso(data, Z.char, group, group.multiplier)
  } else {
    std.Z <- newZG.Unstd.grplasso(data, Z.char, group, group.multiplier)
  }
  Z <- std.Z$std.Z[, , drop = F]
  group <- std.Z$g  
  group.multiplier <- std.Z$m 
  
  delta.obs <- data[, Event.char]
  time <- data[, Time.char]
  p <- ncol(Z)
  nvar.max <- as.integer(nvar.max)
  group.max <- as.integer(group.max)
  
  beta <- rep(0, ncol(Z)) #initial value of beta
  
  if (missing(lambda)) {
    if (nlambda < 2) {
      stop("nlambda must be at least 2", call. = FALSE)
    } else if (nlambda != round(nlambda)){
      stop("nlambda must be a positive integer", call. = FALSE)
    }
    lambda.fit <- set.lambda.cox(delta.obs, Z, time, ID, beta, group, group.multiplier, n.each_prov,
                                 nlambda = nlambda, lambda.min.ratio = lambda.min.ratio)
    lambda.seq <- lambda.fit$lambda.seq
    beta <- lambda.fit$beta
  } else {
    nlambda <- length(lambda)  # Note: lambda can be a single value
    lambda.seq <- as.vector(sort(lambda, decreasing = TRUE))
  }
  
  K <- as.integer(table(group)) #number of features in each group
  K0 <- as.integer(if (min(group) == 0) K[1] else 0)
  K1 <- as.integer(if (min(group) == 0) cumsum(K) else c(0, cumsum(K)))
  
  initial.active.group <- -1
  if (actSet == TRUE){
    if (K0 == 0){
      initial.active.group <- which(K == min(K))[1] - 1
    }
  } else {
    actIter <- max.each.iter
  }

  fit <- StratCox_lasso(delta.obs, Z, n.each_prov, beta, K0, K1, lambda.seq, lambda.early.stop, stop.loss.ratio, 
                        group.multiplier, max.total.iter, max.each.iter, tol, initial.active.group, nvar.max, 
                        group.max, trace.lambda, actSet, actIter, actGroupNum, actSetRemove)
  
  beta <- fit$beta
  eta <- fit$Eta
  df <- fit$Df
  iter <- fit$iter
  loss <- fit$loss
  
  # Eliminate saturated lambda values
  ind <- !is.na(iter)
  lambda <- lambda.seq[ind]
  beta <- beta[, ind, drop = FALSE]
  loss <- loss[ind]
  eta <- eta[, ind, drop = FALSE]
  df <- df[ind]
  iter <- iter[ind]
  
  if (iter[1] == max.total.iter){
    stop("Algorithm failed to converge for any values of lambda", call. = FALSE)
  }
  if (sum(iter) == max.total.iter){
    warning("Algorithm failed to converge for all values of lambda", call. = FALSE)
  }
  
  
  # Original scale
  beta <- unorthogonalize(beta, std.Z$std.Z, group)
  rownames(beta) <- colnames(Z)
  if (std.Z$reorder == TRUE){  # original order of beta
    beta <- beta[std.Z$ord.inv, , drop = F]
  }
  if (standardize == T) {
    original.beta <- matrix(0, nrow = length(std.Z$scale), ncol = ncol(beta))
    original.beta[std.Z$nz, ] <- beta / std.Z$scale[std.Z$nz]
    beta <- original.beta
  }
  
  # Names
  dimnames(beta) <- list(Z.char, round(lambda, digits = 4))
  colnames(eta) <- round(lambda, digits = 4)
  
  result <- structure(list(beta = beta,
                           group = factor(initial.group),
                           lambda = lambda,
                           loss = loss,
                           linear.predictors = eta,  # rescale beta will not change the linear predictors (Z matrix is also standardized)
                           df = df,
                           iter = iter,
                           group.multiplier = group.multiplier),
                      class = "strat_cox")  #define a list for prediction
  if (returnX == TRUE){
    result$returnX <- std.Z
  }
  return(result)
}








