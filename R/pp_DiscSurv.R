#' fit a group penalized generalized regression model
#'
#' Fit a group penalized generalized regression model via coordinate descent method:
#'
#' @param data an `dataframe` or `list` object that contains the variables in the model.
#'
#' @param Event.char name of the event indicator in `data` as a character string.
#'
#' @param prov.char name of provider IDs variable in `data` as a character string.
#'
#' @param Z.char names of covariates in `data` as vector of character strings.
#'
#' @param Time.char name of the observation time in `data` as a character string.
#'
#' @param lambda a user supplied lambda sequence. Typical usage is to have the program compute its own lambda sequence based on `nlambda` and `lambda.min.ratio`.
#'
#' @param nlambda the number of lambda values. Default is 100.
#'
#' @param lambda.min.ratio the fraction of the smallest value for lambda with `lambda.max` (smallest lambda for which all coefficients are zero) on log scale. Default is 1e-04.
#'
#' @param penalize.x  a vector indicates whether the corresponding covariate will be penalized. If equals 0, variable is unpenalized, else is penalized. Default is a vector of 1's (all covariates are penalized).
#'
#' @param penalized.multiplier A vector of values representing multiplicative factors by which each covariate's penalty is to be multiplied. Default is a vector of 1's.
#'
#' @param lambda.early.stop whether the program stop before running the entire sequence of lambda. Early stop based on the ratio of deviance for models under two successive lambda. Default is `FALSE`.
#'
#' @param nvar.max number of maximum selected variables. Default is the number of all covariates.
#'
#' @param stop.dev.ratio if `lambda.early.stop = TRUE`, the ratio of deviance for early stopping. Default is 1e-3.
#
#' @param bound a positive number to avoid inflation of provider effect. Default is 10.
#'
#' @param backtrack for updating the provider effect, whether to use the "backtracking line search" with Newton method.
#'
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
#' @param actSetRemove if `actSet = TRUE`, whether we remove the zero coefficients from the current active set. Default is FALSE.
#'
#' @param returnX whether return the standardized design matrix. Default is FALSE.
#'
#' @param trace.lambda whether display the progress for fitting the entire path. Default is FALSE.
#'
#' @param threads number of cores that are used for parallel computing.
#'
#' @param MM whether we use the "Majorize-Minimization" algorithm to optimize the objective function.
#'
#' @param ... extra arguments to be passed to function.
#'
#'
#' @return An object with S3 class `ppDiscSurv`.
#'
#' \item{beta}{the fitted matrix of covariate coefficients.
#' The number of rows is equal to the number of coefficients,
#' and the number of columns is equal to nlambda.}
#'
#' \item{alpha}{the fitted value of logit-transformed baseline hazard.}
#'
#' \item{gamma}{the fitted value of provider effects.}
#'
#' \item{lambda}{the sequence of `lambda` values in the path.}
#'
#' \item{df}{the estimates of effective number of selected variables all the points along the regularization path.}
#'
#' \item{iter}{the number of iterations until convergence at each value of `lambda`.}
#'
#' @export
#'
#' @seealso \code{\link{coef}} function.
#'
#' @examples
#'\dontrun{
#' data(Surv_Data)
#' Event.char <- 'status'
#' prov.char <- 'Prov.ID'
#' Z.char <- c("Z1", "Z2", "Z3", "Z4", "Z5")
#' Time.char <- "time"
#' fit <- pp.DiscSurv(Surv_Data, Event.char, prov.char, Z.char, Time.char)
#' }
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


pp.DiscSurv <- function(data, Event.char, prov.char, Z.char, Time.char, lambda, nlambda = 100,
                     lambda.min.ratio = 1e-4, penalize.x = rep(1, length(Z.char)), penalized.multiplier,
                     lambda.early.stop = FALSE, nvar.max = p, stop.dev.ratio = 1e-3, bound = 10.0, backtrack = FALSE,
                     tol = 1e-4, max.each.iter = 1e4, max.total.iter = (max.each.iter * nlambda), actSet = TRUE,
                     actIter = max.each.iter, actVarNum = sum(penalize.x == 1), actSetRemove = F, returnX = FALSE,
                     trace.lambda = FALSE, threads = 1, MM = FALSE, ...){

  # Convert the observed time in discrete intervals
  # "time" is converted into order index (integers)
  count.gamma <- length(unique(data[, Time.char]))
  timepoint.increase <- sort(unique(data[, Time.char]))
  # new time start from 1, and time points are {1, 2, 3, ...}
  for (i in 1:count.gamma){
    data[, Time.char][which(data[, Time.char] == timepoint.increase[i])] <- i
  }
  max.timepoint <- length(timepoint.increase)  # the number of gamma that we need

  data <- data[order(factor(data[, prov.char])), ]
  ID <- as.matrix(data[, prov.char]) # ID vector
  colnames(ID) <- prov.char

  pseudo.group <- 1:length(Z.char)    #each variable forms a group
  if (min(penalize.x) == 0){
    pseudo.group[penalize.x == 0] = 0 # set unpenalized variables
  }

  # failure at each time point
  fit <- survival::survfit(survival::Surv(data[, Time.char], data[, Event.char]) ~ 1)
  sum.failure <- fit$n.event
  KM.baseline.hazard <- c(fit$cumhaz[1], fit$cumhaz[2:length(fit$cumhaz)] - fit$cumhaz[1:(length(fit$cumhaz) - 1)])

  failure.each.center <- tapply(data$status, data$Prov.ID, sum)


  # !!! "standardize = T" may cause problems in transforming gamma and alpha back, current we only consider "standardize = F"

  #if (standardize == T){
  #  std.Z <- newZG.Std.grplasso(data, Z.char, pseudo.group, penalized.multiplier)
  #  Z <- std.Z$std.Z[, , drop = F]  # standardized covariate matrix
  #  pseudo.group <- std.Z$g  # new group order
  #  penalized.multiplier <- std.Z$m # new group multiplier
  #} else {
    std.Z <- newZG.Unstd.grplasso(data, Z.char, pseudo.group, penalized.multiplier)
    Z <- std.Z$std.Z[, , drop = F]
    pseudo.group <- std.Z$g
    penalized.multiplier <- std.Z$m
  #}

  delta.obs <- data[, Event.char]  #each observation's event variable
  time <- data[, Time.char]
  p <- ncol(Z)
  nvar.max <- as.integer(nvar.max)

  # gamma start from KM estimator
  KM.baseline.hazard.add.small <- KM.baseline.hazard
  KM.baseline.hazard.add.small[which(KM.baseline.hazard.add.small == 0)] <- 1e-10
  gamma <- log(KM.baseline.hazard.add.small/(1 - KM.baseline.hazard.add.small))
  beta <- rep(0, ncol(Z))
  n.prov <- sapply(split(delta.obs, ID), length)
  alpha.prov <- rep(0, length(n.prov))

  #####---check-----#####
  ##if (missing(lambda)) {
  ##  if (nlambda < 2) {
  ##    stop("nlambda must be at least 2", call. = FALSE)
  ##  } else if (nlambda != round(nlambda)){
  ##    stop("nlambda must be a positive integer", call. = FALSE)
  ##  }
  ##  lambda.fit <- set.lambda.Surv2(delta.obs, Z, time, ID, gamma, beta, alpha.prov, prov.char,
  ##                                 pseudo.group, penalized.multiplier, nlambda = nlambda,
  ##                                 lambda.min.ratio = lambda.min.ratio)
  ##  lambda.seq <- lambda.fit$lambda.seq
  ##  beta <- lambda.fit$beta
  ##  alpha <- lambda.fit$alpha
  ##  gamma <- lambda.fit$gamma
  ##} else {
  ##  nlambda <- length(lambda)  # Note: lambda can be a single value
  ##  lambda.seq <- as.vector(sort(lambda, decreasing = TRUE))
  ##}
  #####---check-----#####

  lambda.seq <- 0

  K <- as.integer(table(pseudo.group)) #number of features in each group
  K0 <- as.integer(if (min(pseudo.group) == 0) K[1] else 0)
  K1 <- as.integer(if (min(pseudo.group) == 0) cumsum(K) else c(0, cumsum(K)))

  initial.active.variable <- -1 ## "-1" means this variable can be any value and we will not use it in our cpp function.
  if (actSet == TRUE){
    if (K0 == 0){ #all variables are penalized
      initial.active.variable <- which(K == min(K))[1] - 1  ## which is actually the first variable
    }
  } else {  ## if we don't use active set method, then the initial active set should contain all penalized variables
    actIter <- max.each.iter
  }

  fit <- pp_DiscSurv_lasso(delta.obs, max.timepoint, Z, n.prov, time, gamma, beta, alpha.prov, K0, K1, sum.failure, failure.each.center,
                           lambda.seq, penalized.multiplier, max.total.iter, max.each.iter, tol, backtrack, MM, bound, initial.active.variable,
                           nvar.max, trace.lambda, threads, actSet, actIter, actVarNum, actSetRemove)

  gamma <- fit$gamma
  beta <- fit$beta
  alpha <- fit$alpha
  eta <- fit$Eta
  df <- fit$Df
  iter <- fit$iter

  # Eliminate saturated lambda values
  ind <- !is.na(iter)
  lambda <- lambda.seq[ind]
  beta <- beta[, ind, drop = FALSE]
  gamma <- gamma[, ind, drop = FALSE]
  eta <- eta[, ind, drop = FALSE]
  df <- df[ind]
  iter <- iter[ind]

  #if (iter[1] == max.total.iter){
  #  stop("Algorithm failed to converge for any values of lambda", call. = FALSE)
  #}
  if (sum(iter) == max.total.iter){
    warning("Algorithm failed to converge for all values of lambda", call. = FALSE)
  }


  # Original scale
  beta <- unorthogonalize(beta, std.Z$std.Z, pseudo.group)
  rownames(beta) <- colnames(Z)
  if (std.Z$reorder == TRUE){  # original order of beta
    beta <- beta[std.Z$ord.inv, , drop = F]
  }
  #if (standardize == T) {
  #  unstandardize.para <- unstandardize(beta, gamma, std.Z)
  #  beta <- unstandardize.para$beta
  #  gamma <- unstandardize.para$gamma
  #}

  # Names
  dimnames(beta) <- list(Z.char, round(lambda, digits = 4))
  if (nrow(gamma) == 1 & length(lambda.seq) == 1){
    gamma <- t(gamma)
    alpha <- t(alpha)
  }

  dimnames(gamma) <- list(paste0("T_", 1:count.gamma), round(lambda, digits = 4))
  dimnames(alpha) <- list(names(n.prov), round(lambda, digits = 4))

  colnames(eta) <- round(lambda, digits = 4)

  result <- structure(list(beta = beta,
                           alpha = gamma,
                           gamma = alpha,
                           lambda = lambda,
                           linear.components = eta,  #eta = alpha +  Z * beta
                           df = df,
                           iter = iter,
                           penalized.multiplier = penalized.multiplier,
                           penalize.x = penalize.x),
                      class = "ppDiscSurv")  #define a list for prediction
  if (returnX == TRUE){
    result$returnX <- std.Z
  }
  return(result)
}
