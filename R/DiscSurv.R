#' Fit a penalized discrete survival model (without provider information)
#'
#' Main function for fitting a penalized discrete survival model without provider information
#'
#' @param data an `dataframe` or `list` object that contains the variables in the model.
#'
#' @param Event.char name of the event indicator in `data` as a character string.

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
#' @return An object with S3 class \code{DiscSurv}.
#'
#' \item{beta}{the fitted matrix of covariate coefficients.
#' The number of rows is equal to the number of coefficients,
#' and the number of columns is equal to nlambda.}
#'
#' \item{alpha}{the fitted value of logit-transformed baseline hazard.}
#'
#' \item{lambda}{the sequence of `lambda` values in the path.}
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
#' data(Surv_Data)
#' data <- Surv_Data$data
#' Event.char <- Surv_Data$Event.char
#' Z.char <- Surv_Data$Z.char
#' Time.char <- Surv_Data$Time.char
#' fit <- DiscSurv(data, Event.char, Z.char, Time.char) # fit a discrete survival model without any given provider information.
#' fit$beta[, 1:5] # covariate coefficient
#' fit$alpha[, 1:5] #time effect
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

DiscSurv <- function(data, Event.char, Z.char, Time.char, lambda, nlambda = 100, lambda.min.ratio = 1e-4, 
                     penalize.x = rep(1, length(Z.char)), penalized.multiplier, lambda.early.stop = FALSE, 
                     nvar.max = p, stop.dev.ratio = 1e-3, bound = 10.0, backtrack = FALSE, tol = 1e-4, 
                     max.each.iter = 1e4, max.total.iter = (max.each.iter * nlambda), actSet = TRUE,
                     actIter = max.each.iter, actVarNum = sum(penalize.x == 1), actSetRemove = F, returnX = FALSE,
                     trace.lambda = FALSE, threads = 1, MM = FALSE, return.transform.data = FALSE, ...){
  count.alpha <- length(unique(data[, Time.char]))
  timepoint.increase <- sort(unique(data[, Time.char]))
  for (i in 1:count.alpha){
    data[, Time.char][which(data[, Time.char] == timepoint.increase[i])] <- i
  }
  
  time.ref <- as.data.frame(cbind(1:count.alpha, timepoint.increase))
  colnames(time.ref) <- c("time.indicator", "time")
  
  max.timepoint <- length(timepoint.increase) 
  
  pseudo.group <- 1:length(Z.char)    
  if (min(penalize.x) == 0){
    pseudo.group[penalize.x == 0] = 0 
  }

  fit <- survival::survfit(survival::Surv(data[, Time.char], data[, Event.char]) ~ 1)
  sum.failure <- fit$n.event
  KM.baseline.hazard <- c(fit$cumhaz[1], fit$cumhaz[2:length(fit$cumhaz)] - fit$cumhaz[1:(length(fit$cumhaz) - 1)])

  std.Z <- newZG.Unstd.grplasso(data, Z.char, pseudo.group, penalized.multiplier)
  Z <- std.Z$std.Z[, , drop = F]
  pseudo.group <- std.Z$g
  penalized.multiplier <- std.Z$m

  delta.obs <- data[, Event.char]
  time <- data[, Time.char]
  p <- ncol(Z)
  nvar.max <- as.integer(nvar.max)
  
  KM.baseline.hazard.add.small <- KM.baseline.hazard
  KM.baseline.hazard.add.small[which(KM.baseline.hazard.add.small == 0)] <- 1e-10
  alpha <- log(KM.baseline.hazard.add.small/(1 - KM.baseline.hazard.add.small))
  beta <- rep(0, ncol(Z))
  
  if (missing(lambda)) {
    if (nlambda < 2) {
      stop("nlambda must be at least 2", call. = FALSE)
    } else if (nlambda != round(nlambda)){
      stop("nlambda must be a positive integer", call. = FALSE)
    }
    lambda.fit <- set.lambda.DiscSurv(delta.obs, Z, time, alpha , beta, pseudo.group, penalized.multiplier, 
                                      nlambda = nlambda, lambda.min.ratio = lambda.min.ratio)
    lambda.seq <- lambda.fit$lambda.seq
    beta <- lambda.fit$beta
    alpha <- lambda.fit$alpha
  } else {
    nlambda <- length(lambda)  
    lambda.seq <- as.vector(sort(lambda, decreasing = TRUE))
  }
  
  K <- as.integer(table(pseudo.group)) 
  K0 <- as.integer(if (min(pseudo.group) == 0) K[1] else 0)
  K1 <- as.integer(if (min(pseudo.group) == 0) cumsum(K) else c(0, cumsum(K)))
  
  initial.active.variable <- -1
  if (actSet == TRUE){
    if (K0 == 0){ 
      initial.active.variable <- which(K == min(K))[1] - 1 
    }
  } else {  ## if we don't use active set method, then the initial active set should contain all penalized variables
    actIter <- max.each.iter
  }
  
  fit <- DiscSurv_lasso(delta.obs, max.timepoint, Z, time, alpha, beta, K0, K1, sum.failure, lambda.seq, 
                        penalized.multiplier, max.total.iter, max.each.iter, tol, backtrack, MM, bound, 
                        initial.active.variable, nvar.max, trace.lambda, threads, actSet, actIter, actVarNum,
                        actSetRemove)
  
  alpha <- fit$gamma # notations in .cpp functions are different
  beta <- fit$beta
  eta <- fit$Eta
  df <- fit$Df
  iter <- fit$iter
  
  # Eliminate saturated lambda values
  ind <- !is.na(iter)
  lambda <- lambda.seq[ind]
  beta <- beta[, ind, drop = FALSE]
  alpha <- alpha[, ind, drop = FALSE] #time effect
  eta <- eta[, ind, drop = FALSE]
  df <- df[ind]
  iter <- iter[ind]
  
  if (sum(iter) == max.total.iter){
    warning("Algorithm failed to converge for all values of lambda", call. = FALSE)
  }
  
  # Original scale
  beta <- unorthogonalize(beta, std.Z$std.Z, pseudo.group)
  rownames(beta) <- colnames(Z)
  if (std.Z$reorder == TRUE){  # original order of beta
    beta <- beta[std.Z$ord.inv, , drop = F]
  }
  
  # Names
  dimnames(beta) <- list(Z.char, round(lambda, digits = 4))
  if (nrow(alpha) == 1 & length(lambda.seq) == 1){
    alpha <- t(alpha)
  }
  
  dimnames(alpha) <- list(paste0("[Time: ", time.ref[, 2], "]"), round(lambda, digits = 4)) #no intercept
  colnames(eta) <- round(lambda, digits = 4)
  
  result <- structure(list(beta = beta,
                           alpha = alpha,
                           lambda = lambda,
                           linear.components = eta,  #eta = gamma +  Z * beta
                           df = df,
                           iter = iter,
                           penalized.multiplier = penalized.multiplier,
                           penalize.x = penalize.x, 
                           time.ref = time.ref),
                      class = "DiscSurv")  #define a list for prediction
  
  if (return.transform.data == TRUE){  #return a dataframe with recoded time
    result$transform.data <- data
  }  
  
  if (returnX == TRUE){
    result$returnX <- std.Z
  }
  return(result)
}
