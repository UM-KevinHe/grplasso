#' Cross-validation for penalized discrete survival model (without provider information)
#'
#' Performs k-fold cross validation for penalized discrete survival model (without provider information) over a grid of values of regularization parameter lambda.
#'
#' @param data an `dataframe` or `list` object that contains the variables in the model.
#'
#' @param Event.char name of the event indicator in `data` as a character string. Event indicator should be a 
#' binary variable with 1 indicating that the event has occurred and 0 indicating (right) censoring.
#'
#' @param Z.char names of covariates in `data` as vector of character strings.
#'
#' @param Time.char name of the follow up time in `data` as a character string.
#'
#' @param penalize.x  a vector indicates whether the corresponding covariate will be penalized, as in \code{pp.DiscSurv} function.
#'
#' @param nfolds the number of cross-validation folds. Default is 10.
#'
#' @param seed the seed of the random number generator in order to obtain reproducible results.
#'
#' @param fold a vector that specifies the fold that observations belongs to. By default the observations are randomly assigned.
#'
#' @param trace.cv \code{cv.DiscSurv} will provide user with the progress of cross validation if `trace.cv = TRUE`. Default is FALSE.
#'
#' @param ... extra arguments to be passed to function.
#'
#' @return An object with S3 class \code{cv.DiscSurv}.
#'
#' \item{cve}{the error for each value of lambda, averaged across the cross-validation folds.}
#'
#' \item{cvse}{the estimated standard error associated with each value of for cve.}
#'
#' \item{lambda}{the sequence of regularization parameter values along which the cross-validation error was calculated.}
#'
#' \item{fit}{the fitted \code{pp.DiscSurv} object for the whole data.}
#'
#' \item{fold}{the fold assignments for cross-validation for each observation}
#'
#' \item{min}{the index of lambda corresponding to lambda.min.}
#'
#' \item{lambda.min}{the value of lambda with the minimum cross-validation error.}
#'
#' @export
#'
#' @examples
#' data(DiscTime)
#' data <- DiscTime$data
#' Event.char <- DiscTime$Event.char
#' Z.char <- DiscTime$Z.char
#' Time.char <- DiscTime$Time.char
#' cv.fit <- cv.DiscSurv(data, Event.char, Z.char, Time.char, nfolds = 10, trace.cv = T)
#' cv.fit$cve
#' cv.fit$lambda.min
#'
#' @references
#' K. He, J. Kalbfleisch, Y. Li, and et al. (2013) Evaluating hospital readmission rates in dialysis facilities; adjusting for hospital effects.
#' \emph{Lifetime Data Analysis}, \strong{19}: 490-512.
#' \cr



cv.DiscSurv <- function(data, Event.char, Z.char, Time.char, penalize.x = rep(1, length(Z.char)),
                        ..., nfolds = 10, seed, fold, trace.cv = FALSE){
  fit.args <- list(...)
  fit.args$data <- data
  fit.args$Event.char <- Event.char
  fit.args$Z.char <- Z.char
  fit.args$Time.char <- Time.char
  fit.args$penalize.x <- penalize.x
  fit.args$returnX <- TRUE
  fit.args$return.transform.data <- TRUE
  
  fit <- do.call("DiscSurv", fit.args) 
  
  data <- fit$transform.data 
  time.ref <- fit$time.ref
  
  # get standardized Z
  ZG <- fit$returnX
  Z <- ZG$std.Z
  delta.obs <- as.matrix(data[, Event.char])
  returnX <- list(...)$returnX
  if (is.null(returnX) || !returnX){ 
    fit$returnX <- NULL
  } 
  
  #setup folds
  if (!missing(seed)) {
    set.seed(seed)
  }
  n.obs <- length(delta.obs)
  original.count.alpha <- length(unique(data[, Time.char]))
  
  try.times <- 100
  for (s in 1:try.times){
    if (missing(fold)){ # make sure each fold contains same proportion of censor and failure
      ind1 <- which(delta.obs == 1)
      ind0 <- which(delta.obs == 0)
      n1 <- length(ind1)
      n0 <- length(ind0)
      fold1 <- 1:n1 %% nfolds  # assign observations to different folds
      fold0 <- (n1 + 1:n0) %% nfolds
      fold1[fold1 == 0] <- nfolds # set fold "0" to fold max
      fold0[fold0 == 0] <- nfolds
      fold <- integer(n.obs) 
      fold[delta.obs == 1] <- sample(fold1)
      fold[delta.obs == 0] <- sample(fold0)
    } else {
      nfolds <- max(fold)
    }
    
    # make sure each training data (9/10 data) contains all timepoints
    len <- rep(NA, nfolds)  
    for (i in 1:nfolds){
      len[i] <- length(unique(data[fold != i, Time.char])) != original.count.alpha
    }
    
    if (sum(len) == 0) {
      break
    }
    
    if (s == try.times){
      stop("Having too many time points can impede the proper functioning of the cross-validation procedure.
           Please attempt to merge some adjacent time points.", call. = FALSE)
    }
  }
  
  data.small <- cbind(data[, c(Time.char, Event.char)], fold)
  expand.fold <- discSurv::dataLong(dataShort = data.small, timeColumn = Time.char, 
                                    eventColumn = Event.char, timeAsFactor = TRUE)$fold
  
  # Do Cross-Validation
  E <- Y <- matrix(NA, nrow = sum(data[, Time.char]), ncol = length(fit$lambda)) # stored as expanded matrix
  cv.args <- list(...)
  cv.args$lambda <- fit$lambda
  cv.args$group <- ZG$g
  cv.args$penalized.multiplier <- ZG$m
  
  for (i in 1:nfolds) {
    if (trace.cv == TRUE){
      cat("Starting CV fold #", i, sep = "", "...\n")
    }
    res <- cvf.DiscSurv(i, data, Event.char, Z.char, Time.char, fold, cv.args)
    Y[expand.fold == i, 1:res$nl] <- res$yhat
    E[expand.fold == i, 1:res$nl] <- res$loss
  }
  
  # Eliminate saturated lambda values, if any
  ind <- which(apply(is.finite(E), 2, all))
  E <- E[, ind, drop=FALSE]
  Y <- Y[, ind]
  lambda <- fit$lambda[ind]
  
  cve <- apply(E, 2, mean)  # mean cross entropy loss within each lambda value
  cvse <- apply(E, 2, sd) / sqrt(n.obs) # standardized loss
  min <- which.min(cve)  #find index of lambda with minimum cve
  
  result <- structure(list(cve = cve, 
                           cvse = cvse, 
                           lambda = lambda, 
                           fit = fit, #model with entire data
                           fold = fold,
                           min = min, 
                           lambda.min = lambda[min]),
                      class = "cv.DiscSurv")
  return(result)
}


cvf.DiscSurv <- function(i, data, Event.char, Z.char, Time.char, fold, cv.args){
  cv.args$data <- data[fold != i, , drop = FALSE]
  cv.args$Event.char <- Event.char
  cv.args$Z.char <- Z.char
  cv.args$Time.char <- Time.char
  
  
  fit.i <- do.call("DiscSurv", cv.args)  
  data.i <- data[fold == i, , drop = FALSE]
  yhat.i <- predict(fit.i, data.i, Event.char, Z.char, Time.char, 
                    lambda = fit.i$lambda, type = "response", return.Array = FALSE)
  
  data.small <- data.i[, c(Event.char, Time.char)]
  Y.i <- discSurv::dataLong(dataShort = data.small, timeColumn = Time.char, 
                            eventColumn = Event.char, timeAsFactor = TRUE)$y
  
  loss <- loss.Disc.Surv(Y.i, yhat.i)  
  list(loss = loss, nl = length(fit.i$lambda), yhat = yhat.i)
}

