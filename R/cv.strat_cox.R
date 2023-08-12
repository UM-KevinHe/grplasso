#' Cross-validation for penalized stratified cox model
#'
#' Performs k-fold cross validation for a penalized stratified cox model over a grid of values of regularization parameter lambda.
#'
#' @param data an `dataframe` or `list` object that contains the variables in the model.
#'
#' @param Event.char name of the event indicator in `data` as a character string. Event indicator should be a 
#' binary variable with 1 indicating that the event has occurred and 0 indicating (right) censoring.
#'
#' @param prov.char name of provider IDs variable in `data` as a character string. (can be seen as "stratum")
#'
#' @param Z.char names of covariates in `data` as vector of character strings.
#'
#' @param Time.char name of the follow up time in `data` as a character string.
#'
#' @param group a vector describing the grouping of the coefficients. If there are coefficients to be included in the model without being penalized, assign them to group 0 (or "0").
#'
#' @param nfolds the number of cross-validation folds. Default is 10.
#'
#' @param seed the seed of the random number generator in order to obtain reproducible results.
#'
#' @param fold a vector that specifies the fold that observations belongs to. By default the observations are randomly assigned.
#'
#' @param trace.cv \code{cv.strat_cox} will provide user with the progress of cross validation if `trace.cv = TRUE`. Default is FALSE.
#'
#' @param ... extra arguments to be passed to function.
#'
#'
#' @return An object with S3 class \code{cv.strat_cox}.
#'
#' \item{cve}{the error for each value of lambda, averaged across the cross-validation folds.}
#'
#' \item{cvse}{the estimated standard error associated with each value of for cve.}
#'
#' \item{lambda}{the sequence of regularization parameter values along which the cross-validation error was calculated.}
#'
#' \item{fit}{the fitted \code{strat_cox} object for the whole data.}
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
#' data(Cox_Data)
#' data <- Cox_Data$data
#' Event.char <- Cox_Data$Event.char
#' prov.char <- Cox_Data$prov.char
#' Z.char <- Cox_Data$Z.char
#' Time.char <- Cox_Data$Time.char
#' cv.fit <- cv.strat_cox(data, Event.char, prov.char, Z.char, Time.char, group = c(1, 2, 2, 3, 3), nfolds = 10, se = "quick")
#' # the best lambda using cross validation
#' cv.fit$lambda.min
#'
#' @references
#' K. He, J. Kalbfleisch, Y. Li, and et al. (2013) Evaluating hospital readmission rates in dialysis facilities; adjusting for hospital effects.
#' \emph{Lifetime Data Analysis}, \strong{19}: 490-512.
#' \cr

cv.strat_cox <- function(data, Event.char, prov.char, Z.char, Time.char, group = 1:length(Z.char), 
                         se=c('quick', 'bootstrap'), ...,  nfolds = 10, seed, fold, trace.cv = FALSE){
  if (missing(prov.char)){
    stop("stratum information must be provided!", call. = FALSE)
  }
  data <- data[order(data[, prov.char], data[, Time.char]), ]
  
  # "...": additional arguments to "grp.lasso"
  # "fold": a vector that specifies the fold that observations belongs to
  se <- match.arg(se)
  fit.args <- list(...)
  fit.args$data <- data
  fit.args$Event.char <- Event.char
  fit.args$prov.char <- prov.char
  fit.args$Z.char <- Z.char
  fit.args$Time.char <- Time.char
  fit.args$group <- group
  fit.args$returnX <- TRUE
  
  fit <- do.call("Strat.cox", fit.args)  #fit the entire model
  
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
  
  # we need each training data (9/10 data) contains all strata?
  # original.count.stratum <- length(unique(data[, prov.char]))

  if (missing(fold)){ # make sure each fold contains same proportion of Y = 0 and Y = 1
    ind1 <- which(delta.obs == 1)
    ind0 <- which(delta.obs == 0)
    n1 <- length(ind1)
    n0 <- length(ind0)
    fold1 <- 1:n1 %% nfolds  # assign observations to different folds
    fold0 <- (n1 + 1:n0) %% nfolds
    fold1[fold1==0] <- nfolds # set fold "0" to fold max
    fold0[fold0==0] <- nfolds
    fold <- integer(n.obs)
    fold[delta.obs == 1] <- sample(fold1)
    fold[delta.obs == 0] <- sample(fold0)
  } else {
    nfolds <- max(fold)
  }
  
  # Cross-Validation
  Y <- matrix(NA, nrow = n.obs, ncol = length(fit$lambda))
  cv.args <- list(...)
  cv.args$lambda <- fit$lambda
  cv.args$group <- ZG$g
  cv.args$group.multiplier <- ZG$m
  
  for (i in 1:nfolds) {
    if (trace.cv == TRUE){
      cat("Starting CV fold #", i, sep="","...\n")
    }
    res <- cvf.strat_cox(i, data, Event.char, Z.char, prov.char, Time.char, fold, cv.args)
    Y[fold == i, 1:res$nl] <- res$yhat # predicted "eta"
  }
  
  # Eliminate saturated lambda values, if any
  ind <- which(apply(is.finite(Y), 2, all))
  Y <- Y[, ind] # will keep the increasing order of "time" within each provider
  lambda <- fit$lambda[ind]
  
  #define stratum information
  unique.prov <- unique(data[, prov.char])
  prov.ref <- cbind(1:length(unique.prov), unique.prov)
  colnames(prov.ref) <- c("New.ID", prov.char)
  ID <- as.matrix(data[, prov.char])  #stratum indicator
  colnames(ID) <- prov.char
  ID <- merge(ID, prov.ref, by = prov.char)[, 2, drop = F] 
  colnames(ID) <- prov.char
  
  if (se == "quick") {
    L <- loss.strat_cox(delta.obs, Y, ID, total=FALSE)
    cve <- apply(L, 2, sum)/sum(delta.obs)
    cvse <- apply(L, 2, sd)*sqrt(nrow(L))/sum(delta.obs)
  } else {
    cve <- as.double(loss.strat_cox(delta.obs, Y, ID, total = TRUE))/sum(delta.obs)
    cvse <- se.strat_cox(delta.obs, Y, ID)/sum(delta.obs)
  }

  min <- which.min(cve)  #find index of lambda with minimum cve
  
  result <- structure(list(cve = cve,
                           cvse = cvse,
                           lambda = lambda,
                           fit = fit,
                           fold = fold,
                           min = min,
                           lambda.min = lambda[min]),
                      class = "cv.strat_cox")
  return(result)
}

cvf.strat_cox <- function(i, data, Event.char, Z.char, prov.char, Time.char, fold, cv.args){
  cv.args$data <- data[fold != i, , drop = FALSE]  #will not change the (increasing) order of time within each provider
  cv.args$Event.char <- Event.char
  cv.args$prov.char <- prov.char
  cv.args$Z.char <- Z.char
  cv.args$Time.char <- Time.char
  
  fit.i <- do.call("Strat.cox", cv.args)  #fit the discrete survival model using one training data set (9/10 data)
  data.i <- data[fold == i, , drop = FALSE]  #current validation data
  yhat.i <- predict(fit.i, data.i, Event.char, prov.char, Z.char, Time.char, 
                    lambda = fit.i$lambda, type = "link") # exp(eta) over all given lambda
  list(nl = length(fit.i$lambda), 
       yhat = yhat.i)
}

se.strat_cox <-function(delta.obs, y.hat, ID, B = 100) {
  cve <- matrix(NA, B, ncol(y.hat))
  for (b in 1:B) {
    ind <- sort(sample(1:nrow(y.hat), replace=TRUE)) #"sort" to keep the time order
    cve[b,] <- loss.strat_cox(delta.obs[ind, , drop = F], y.hat[ind, , drop = F], ID[ind, , drop = F])
  }
  return(apply(cve, 2, sd))
}
