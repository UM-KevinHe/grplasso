#' cross-validation for pp.lasso
#'
#' performs k-fold cross validation for penalized regression models over a grid of values of regularization parameter lambda.
#'
#' @param data an `dataframe` or `list` object that contains the variables in the model.
#'
#' @param Y.char name of the response variable from `data` as a character string, as in \code{pp.lasso} function.
#'
#' @param Z.char names of covariates from `data` as vector of character strings, as in \code{pp.lasso} function.
#'
#' @param prov.char name of provider IDs variable from `data` as a character string, as in \code{pp.lasso} function.
#'
#' @param penalize.x  a vector indicates whether the corresponding covariate will be penalized, as in \code{pp.lasso} function.
#'
#' @param nfolds the number of cross-validation folds. Default is 10.
#'
#' @param seed the seed of the random number generator in order to obtain reproducible results.
#'
#' @param fold a vector that specifies the fold that observations belongs to. By default the observations are randomly assigned.
#'
#' @param trace.cv \code{cv.pp.lasso} will provide user with the progress of cross validation if `trace.cv = TRUE`. Default is FALSE.
#'
#' @param ... extra arguments to be passed to function.
#'
#' @return An object with S3 class \code{cv.ppLasso}.
#'
#' \item{cve}{the error for each value of lambda, averaged across the cross-validation folds.}
#'
#' \item{cvse}{the estimated standard error associated with each value of for cve.}
#'
#' \item{lambda}{the sequence of regularization parameter values along which the cross-validation error was calculated.}
#'
#' \item{fit}{the fitted \code{pp.lasso} object for the whole data.}
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
#' data(GLM_Data)
#' data <- GLM_Data$data
#' Y.char <- GLM_Data$Y.char
#' prov.char <- GLM_Data$prov.char
#' Z.char <- GLM_Data$Z.char
#' cv.fit <- cv.pp.lasso(data, Y.char, Z.char, prov.char, nfolds = 10)
#' # the best lambda using cross validation
#' cv.fit$lambda.min
#'
#' @references
#' K. He, J. Kalbfleisch, Y. Li, and et al. (2013) Evaluating hospital readmission rates in dialysis facilities; adjusting for hospital effects.
#' \emph{Lifetime Data Analysis}, \strong{19}: 490-512.
#' \cr

cv.pp.lasso <- function(data, Y.char, Z.char, prov.char, penalize.x = rep(1, length(Z.char)), ..., nfolds = 10,
                        seed, fold, trace.cv = FALSE){
  if (missing(prov.char)){ #single intercept
    prov.char <- "intercept"
    data$intercept <- matrix(1, nrow = nrow(data))
  }

  # "...": additional arguments to "pp.lasso"
  # "fold": a vector that specifies the fold that observations belongs to
  fit.args <- list(...)
  fit.args$data <- data
  fit.args$Y.char <- Y.char
  fit.args$Z.char <- Z.char
  fit.args$prov.char <- prov.char
  fit.args$penalize.x <- penalize.x
  fit.args$returnX <- TRUE

  fit <- do.call("pp.lasso", fit.args)  #fit the entire model

  # get standardized Z,
  ZG <- fit$returnX
  Z <- ZG$std.Z
  y <- as.matrix(data[, Y.char])
  returnX <- list(...)$returnX
  if (is.null(returnX) || !returnX){
    fit$returnX <- NULL
  }

  #setup folds
  if (!missing(seed)) {
    set.seed(seed)
  }

  n.obs <- length(y)

  if (missing(fold)){ # make sure each fold contains same proportion of Y = 0 and Y = 1
    ind1 <- which(y == 1)
    ind0 <- which(y == 0)
    n1 <- length(ind1)
    n0 <- length(ind0)
    fold1 <- 1:n1 %% nfolds  # assign observations to different folds
    fold0 <- (n1 + 1:n0) %% nfolds
    fold1[fold1==0] <- nfolds # set fold "0" to fold max
    fold0[fold0==0] <- nfolds
    fold <- integer(n.obs)
    fold[y == 1] <- sample(fold1)
    fold[y == 0] <- sample(fold0)
  } else {
    nfolds <- max(fold)
  }

  # Do Cross-Validation
  E <- Y <- matrix(NA, nrow = length(y), ncol = length(fit$lambda))
  PE <- E
  cv.args <- list(...)
  cv.args$lambda <- fit$lambda
  cv.args$group <- ZG$g
  cv.args$penalized.multiplier <- ZG$m

  for (i in 1:nfolds) {
    if (trace.cv == TRUE){
      cat("Starting CV fold #", i, sep="","...\n")
    }

    res <- cvf.pplasso(i, data, Y.char, Z.char, prov.char, fold, cv.args)
    Y[fold == i, 1:res$nl] <- res$yhat
    E[fold == i, 1:res$nl] <- res$loss
    PE[fold == i, 1:res$nl] <- res$pe  # wrong predict class
  }

  # Eliminate saturated lambda values, if any
  ind <- which(apply(is.finite(E), 2, all))
  E <- E[, ind, drop=FALSE]
  Y <- Y[, ind]
  lambda <- fit$lambda[ind]

  cve <- apply(E, 2, mean)  # mean cross entropy loss within each lambda value
  cvse <- apply(E, 2, sd) / sqrt(n.obs) # standardized loss
  min <- which.min(cve)  #find index of lambda with minimum cve
  pe <- apply(PE[, ind], 2, mean)  # mean predict class error

  result <- structure(list(cve = cve,
                           cvse = cvse,
                           lambda = lambda,
                           fit = fit, #model with entire data
                           pe = pe,
                           fold = fold,
                           min = min,
                           lambda.min = lambda[min]),
                      class = "cv.ppLasso")
  return(result)
}
