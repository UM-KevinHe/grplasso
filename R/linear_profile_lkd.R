#' Fit the linear model using profile likelihood
#'
#' Fitting a linear model for data with high-dimensional center effect through profile likelihood estimation
#'
#' @param data an `dataframe` or `list` object that contains the variables in the model.
#'
#' @param Y.char name of the response variable in `data` as a character string.
#'
#' @param Z.char names of covariates in `data` as vector of character strings.
#'
#' @param prov.char name of provider IDs variable in `data` as a character string. If "prov.char" is not specified, all observations are are considered to be from the same provider.
#'
#' @param ... extra arguments to be passed to function.
#'
#' @return An object with S3 class \code{proflkd.linear}
#'
#' \item{beta}{the fitted matrix of covariate coefficients.}
#' 
#' \item{gamma}{the fitted matrix of provider effects.}
#'
#' @export
#'
#' @examples
#' data(linear_data)
#' data <- linear_data$data
#' Y.char <- linear_data$Y.char
#' prov.char <- linear_data$prov.char
#' Z.char <- linear_data$Z.char
#' fit <- prof_lkd.linear(data, Y.char, Z.char, prov.char)
#' fit$beta[1:5, ]
#' fit$gamma[1:5, ]
#'
#' @importFrom Rcpp evalCpp
#'

prof_lkd.linear <- function(data, Y.char, Z.char, prov.char, ...){
  data <- data[order(factor(data[, prov.char])), ]
  ID <- as.matrix(data[, prov.char])
  Y <- data[, Y.char]
  n.prov <- sapply(split(Y, ID), length)
  m <- length(n.prov)
  sum.first.term <- matrix(0, nrow = length(Z.char), ncol = length(Z.char))
  sum.second.term <- matrix(0, nrow = length(Z.char), ncol = 1)
  for (j in names(n.prov)){
    temp.X <- as.matrix(data[which(data$prov.ID == j), Z.char])
    temp.Y <- as.matrix(data[which(data$prov.ID == j), Y.char, drop = F])
    n <- length(temp.Y)
    Qn <- diag(1, nrow = n, ncol = n) - matrix(1, nrow = n, ncol = n)/n
    sum.first.term <- sum.first.term + t(temp.X) %*% Qn %*% temp.X
    sum.second.term <- sum.second.term + t(temp.X) %*% Qn %*% temp.Y
  }
  beta <- solve(sum.first.term) %*% sum.second.term
  colnames(beta) <- "coef"
  
  gamma <- c()
  for (j in names(n.prov)){
    temp.y_bar <- mean(as.matrix(data[which(data$prov.ID == j), Y.char, drop = F]))
    temp.x_bar <- as.matrix(colMeans(as.matrix(data[which(data$prov.ID == j), Z.char])))
    gamma <- c(gamma, temp.y_bar - t(temp.x_bar) %*% beta)
  }
  gamma <- as.matrix(gamma)
  colnames(gamma) <- "coef"
  rownames(gamma)<- names(n.prov)
  
  result <- structure(list(beta = beta,
                           gamma = gamma),
                      class = "proflkd.linear")  #define a list for prediction
  return(result)
}

