% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stratified_cox.R
\name{Strat.cox}
\alias{Strat.cox}
\title{Fit a penalized stratified cox model}
\usage{
Strat.cox(
  data,
  Event.char,
  Z.char,
  Time.char,
  prov.char,
  group = 1:length(Z.char),
  group.multiplier,
  standardize = T,
  lambda,
  nlambda = 100,
  lambda.min.ratio = 0.001,
  lambda.early.stop = FALSE,
  nvar.max = p,
  group.max = length(unique(group)),
  stop.loss.ratio = 0.001,
  tol = 1e-04,
  max.each.iter = 10000,
  max.total.iter = (max.each.iter * nlambda),
  actSet = TRUE,
  actIter = max.each.iter,
  actGroupNum = sum(unique(group) != 0),
  actSetRemove = F,
  returnX = FALSE,
  trace.lambda = FALSE,
  ...
)
}
\arguments{
\item{data}{an \code{dataframe} or \code{list} object that contains the variables in the model.}

\item{Event.char}{name of the event indicator in \code{data} as a character string. Event indicator should be a
binary variable with 1 indicating that the event has occurred and 0 indicating (right) censoring.}

\item{Z.char}{names of covariates in \code{data} as vector of character strings.}

\item{Time.char}{name of the follow up time in \code{data} as a character string.}

\item{prov.char}{name of stratum indicator in \code{data} as a character string.
If "prov.char" is not specified, all observations are are considered to be from the same stratum.}

\item{group}{a vector describing the grouping of the coefficients. If there are coefficients to be included in the model without being penalized, assign them to group 0 (or "0").}

\item{group.multiplier}{A vector of values representing multiplicative factors by which each covariate's penalty is to be multiplied. Default is a vector of 1's.}

\item{standardize}{logical flag for x variable standardization, prior to fitting the model sequence. The coefficients are always returned on the original scale. Default is \code{standardize=TRUE}.}

\item{lambda}{a user supplied lambda sequence. Typical usage is to have the program compute its own lambda sequence based on \code{nlambda} and \code{lambda.min.ratio}.}

\item{nlambda}{the number of lambda values. Default is 100.}

\item{lambda.min.ratio}{the fraction of the smallest value for lambda with \code{lambda.max} (smallest lambda for which all coefficients are zero) on log scale. Default is 1e-03.}

\item{lambda.early.stop}{whether the program stop before running the entire sequence of lambda. Early stop based on the ratio of deviance for models under two successive lambda. Default is \code{FALSE}.}

\item{nvar.max}{number of maximum selected variables. Default is the number of all covariates.}

\item{group.max}{number of maximum selected groups. Default is the number of all groups.}

\item{stop.loss.ratio}{if \code{lambda.early.stop = TRUE}, the ratio of loss for early stopping. Default is 1e-3.}

\item{tol}{convergence threshold. For each lambda, the program will stop if the maximum change of covariate coefficient is smaller than \code{tol}. Default is 1e-4.}

\item{max.each.iter}{maximum number of iterations for each lambda. Default is 1e4.}

\item{max.total.iter}{maximum number of iterations for entire path. Default is \code{max.each.iter} * \code{nlambda}.}

\item{actSet}{whether to use the active method for variable selection. Default is TRUE.}

\item{actIter}{if \code{actSet = TRUE}, the maximum number of iterations for a new updated active set. Default is \code{max.each.iter} (i.e. we will update the current active set until convergence ).}

\item{actGroupNum}{if \code{actSet = TRUE}, the maximum number of variables that can be selected into the new active set for each time when the active set is updated. Default is number of groups.}

\item{actSetRemove}{if \code{actSet = TRUE}, whether we remove the zero coefficients from the current active set. Default is FALSE.}

\item{returnX}{whether return the standardized design matrix. Default is FALSE.}

\item{trace.lambda}{whether display the progress for fitting the entire path. Default is FALSE.}

\item{...}{extra arguments to be passed to function.}
}
\value{
An object with S3 class \code{strat_cox}.

\item{beta}{the fitted matrix of covariate coefficients.
The number of rows is equal to the number of coefficients,
and the number of columns is equal to nlambda.}

\item{group}{a vector describing the grouping of the coefficients.}

\item{lambda}{the sequence of \code{lambda} values in the path.}

\item{loss}{the likelihood of the fitted model at each value of \code{lambda}.}

\item{linear.predictors}{the linear predictors of the fitted model at each value of \code{lambda}.}

\item{df}{the estimates of effective number of selected variables all the points along the regularization path.}

\item{iter}{the number of iterations until convergence at each value of \code{lambda}.}
}
\description{
Main function for fitting a penalized stratified cox model.
}
\details{
The model is fit by Newton method and coordinate descent method.
}
\examples{
data(ContTime)
data <- ContTime$data
Event.char <- ContTime$Event.char
prov.char <- ContTime$prov.char
Z.char <- ContTime$Z.char
Time.char <- ContTime$Time.char
fit <- Strat.cox(data, Event.char, Z.char, Time.char, prov.char, group = c(1, 2, 2, 3, 3))
fit$beta[, 1:5]

}
\references{
K. He, J. Kalbfleisch, Y. Li, and et al. (2013) Evaluating hospital readmission rates in dialysis facilities; adjusting for hospital effects.
\emph{Lifetime Data Analysis}, \strong{19}: 490-512.
\cr
}
\seealso{
\code{\link{coef}}, \code{\link{plot}} function.
}
