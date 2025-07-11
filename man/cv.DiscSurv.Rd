% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cv.DiscSurv.R
\name{cv.DiscSurv}
\alias{cv.DiscSurv}
\title{Cross-validation for penalized discrete survival model (without provider information)}
\usage{
cv.DiscSurv(
  data,
  Event.char,
  Z.char,
  Time.char,
  penalize.x = rep(1, length(Z.char)),
  ...,
  nfolds = 10,
  seed,
  fold,
  trace.cv = FALSE
)
}
\arguments{
\item{data}{an \code{dataframe} or \code{list} object that contains the variables in the model.}

\item{Event.char}{name of the event indicator in \code{data} as a character string. Event indicator should be a
binary variable with 1 indicating that the event has occurred and 0 indicating (right) censoring.}

\item{Z.char}{names of covariates in \code{data} as vector of character strings.}

\item{Time.char}{name of the follow up time in \code{data} as a character string.}

\item{penalize.x}{a vector indicates whether the corresponding covariate will be penalized, as in \code{pp.DiscSurv} function.}

\item{...}{extra arguments to be passed to function.}

\item{nfolds}{the number of cross-validation folds. Default is 10.}

\item{seed}{the seed of the random number generator in order to obtain reproducible results.}

\item{fold}{a vector that specifies the fold that observations belongs to. By default the observations are randomly assigned.}

\item{trace.cv}{\code{cv.DiscSurv} will provide user with the progress of cross validation if \code{trace.cv = TRUE}. Default is FALSE.}
}
\value{
An object with S3 class \code{cv.DiscSurv}.

\item{cve}{the error for each value of lambda, averaged across the cross-validation folds.}

\item{cvse}{the estimated standard error associated with each value of for cve.}

\item{lambda}{the sequence of regularization parameter values along which the cross-validation error was calculated.}

\item{fit}{the fitted \code{pp.DiscSurv} object for the whole data.}

\item{fold}{the fold assignments for cross-validation for each observation}

\item{min}{the index of lambda corresponding to lambda.min.}

\item{lambda.min}{the value of lambda with the minimum cross-validation error.}
}
\description{
Performs k-fold cross validation for penalized discrete survival model (without provider information) over a grid of values of regularization parameter lambda.
}
\examples{
data(DiscTime)
data <- DiscTime$data
Event.char <- DiscTime$Event.char
Z.char <- DiscTime$Z.char
Time.char <- DiscTime$Time.char
cv.fit <- cv.DiscSurv(data, Event.char, Z.char, Time.char, nfolds = 10, trace.cv = T)
cv.fit$cve
cv.fit$lambda.min

}
\references{
K. He, J. Kalbfleisch, Y. Li, and et al. (2013) Evaluating hospital readmission rates in dialysis facilities; adjusting for hospital effects.
\emph{Lifetime Data Analysis}, \strong{19}: 490-512.
\cr
}
