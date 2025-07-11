% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.cve.ppLasso.R
\name{plot.cv.ppLasso}
\alias{plot.cv.ppLasso}
\alias{plot.cv.gr_ppLasso}
\title{Plot the cross entropy loss from a cv.ppLasso or cv.gr_ppLasso object}
\usage{
\method{plot}{cv.ppLasso}(
  fit,
  log.x = T,
  vertical.line = T,
  col.vertical.line = "blue",
  col.dot = "red"
)

\method{plot}{cv.gr_ppLasso}(
  fit,
  log.x = T,
  vertical.line = T,
  col.vertical.line = "blue",
  col.dot = "red"
)
}
\arguments{
\item{fit}{a \code{cv.gr_ppLasso} object.}

\item{log.x}{whether the horizontal axis be on the log scale.}

\item{vertical.line}{whether draws a vertical line at the value where cross-validaton error is minimized.}
}
\description{
Return the plot of the cross entropy loss from a \code{cv.ppLasso} or \code{cv.gr_ppLasso} object
}
\examples{
data(BinaryData)
data <- BinaryData$data
Y.char <- BinaryData$Y.char
prov.char <- BinaryData$prov.char
Z.char <- BinaryData$Z.char
cv.fit.pplasso <- cv.pp.lasso(data, Y.char, Z.char, prov.char, nfolds = 10)
plot(cv.fit.pplasso)
data(BinaryData)
data <- BinaryData$data
Y.char <- BinaryData$Y.char
prov.char <- BinaryData$prov.char
Z.char <- BinaryData$Z.char
group <- BinaryData$group
cv.fit.grplasso <- cv.grp.lasso(data, Y.char, Z.char, prov.char, group = group, nfolds = 10)
plot(cv.fit.grplasso)

}
