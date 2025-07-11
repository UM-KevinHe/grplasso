% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Data.R
\docType{data}
\name{ContTime}
\alias{ContTime}
\title{Example data for Cox models}
\format{
A list containing the following elements:
\describe{
\item{data}{example data.\code{time} represents the observation time; \code{status} is the event indicator; \code{Prov.ID} is the center indicator (include 5 centers); \code{Z1}, ..., \code{Z5} are 5 continuous covariates.}
\item{Event.char}{variable name of the event indicator.}
\item{prov.char}{variable name of the provider indicator.}
\item{Z.char}{variable names of covariates.}
\item{Time.char}{variable name of the observation time.}
}
}
\usage{
data(ContTime)
}
\description{
A simulated data set containing observation time, event indicator, provider information and 5 covariates.
}
\examples{
data(ContTime)
data <- ContTime$data
head(data)
}
\keyword{datasets}
