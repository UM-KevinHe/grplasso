#' Example data for generalized linear models
#'
#' A simulated data set containing response variable, provider information and 5 covariates.
#' @name BinaryData
#' @docType data
#' @usage data(BinaryData)
#' @keywords datasets
#'
#' @format A list containing the following elements:
#' \describe{
#'   \item{data}{example data. \code{Y} is the response variable; \code{Prov.ID} is the center indicator (include 20 centers); \code{Z1}, ..., \code{Z5} are 5 continuous covariates.}
#'   \item{Y.char}{variable name of the response variable.}
#'   \item{prov.char}{variable name of the provider indicator.}
#'   \item{Z.char}{variable names of covariates.}
#'   \item{group}{vector describing how the covariates are grouped.}
#'
#' }
#' 
#' @examples
#' data(BinaryData)
#' data <- BinaryData$data
#' head(data)
#' 
"BinaryData"

#' Example data for discrete survival models
#'
#' A simulated data set containing observation time, event indicator, provider information and 5 covariates.
#' @name DiscTime
#' @docType data
#' @usage data(DiscTime)
#' @keywords datasets
#' 
#' @format A list containing the following elements:
#' \describe{
#'   \item{data}{example data.\code{time} represents the observation time; \code{status} is the event indicator; \code{Prov.ID} is the center indicator (include 5 centers); \code{Z1}, ..., \code{Z5} are 5 continuous covariates.}
#'   \item{Event.char}{variable name of the event indicator.}
#'   \item{prov.char}{variable name of the provider indicator.}
#'   \item{Z.char}{variable names of covariates.}
#'   \item{Time.char}{variable name of the observation time.}
#' }
#' 
#' @examples
#' data(DiscTime)
#' data <- DiscTime$data
#' head(data)
"DiscTime"

#' Example data for Cox models
#'
#' A simulated data set containing observation time, event indicator, provider information and 5 covariates.
#' @name ContTime
#' @docType data
#' @usage data(ContTime)
#' @keywords datasets
#' 
#' @format A list containing the following elements:
#' \describe{
#'   \item{data}{example data.\code{time} represents the observation time; \code{status} is the event indicator; \code{Prov.ID} is the center indicator (include 5 centers); \code{Z1}, ..., \code{Z5} are 5 continuous covariates.}
#'   \item{Event.char}{variable name of the event indicator.}
#'   \item{prov.char}{variable name of the provider indicator.}
#'   \item{Z.char}{variable names of covariates.}
#'   \item{Time.char}{variable name of the observation time.}
#' }
#' 
#' @examples
#' data(ContTime)
#' data <- ContTime$data
#' head(data)
"ContTime"


