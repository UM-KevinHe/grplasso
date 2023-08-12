#' Example data for generalize linear model
#'
#' A simulated data set containing response variable, provider information and 5 covariates.
#' @name GLM_Data
#' @docType data
#' @usage data(GLM_Data)
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
#' data(GLM_Data)
#' data <- GLM_Data$data
#' head(data)
#' 
"GLM_Data"

#' Example data for discrete survival model
#'
#' A simulated data set containing observation time, event indicator, provider information and 5 covariates.
#' @name Surv_Data
#' @docType data
#' @usage data(Surv_Data)
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
#' data(Surv_Data)
#' data <- Surv_Data$data
#' head(data)
"Surv_Data"

#' Example data for Cox model
#'
#' A simulated data set containing observation time, event indicator, provider information and 5 covariates.
#' @name Cox_Data
#' @docType data
#' @usage data(Cox_Data)
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
#' data(Cox_Data)
#' data <- Cox_Data$data
#' head(data)
"Cox_Data"
