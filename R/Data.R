#' example data for pp.lasso and grp.lasso
#'
#' A simulated data set containing response variable, provider information and 5 covariates.
#' @name GLM_Data
#' @docType data
#' @usage data(GLM_Data)
#'
#' @format A list containing the following elements:
#' \describe{
#'   \item{data}{example data.`Y` is the response variable; `Prov.ID` is the provider indicator; `Z1`, ..., `Z5` are 5 continuous covariates.}
#'   \item{Y.char}{variable name of the response variable.}
#'   \item{prov.char}{variable name of the provider indicator.}
#'   \item{Z.char}{variable names of covariates.}
#'   \item{group}{vector describing how the covariates are grouped.}
#'
#' }
"GLM_Data"

#' example data for pp.DiscSurv
#'
#' A simulated data set containing observation time, event indicator, provider information and 5 covariates.
#' @name Surv_Data
#' @docType data
#' @usage data(Surv_Data)
#'
#' @format A list containing the following elements:
#' \describe{
#'   \item{data}{example data.`time` represents the observation time; `status` is the event indicator; `Prov.ID` is the provider indicator; `Z1`, ..., `Z5` are 5 continuous covariates.}
#'   \item{Event.char}{variable name of the event indicator.}
#'   \item{prov.char}{variable name of the provider indicator.}
#'   \item{Z.char}{variable names of covariates.}
#'   \item{Time.char}{variable name of the observation time.}
#' }
"Surv_Data"
