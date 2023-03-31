#' example data for `pp.lasso` and `grp.lasso`
#'
#' A simulated data set containing response variable, provider information and 5 covariates.
#' @name GLM_Data
#' @docType data
#' @usage data(GLM_Data)
#'
#' @format A list containing the following elements:
#' \describe{
#'   \item{Y}{response variable.}
#'   \item{Prov.ID}{provider information.}
#'   \item{Z1, ..., Z5}{5 continuous covariates.}
#' }
"GLM_Data"

#' example data for `pp.DiscSurv`
#'
#' A simulated data set containing time, status, provider information and 5 covariates.
#' @name Surv_Data
#' @docType data
#' @usage data(Surv_Data)
#'
#' @format A list containing the following elements:
#' \describe{
#'   \item{time}{observation time.}
#'   \item{status}{event indicator.}
#'   \item{Prov.ID}{provider information.}
#'   \item{Z1, ..., Z5}{5 continuous covariates.}
#' }
"Surv_Data"
