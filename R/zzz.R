# Internal helpers are deliberately NOT exported.  This package used to carry a
# blanket exportPattern of "^[[:alpha:]]+", which exported every internal
# function and left 44 objects exported-but-undocumented in R CMD check.  The
# user-facing API now carries an explicit @export on each function.

#' @useDynLib grplasso, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats approx as.formula complete.cases lm plogis predict rbinom
#'   relevel rnorm rpois runif sd
#' @importFrom utils tail
NULL

# `Individual` is a column name used inside tidyr/dplyr pipelines, so R CMD
# check sees it as an undefined global.  Declare it.
utils::globalVariables("Individual")
