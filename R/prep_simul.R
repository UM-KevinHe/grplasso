# prepare simulation data
#' Screen and prepare clustered data before fitting
#'
#' Checks a data frame for missingness, providers below a size cutoff,
#' near-zero-variance covariates, pairwise correlation and multicollinearity,
#' and returns the data with an `included` indicator.
#'
#' @param data a data frame containing the outcome, the covariates and the
#'   provider identifier.
#' @param Y.char name of the outcome column.
#' @param Z.char names of the covariate columns.
#' @param prov.char name of the provider identifier column.
#' @param cutoff providers with fewer than `cutoff` observations are flagged.
#' @param check whether to run the screening checks. Default `TRUE`.
#'
#' @return The input data with an added `included` column.
#'
#' @section Optional dependencies:
#' The near-zero-variance and VIF checks use \pkg{caret} and \pkg{olsrr},
#' which are in `Suggests` rather than `Imports`. If either is not installed the
#' corresponding check is skipped with a message rather than failing.
#'
#' @examples
#' set.seed(5)
#' sim <- Simulation_data_GroupLasso(
#'   list(m = 10, n.beta = 6, n.groups = 3, prop.NonZero.group = 0.5,
#'        prop.outlier = 0, rho = 0.5),
#'   prov.size.mean = 30)
#' d <- sim$sim.data
#' prepped <- fe.data.prep(d, "Y", colnames(d)[3:8], "Prov.ID",
#'                         cutoff = 5, check = FALSE)
#' table(prepped$included)
#'
#' @export
fe.data.prep <- function(data, Y.char, Z.char, prov.char, cutoff=10, check=TRUE) {
  if (check) {
    message("Checking absence of variables ... ")
    Y.ind <- match(Y.char, names(data))
    if (is.na(Y.ind)) {
      stop(paste("Response variable '", Y.char, "' NOT found!", sep=""),call.=FALSE)
    }
    Z.ind <- match(Z.char, names(data))
    if (sum(is.na(Z.ind)) > 0) {
      stop(paste("Covariate(s) '", paste(Z.char[is.na(Z.ind)], collapse="', '"), "' NOT found!", sep=""),call.=FALSE)
    }
    prov.ind <- match(prov.char, names(data))
    if (is.na(prov.ind)) {
      stop(paste("Provider ID '", prov.char, "' NOT found!", sep=""),call.=FALSE)
    }
    message("Checking absence of variables completed!")
    message("Checking missingness of variables ... ")
    if (sum(complete.cases(data[,c(Y.char,Z.char,prov.char)]))==NROW(data)) {
      message("Missing values NOT found. Checking missingness of variables completed!")
    } else {
      check.na <- function(name) {
        if (sum(is.na(data[,name])) > 0) {
          warning(sum(is.na(data[,name]))," out of ",NROW(data[,name])," in '",name,"' missing!",immediate.=TRUE,call.=FALSE)
        }
      }
      invisible(sapply(c(Y.char,Z.char,prov.char), check.na))
      missingness <- (1 - sum(complete.cases(data[,c(Y.char,Z.char,prov.char)])) / NROW(data)) * 100
      stop(paste(round(missingness,2), "% of all observations are missing!",sep=""),call.=FALSE)
    }
    ## caret and olsrr are in Suggests, not Imports: they are only needed for
    ## these optional screening checks, and caret alone pulls in dozens of
    ## packages.  Skip with a message rather than fail if they are absent.
    if (!requireNamespace("caret", quietly = TRUE)) {
      message("Package 'caret' not installed; skipping the near-zero-variance check.")
      nzv <- NULL
    } else {
    message("Checking variation in covariates ... ")
    nzv <- caret::nearZeroVar(data[,Z.char], saveMetrics=TRUE)
    if (sum(nzv$zeroVar==TRUE) > 0) {
      stop("Covariate(s) '", paste(row.names(nzv[nzv$zeroVar==TRUE,]), collapse="', '"),
           "' with zero variance(s)!", call.=FALSE)
    } else if (sum(nzv$nzv==TRUE) > 0) {
      warning("Covariate(s) '",paste(row.names(nzv[nzv$nzv==TRUE,]), collapse="', '"),
              "' with near zero variance(s)!",immediate.=TRUE,call.=FALSE)
    }
    message("Checking variation in covariates completed!")
    }
    message("Checking pairwise correlation among covariates ... ")
    cor <- cor(data[,Z.char])
    threshold.cor <- 0.9
    if (sum(abs(cor[upper.tri(cor)])>threshold.cor) > 0) {
      cor[lower.tri(cor,diag=TRUE)] <- 0
      ind <- which(abs(cor)>threshold.cor)
      pairs <- sapply(ind, function(ind) c(rownames(cor)[ind%%NROW(cor)],
                                           colnames(cor)[ind%/%NROW(cor)+1]))
      warning("The following ", NCOL(pairs),
              " pair(s) of covariates are highly correlated (correlation > ",
              threshold.cor,"): ", immediate.=TRUE, call.=FALSE)
      invisible(apply(pairs,2,function(col) message('("',paste(col, collapse='", "'),'")')))
    }
    message("Checking pairwise correlation among covariates completed!")
    ## check VIF
    if (!requireNamespace("olsrr", quietly = TRUE)) {
      message("Package 'olsrr' not installed; skipping the VIF check.")
    } else {
    message("Checking VIF of covariates ... ")
    m.lm <- lm(as.formula(paste(Y.char,"~",paste(Z.char, collapse="+"))), data=data)
    vif <- olsrr::ols_vif_tol(m.lm)
    if(sum(vif$VIF >= 10) > 0){
      warning("Covariate(s) '",
              paste(as.data.frame(vif)[vif$VIF>=10,"Variables"], collapse="', '"),
              "' with serious multicollinearity!",immediate.=TRUE,call.=FALSE)
    }
    message("Checking VIF of covariates completed!")
    }
  }
  data <- data[order(factor(data[,prov.char])),]
  prov.size <- as.integer(table(data[,prov.char]))
  prov.size.long <- rep(prov.size,prov.size)
  data$included <- 1 * (prov.size.long > cutoff)
  warning(sum(prov.size<=cutoff)," out of ",length(prov.size),
          " providers considered small and filtered out!",immediate.=TRUE,call.=FALSE)
  prov.list <- unique(data[data$included==1,prov.char])
  prov.no.readm <-
    prov.list[sapply(split(data[data$included==1,Y.char], factor(data[data$included==1,prov.char])),sum)==0]
  data$no.readm <- 0
  data$no.readm[data[,prov.char] %in% c(prov.no.readm)] <- 1
  message(paste(length(prov.no.readm),"out of",length(prov.list),
                "remaining providers with no readmission within 30 days."))
  prov.all.readm <-
    prov.list[sapply(split(1-data[data$included==1,Y.char],factor(data[data$included==1,prov.char])),sum)==0]
  data$all.readm <- 0
  data$all.readm[data[,prov.char]%in%c(prov.all.readm)] <- 1
  message(paste(length(prov.all.readm),"out of",length(prov.list),
                "remaining providers with all readmissions within 30 days."))
  message(paste0("After screening, ", round(sum(data[data$included==1,Y.char])/length(data[data$included==1,Y.char])*100,2),
                 "% of all discharges were readmitted within 30 days."))
  return(data)
}

