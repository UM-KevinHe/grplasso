# prepare simulation data
fe.data.prep <- function(data, Y.char, Z.char, prov.char, cutoff=10, check=TRUE) {
  if (check) {
    message("Checking absence of variables ... ")
    Y.ind <- match(Y.char, names(data))
    if (is.na(Y.ind)) {
      stop(paste("Response variable '", Y.char, "' NOT found!", sep=""),call.=F)
    }
    Z.ind <- match(Z.char, names(data))
    if (sum(is.na(Z.ind)) > 0) {
      stop(paste("Covariate(s) '", paste(Z.char[is.na(Z.ind)], collapse="', '"), "' NOT found!", sep=""),call.=F)
    }
    prov.ind <- match(prov.char, names(data))
    if (is.na(prov.ind)) {
      stop(paste("Provider ID '", prov.char, "' NOT found!", sep=""),call.=F)
    }
    message("Checking absence of variables completed!")
    message("Checking missingness of variables ... ")
    if (sum(complete.cases(data[,c(Y.char,Z.char,prov.char)]))==NROW(data)) {
      message("Missing values NOT found. Checking missingness of variables completed!")
    } else {
      check.na <- function(name) {
        if (sum(is.na(data[,name])) > 0) {
          warning(sum(is.na(data[,name]))," out of ",NROW(data[,name])," in '",name,"' missing!",immediate.=T,call.=F)
        }
      }
      invisible(sapply(c(Y.char,Z.char,prov.char), check.na))
      missingness <- (1 - sum(complete.cases(data[,c(Y.char,Z.char,prov.char)])) / NROW(data)) * 100
      stop(paste(round(missingness,2), "% of all observations are missing!",sep=""),call.=F)
    }
    message("Checking variation in covariates ... ")
    nzv <- caret::nearZeroVar(data[,Z.char], saveMetrics=T)
    if (sum(nzv$zeroVar==T) > 0) {
      stop("Covariate(s) '", paste(row.names(nzv[nzv$zeroVar==T,]), collapse="', '"),
           "' with zero variance(s)!", call.=F)
    } else if (sum(nzv$nzv==T) > 0) {
      warning("Covariate(s) '",paste(row.names(nzv[nzv$nzv==T,]), collapse="', '"),
              "' with near zero variance(s)!",immediate.=T,call.=F)
    }
    message("Checking variation in covariates completed!")
    message("Checking pairwise correlation among covariates ... ")
    cor <- cor(data[,Z.char])
    threshold.cor <- 0.9
    if (sum(abs(cor[upper.tri(cor)])>threshold.cor) > 0) {
      cor[lower.tri(cor,diag=T)] <- 0
      ind <- which(abs(cor)>threshold.cor)
      pairs <- sapply(ind, function(ind) c(rownames(cor)[ind%%NROW(cor)],
                                           colnames(cor)[ind%/%NROW(cor)+1]))
      warning("The following ", NCOL(pairs),
              " pair(s) of covariates are highly correlated (correlation > ",
              threshold.cor,"): ", immediate.=T, call.=F)
      invisible(apply(pairs,2,function(col) message('("',paste(col, collapse='", "'),'")')))
    }
    message("Checking pairwise correlation among covariates completed!")
    ## check VIF
    message("Checking VIF of covariates ... ")
    m.lm <- lm(as.formula(paste(Y.char,"~",paste(Z.char, collapse="+"))), data=data)
    vif <- olsrr::ols_vif_tol(m.lm)
    if(sum(vif$VIF >= 10) > 0){
      warning("Covariate(s) '",
              paste(as.data.frame(vif)[vif$VIF>=10,"Variables"], collapse="', '"),
              "' with serious multicollinearity!",immediate.=T,call.=F)
    }
    message("Checking VIF of covariates completed!")
  }
  data <- data[order(factor(data[,prov.char])),]
  prov.size <- as.integer(table(data[,prov.char]))
  prov.size.long <- rep(prov.size,prov.size)
  data$included <- 1 * (prov.size.long > cutoff)
  warning(sum(prov.size<=cutoff)," out of ",length(prov.size),
          " providers considered small and filtered out!",immediate.=T,call.=F)
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

