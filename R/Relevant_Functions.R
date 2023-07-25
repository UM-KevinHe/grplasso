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

# "SerBIN" for computing response residuals
SerBIN.residuals <- function(Y, Z, n.prov, gamma.prov, beta){
  fit <- SerBIN(Y, Z, n.prov, gamma.prov, beta)
  gamma.prov <- as.numeric(fit$gamma);
  beta <- as.numeric(fit$beta)
  gamma.obs <- rep(gamma.prov, n.prov)
  Pred <- as.numeric(plogis(gamma.obs+Z %*% beta))
  response.residual <- as.matrix(Y - Pred, ncol = 1)
  colnames(response.residual) <- "response residuals"
  ls <- list(beta = beta, gamma = gamma.prov, residual = response.residual)
  return(ls)
}

# construct lambda sequence & re-initialize beta and gamma
set.lambda.grplasso <- function(Y, Z, ID, group, n.prov, gamma.prov, beta, group.multiplier,
                                nlambda = 100, lambda.min.ratio = 1e-04){
  n <- nrow(Z)
  K <- table(group)
  K1 <- if (min(group) == 0) cumsum(K) else c(0, cumsum(K))
  storage.mode(K1) <- "integer"
  if (K1[1] != 0) {
    fit <- SerBIN.residuals(Y, Z[, group == 0, drop = F], n.prov, gamma.prov, beta[1:sum(group == 0)])
    r <- fit$residual
    beta.initial <- c(fit$beta, rep(0, length(beta) - length(fit$beta)))
    gamma.initial <- fit$gamma
  } else {
    mean.Y <- sapply(split(Y, ID), mean)
    n.prov <- sapply(split(Y, ID), length)
    r <- Y - rep(mean.Y, n.prov)
    beta.initial <- beta
    gamma.initial <- gamma.prov
  }
  lambda.max <- Z_max_grLasso(Z, r, K1, as.double(group.multiplier))/n
  lambda.seq <- exp(seq(log(lambda.max), log(lambda.min.ratio * lambda.max), length = nlambda))
  lambda.seq[1] <- lambda.seq[1] + 1e-5
  ls <- list(beta = beta.initial, gamma = gamma.initial, lambda.seq = lambda.seq)
  return(ls)
}

# "DiscSurv.residuals" is used for computing the sum of response residuals within each individual
DiscSurv.residuals <- function(delta.obs, new.Z, time, alpha, beta.new, n.true_beta){
  fit <- NR_residuals(as.matrix(time), as.matrix(new.Z), as.matrix(delta.obs), as.matrix(alpha),
                      as.matrix(beta.new), tol = 1e-4, max_iter = 1e4)  #.cpp functions
  alpha <- as.numeric(fit$alpha)
  beta <- as.numeric(fit$beta)  #beta here includes original beta and dummy "center effect"
  eta <- as.matrix(new.Z) %*% as.matrix(beta)
  residuals <- DiscSurv_residuals(nrow(new.Z), delta.obs, time, alpha, eta)
  colnames(residuals) <- "Within person residuals"
  beta.out <- ifelse(n.true_beta != 0, 
                 beta[1:n.true_beta],
                 NA)
  ls <- list(beta = beta.out, 
             alpha = alpha, 
             gamma = beta[(1 + n.true_beta):length(beta)],
             residual = residuals)
  return(ls)
}

set.lambda.Surv <- function(delta.obs, Z, time, ID, alpha, beta, gamma.prov, prov.char, group, group.multiplier,
                             nlambda = 100, lambda.min.ratio = 1e-04){
  n <- nrow(Z)
  K <- table(group)
  K1 <- if (min(group) == 0) cumsum(K) else c(0, cumsum(K))
  storage.mode(K1) <- "integer"
  if (K1[1] != 0) {  ## if some beta are unpenalized
    dummy_center <- fastDummies::dummy_cols(ID, select_columns = prov.char, remove_selected_columns = TRUE, remove_first_dummy = TRUE)
    new.Z <- cbind(Z[, group == 0, drop = F], dummy_center) #create new dummy variables for center effect (remove 1st center, which is treated as reference center)
    beta.new <- c(beta[1:sum(group == 0)], gamma.prov[2:(ncol(dummy_center) + 1)])
    n.true_beta <- sum(group == 0)
    fit <- DiscSurv.residuals(delta.obs, new.Z, time, alpha, beta.new, n.true_beta)
    r <- fit$residual

    beta.initial <- c(fit$beta[1:sum(group == 0)], rep(0, length(beta) - sum(group == 0))) 
    gamma.initial <- c(gamma.prov[1], fit$gamma) # the first center effect is not estimated
    alpha.initial <- fit$alpha  #initial time effect
  } else {  ## only time and center
    new.Z <- fastDummies::dummy_cols(ID, select_columns = prov.char, remove_selected_columns = TRUE, remove_first_dummy = TRUE)
    beta.new <- gamma.prov[2:(ncol(new.Z) + 1)]
    n.true_beta <- 0
    fit <- DiscSurv.residuals(delta.obs, new.Z, time, alpha, beta.new, n.true_beta)
    r <- fit$residual

    beta.initial <- beta  #use original beta (all 0)
    gamma.initial <- c(gamma.prov[1], fit$gamma)
    alpha.initial <- fit$alpha  #initial time effect
  }
  lambda.max <- Z_max_grLasso(Z, r, K1, as.double(group.multiplier))/n
  lambda.seq <- exp(seq(log(lambda.max), log(lambda.min.ratio * lambda.max), length = nlambda))
  lambda.seq[1] <- lambda.seq[1] + 1e-5
  ls <- list(beta = beta.initial, alpha = alpha.initial, gamma = gamma.initial, lambda.seq = lambda.seq)
  return(ls)
}


# set up group information
# m: group multiplier, default (missing) is the square root of group size of remaining features
setupG <- function(group, m){
  group.factor <- factor(group)
  if (any(levels(group.factor) == '0')) {
    g <- as.integer(group.factor) - 1
    lev <- levels(group.factor)[levels(group.factor) != '0']
  } else {
    g <- as.integer(group.factor)
    lev <- levels(group.factor)
  }
  if (is.numeric(group) | is.integer(group)) {
    lev <- paste0("G", lev)
  }
  if (missing(m)) {
    m <- rep(NA, length(lev))
    names(m) <- lev
  } else {
    TRY <- try(as.integer(group) == g)
    if (inherits(TRY, 'try-error') || any(!TRY)) stop('Attempting to set group.multiplier is ambiguous if group is not a factor', call. = FALSE)
    if (length(m) != length(lev)) stop("Length of group.multiplier must equal number of penalized groups", call. = FALSE)
    if (storage.mode(m) != "double") storage.mode(m) <- "double"
    if (any(m < 0)) stop('group.multiplier cannot be negative', call.=FALSE)
  }
  # "g" contains the group index of each column, but convert "character" group name into integer
  structure(g, levels = lev, m = m)
}

# remove constant columns if necessary
subsetG <- function(g, nz) { # nz: index of non-constant features
  lev <- attr(g, 'levels')
  m <- attr(g, 'm')
  new <- g[nz] # only include non-constant columns
  dropped <- setdiff(g, new)  # If the entire group has been dropped
  if (length(dropped) > 0) {
    lev <- lev[-dropped] # remaining group
    m <- m[-dropped]
    group.factor <- factor(new) #remaining group factor
    new <- as.integer(group.factor) - 1 * any(levels(group.factor) == '0')  #new group index
  }
  structure(new, levels = lev, m = m)
}

# reorder group index of features
reorderG <- function(g, m) {
  og <- g
  lev <- attr(g, 'levels')
  m <- attr(g, 'm')
  if (any(g == 0)) {
    g <- as.integer(relevel(factor(g), "0")) - 1
  }
  if (any(order(g) != 1:length(g))) {
    reorder <- TRUE
    gf <- factor(g)
    if (any(levels(gf) == "0")) {
      gf <- relevel(gf, "0")
      g <- as.integer(gf) - 1
    } else {
      g <- as.integer(gf)
    }
    ord <- order(g)
    ord.inv <- match(1:length(g), ord)
    g <- g[ord]
  } else {
    reorder <- FALSE
    ord <- ord.inv <- NULL
  }
  structure(g, levels = lev, m = m, ord = ord, ord.inv = ord.inv, reorder = reorder)
}

#Feather-level standardization
standardize.Z <- function(Z){
  mysd <- function(z){
    sqrt(sum((z - mean(z))^2)/length(z))
  }
  new.Z <- scale(as.matrix(Z), scale = apply(as.matrix(Z), 2, mysd))
  center.Z <- attributes(new.Z)$`scaled:center`
  scale.Z <- attributes(new.Z)$`scaled:scale`
  new.Z <- new.Z[, , drop = F]
  res <- list(new.Z = new.Z, center.Z = center.Z, scale.Z = scale.Z)
  return(res)
}

## converting standardized betas back to original variables
unstandardize <- function(beta, gamma, std.Z){
  original.beta <- matrix(0, nrow = length(std.Z$scale), ncol = ncol(beta))
  original.beta[std.Z$nz, ] <- beta / std.Z$scale[std.Z$nz]  # modified beta
  original.gamma <- t(apply(gamma, 1, function(x) x - crossprod(std.Z$center, original.beta))) # modified intercepts (gamma)
  return(list(gamma = original.gamma, beta = original.beta))
}


# Group-level orthogonalization (column in new order, from group_0 to group_max)
orthogonalize <- function(Z, group) {
  z.names <- colnames(Z)
  n <- nrow(Z)
  J <- max(group)
  QL <- vector("list", J)
  orthog.Z <- matrix(0, nrow = nrow(Z), ncol = ncol(Z))
  colnames(orthog.Z) <- z.names
  # unpenalized group will not be orthogonalized
  orthog.Z[, which(group == 0)] <- Z[, which(group == 0)]

  # SVD and generate orthogonalized X
  for (j in seq_along(integer(J))) {
    ind <- which(group == j)
    if (length(ind) == 0) { # skip 0-length group
      next
    }
    SVD <- svd(Z[, ind, drop = FALSE], nu = 0)  # Q matrix (orthonormal matrix of eigenvectors)
    r <- which(SVD$d > 1e-10)  #remove extremely small singular values
    QL[[j]] <- sweep(SVD$v[, r, drop = FALSE], 2, sqrt(n)/SVD$d[r], "*") # Q * Lambda^{-1/2}
    orthog.Z[, ind[r]] <- Z[, ind] %*% QL[[j]]  # group orthogonalized X, where (X^T * X)/n = I
  }
  nz <- !apply(orthog.Z == 0, 2, all)  #find all zero
  orthog.Z <- orthog.Z[, nz, drop = FALSE]
  attr(orthog.Z, "QL") <- QL
  attr(orthog.Z, "group") <- group[nz]
  return(orthog.Z)
}

# convert orthogonalized beta back to original scales
unorthogonalize <- function(beta, Z, group) {
  ind <- !sapply(attr(Z, "QL"), is.null)
  QL <- Matrix::bdiag(attr(Z, "QL")[ind]) #block diagonal matrix
  if (sum(group == 0) > 0){ #some groups are unpenalized
    ind0 <- which(group==0)
    original.beta <- as.matrix(rbind(beta[ind0, , drop = FALSE], QL %*% beta[-ind0, , drop = FALSE]))
  } else {  # all groups are penalized
    original.beta <- as.matrix(QL %*% beta)
  }
  return(original.beta)
}


# standardize + orthogonalize covariate matrix
newZG.Std.grplasso <- function(data, Z.char, g, m){
  Z <- as.matrix(data[, Z.char, drop = F])
  if (any(is.na(Z))){
    stop("Missing data (NA's) detected in covariate matrix!", call. = FALSE)
  }
  if (length(g) != ncol(Z)) {
    stop ("Dimensions of group is not compatible with Z", call. = FALSE)
  }
  G <- setupG(g, m) # setup group
  std <- standardize.Z(Z)
  std.Z <- std[[1]]
  center <- std[[2]]
  scale <- std[[3]]
  nz <- which(scale > 1e-6)   # non-constant columns
  if (length(nz) != ncol(Z)) {
    std.Z <- std.Z[, nz, drop = F]
    G <- subsetG(G, nz)
  }
  # Reorder groups
  G <- reorderG(G, attr(G, 'm'))
  if (attr(G, 'reorder')){
    std.Z <- std.Z[, attr(G, 'ord')]
  }
  # Group-level orthogonalization
  std.Z <- orthogonalize(std.Z, G)
  g <- attr(std.Z, "group")
  # Set group multiplier if missing
  m <- attr(G, 'm')
  if (all(is.na(m))) {
    m <- sqrt(table(g[g != 0]))
  }
  res <- list(std.Z = std.Z, g = g, m = m, reorder = attr(G, 'reorder'), nz = nz,
              ord.inv = attr(G, 'ord.inv'), center = center, scale = scale)
  return(res)
}

# Only orthogonalize covariate matrix
newZG.Unstd.grplasso <- function(data, Z.char, g, m){
  Z <- as.matrix(data[, Z.char, drop = F])
  if (any(is.na(Z))){
    stop("Missing data (NA's) detected in covariate matrix!", call. = FALSE)
  }
  if (length(g) != ncol(Z)) {
    stop ("Dimensions of group is not compatible with Z", call. = FALSE)
  }
  G <- setupG(g, m)
  mysd <- function(x){
    sqrt(sum((x - mean(x))^2)/length(x))
  }
  scale <- apply(as.matrix(Z), 2, mysd)
  nz <- which(scale > 1e-6) #remove constant columns
  if (length(nz) != ncol(Z)) {
    std.Z <- Z[, nz, drop = F]
    G <- subsetG(G, nz)
  } else {
    std.Z <- Z
  }
  G <- reorderG(G, attr(G, 'm'))
  if (attr(G, 'reorder')){
    std.Z <- std.Z[, attr(G, 'ord')]
  }
  std.Z <- orthogonalize(std.Z, G)
  g <- attr(std.Z, "group")
  # Set group multiplier if missing
  m <- attr(G, 'm')
  if (all(is.na(m))) {
    m <- sqrt(table(g[g != 0]))
  }
  res <- list(std.Z = std.Z, g = g, m = m, reorder = attr(G, 'reorder'),
              ord.inv = attr(G, 'ord.inv'), nz = nz)
  return(res)
}


## handle multiple columns of Y
newY <- function(data, Y.char){
  y <- data[, Y.char]
  if (is.data.frame(y)){
    y <- as.matrix(y)
  }
  if (is.matrix(y)) {
    d <- dim(y)
    y <- t(y)
  } else {
    d <- c(length(y), 1)
  }

  # Convert fuzzy binomial data
  if (typeof(y) != "logical") {
    tab <- table(y)
    if (length(tab) > 2) stop("Outcome is not binary", call.=FALSE)
    if (!identical(names(tab), c("0", "1"))) {
      message(paste0("Logistic regression modeling Pr(Y = ", names(tab)[2], ")"))
      y <- as.double(as.character(y) == names(tab)[2])  #convert to 0 & 1
      if (d[2] > 1) {
        attr(y, "dim") <- d
      }
    }
  }

  # Convert to double, if necessary
  if (typeof(y) != "double") {
    tryCatch(storage.mode(y) <- "double", warning = function(w) {stop("Y must be numeric or able to be coerced to numeric", call. = FALSE)})
  }
  if (any(is.na(y))){
    stop("Missing data (NA's) detected in outcome Y!", call. = FALSE)
  }

  # Handle multi
  if (is.matrix(y)) {
    if (ncol(y) > 1) {
      if (is.null(colnames(y))){
        paste("Y", 1:ncol(y), sep = "")
      }
    }
    attributes(y) <- NULL
  }

  attr(y, "m") <- d[2]
  return(y)
}



cvf.pplasso <- function(i, data, Y.char, Z.char, prov.char, fold, cv.args){
  cv.args$data <- data[fold != i, , drop = FALSE]
  cv.args$Y.char <- Y.char
  cv.args$Z.char <- Z.char
  cv.args$prov.char <- prov.char

  fit.i <- do.call("pp.lasso", cv.args)

  data.i <- data[fold == i, , drop = FALSE]
  Y.i <- data.i[, Y.char]
  yhat.i <- matrix(predict(fit.i, data.i, Z.char, prov.char, type = "response"), nrow(data.i)) #y-hat matrix across all given lambda
  loss <- loss.grp.lasso(Y.i, yhat.i)  ## cross entropy loss
  pe <- (yhat.i < 0.5) == Y.i  # wrong prediction
  list(loss = loss, pe = pe, nl = length(fit.i$lambda), yhat = yhat.i)
}



cvf.grplasso <- function(i, data, Y.char, Z.char, prov.char, fold, cv.args){
  cv.args$data <- data[fold != i, , drop = FALSE]
  cv.args$Y.char <- Y.char
  cv.args$Z.char <- Z.char
  cv.args$prov.char <- prov.char

  fit.i <- do.call("grp.lasso", cv.args)

  data.i <- data[fold == i, , drop = FALSE]
  Y.i <- data.i[, Y.char]
  yhat.i <- matrix(predict(fit.i, data.i, Z.char, prov.char, type = "response"), nrow(data.i)) #y-hat matrix across all given lambda
  loss <- loss.grp.lasso(Y.i, yhat.i)  ## cross entropy loss
  pe <- (yhat.i < 0.5) == Y.i  # wrong prediction
  list(loss = loss, pe = pe, nl = length(fit.i$lambda), yhat = yhat.i)
}



cvf.ppDiscSurv <- function(i, data, Event.char, prov.char, Z.char, Time.char, fold, original.count.alpha, cv.args){
  cv.args$data <- data[fold != i, , drop = FALSE]  #current training data (should have all time points
  cv.args$Event.char <- Event.char
  cv.args$prov.char <- prov.char
  cv.args$Z.char <- Z.char
  cv.args$Time.char <- Time.char
  
  # some time point might lost within one sub data set, so we need to reorder time
  #count.gamma <- length(unique(cv.args$data[, Time.char]))
  #if (count.gamma != original.count.gamma){
  #  timepoint.increase <- sort(unique(data[, Time.char]))
  #  # new time start from 1, and time points are {1, 2, 3, ...}
  #  for (i in 1:count.gamma){
  #    data[, Time.char][which(data[, Time.char] == timepoint.increase[i])] <- i
  #  }
  #}
  
  fit.i <- do.call("pp.DiscSurv", cv.args)  #fit the discrete survival model using one training data set (9/10 data)
  data.i <- data[fold == i, , drop = FALSE]  #current validation data
  yhat.i <- predict(fit.i, data.i, Z.char, Time.char, type = "response", return.Array = FALSE) # y-hat matrix across all given lambda; data has been expanded
  
  data.small <- data.i[, c(Event.char, Time.char)]
  Y.i <- discSurv::dataLong(dataShort = data.small, timeColumn = Time.char, 
                            eventColumn = Event.char, timeAsFactor = TRUE)$y
  
  loss <- loss.Disc.Surv(Y.i, yhat.i)  ## cross entropy loss
  list(loss = loss, nl = length(fit.i$lambda), yhat = yhat.i)
}



cvf.DiscSurv <- function(i, data, Event.char, Z.char, Time.char, fold, original.count.gamma, cv.args){
  cv.args$data <- data[fold != i, , drop = FALSE]
  cv.args$Event.char <- Event.char
  cv.args$Z.char <- Z.char
  cv.args$Time.char <- Time.char

  # some time point might lost within one sub data set, so we need to reorder time
  #count.gamma <- length(unique(cv.args$data[, Time.char]))
  #if (count.gamma != original.count.gamma){
  #  timepoint.increase <- sort(unique(data[, Time.char]))
  #  # new time start from 1, and time points are {1, 2, 3, ...}
  #  for (i in 1:count.gamma){
  #    data[, Time.char][which(data[, Time.char] == timepoint.increase[i])] <- i
  #  }
  #}

  fit.i <- do.call("Disc.Surv", cv.args)  #fit the discrete survival model using one training data set
  data.i <- data[fold == i, , drop = FALSE]
  yhat.i <- predict(fit.i, data.i, Z.char, Time.char, type = "response") # y-hat matrix across all given lambda; data has been expanded

  data.small <- data.i[, c(Event.char, Time.char)]
  Y.i <- discSurv::dataLong(dataShort = data.small, timeColumn = Time.char, 
                            eventColumn = Event.char, timeAsFactor = TRUE)$y

  loss <- loss.Disc.Surv(Y.i, yhat.i)  ## cross entropy loss
  list(loss = loss, nl = length(fit.i$lambda), yhat = yhat.i)
}


loss.grp.lasso <- function(Y.i, yhat.i){
  yhat.i[yhat.i < 0.00001] <- 0.00001
  yhat.i[yhat.i > 0.99999] <- 0.99999

  if (is.matrix(yhat.i) == T) {
    val <- matrix(NA, nrow = nrow(yhat.i), ncol = ncol(yhat.i))
    if (sum(Y.i == 1)) {  # if all 1 or all zero, then we only need calculate one of the following "if"
      val[Y.i == 1,] <- - 2 * log(yhat.i[Y.i == 1, , drop = FALSE])
    }
    if (sum(Y.i == 0)){
      val[Y.i == 0,] <- -2 * log(1 - yhat.i[Y.i == 0, , drop = FALSE])
    }
  } else {
    val <- double(length(Y.i))  # a zero vector
    if (sum(Y.i == 1)) {  # if all 1 or all zero, then we only need calculate one of the following "if"
      val[Y.i == 1] <- - 2 * log(yhat.i[Y.i == 1])
    }
    if (sum(Y.i == 0)){
      val[Y.i == 0] <- -2 * log(1 - yhat.i[Y.i == 0])
    }
  }
  return(val)
}


loss.Disc.Surv <- function(Y.i, yhat.i){
  yhat.i[yhat.i < 0.00001] <- 0.00001
  yhat.i[yhat.i > 0.99999] <- 0.99999

  if (is.matrix(yhat.i) == T) {
    val <- matrix(NA, nrow = nrow(yhat.i), ncol = ncol(yhat.i))
    if (sum(Y.i == 1)) {  # if all 1 or all zero, then we only need calculate one of the following "if"
      val[Y.i == 1,] <- - 2 * log(yhat.i[Y.i == 1, , drop = FALSE])
    }
    if (sum(Y.i == 0)){
      val[Y.i == 0,] <- -2 * log(1 - yhat.i[Y.i == 0, , drop = FALSE])
    }
  } else {
    val <- double(length(Y.i))  # a zero vector
    if (sum(Y.i == 1)) {  # if all 1 or all zero, then we only need calculate one of the following "if"
      val[Y.i == 1] <- - 2 * log(yhat.i[Y.i == 1])
    }
    if (sum(Y.i == 0)){
      val[Y.i == 0] <- -2 * log(1 - yhat.i[Y.i == 0])
    }
  }
  return(val)
}

