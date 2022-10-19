grp.lasso <- function(data, Y.char, Z.char, prov.char, group = 1:length(Z.char), group.multiplier, 
                      standardize = T, lambda, nlambda = 100, lambda.min.ratio = 1e-4, lambda.early.stop = FALSE, 
                      nvar.max = p, group.max = length(unique(group)), stop.dev.ratio = 1e-3, bound = 10.0, 
                      backtrack = FALSE, tol = 1e-4, max.each.iter = 1e4, max.total.iter = (max.each.iter * nlambda), 
                      actSet = TRUE, actIter = max.each.iter, actGroupNum = sum(unique(group) != 0), actSetRemove = F,
                      returnX = FALSE, trace.lambda = FALSE, threads = 1, ...){
  if (!is.null(data$included)){  # data after using preparation function
    data <- data[data$included == 1, ]
  }
  
  if (missing(prov.char)){ #single intercept
    ID <- matrix(1, nrow = nrow(data))
    colnames(ID) <- "intercept"
    prov.char <- "intercept"
    single.intercept <- TRUE
  } else {
    data <- data[order(factor(data[, prov.char])),] #data should be sorted by "prov.char"
    ID <- as.matrix(data[, prov.char])
    colnames(ID) <- prov.char 
    single.intercept <- FALSE
  }
  
  initial.group <- group
  if (standardize == T){
    std.Z <- newZG.Std.grplasso(data, Z.char, group, group.multiplier)
    Z <- std.Z$std.Z[, , drop = F]  # standardized covariate matrix
    group <- std.Z$g  # new group order
    group.multiplier <- std.Z$m # new group multiplier
  } else {
    std.Z <- newZG.Unstd.grplasso(data, Z.char, group, group.multiplier)
    Z <- std.Z$std.Z[, , drop = F] 
    group <- std.Z$g 
    group.multiplier <- std.Z$m 
  }
  Y <- newY(data, Y.char)
  
  p <- ncol(Z)
  nvar.max <- as.integer(nvar.max)
  group.max <- as.integer(group.max)
  
  n.prov <- sapply(split(Y, ID), length) 
  gamma.prov <- rep(log(mean(Y)/(1 - mean(Y))), length(n.prov))
  beta <- rep(0, ncol(Z))
  
  if (missing(lambda)) {
    if (nlambda < 2) {
      stop("nlambda must be at least 2", call. = FALSE)
    } else if (nlambda != round(nlambda)){
      stop("nlambda must be a positive integer", call. = FALSE)
    } 
    lambda.fit <- set.lambda.grplasso(Y, Z, ID, group, n.prov, gamma.prov, beta, group.multiplier, 
                                      nlambda = nlambda, lambda.min.ratio = lambda.min.ratio)
    lambda.seq <- lambda.fit$lambda.seq
    beta <- lambda.fit$beta
    gamma.prov <- lambda.fit$gamma
  } else {
    nlambda <- length(lambda)  # Note: lambda can be a single value
    lambda.seq <- as.vector(sort(lambda, decreasing = TRUE))
  }
  
  # nullDev: theoretical deviance of the model that only contains provider effects
  mean.Y.obs <- rep(sapply(split(Y, ID), mean), n.prov) 
  nullDev <- Deviance(Y, mean.Y.obs)
  
  K <- as.integer(table(group)) #number of features in each group
  K0 <- as.integer(if (min(group) == 0) K[1] else 0) 
  K1 <- as.integer(if (min(group) == 0) cumsum(K) else c(0, cumsum(K)))  
  
  initial.active.group <- -1
  if (actSet == TRUE){
    if (K0 == 0){
      initial.active.group <- which(K == min(K))[1] - 1  
    }
  } else {
    actIter <- max.each.iter
  }
  
  # main algorithm
  fit <- grp_lasso(Y, Z, n.prov, gamma.prov, beta, K0, K1, lambda.seq, lambda.early.stop, stop.dev.ratio, group.multiplier, 
                   max.total.iter, max.each.iter, tol, nullDev, backtrack, bound, initial.active.group, nvar.max, group.max, 
                   trace.lambda, single.intercept, threads, actSet, actIter, actGroupNum, actSetRemove)
  
  gamma <- fit$gamma
  beta <- fit$beta
  loss <- fit$Deviance
  eta <- fit$Eta
  df <- fit$Df
  iter <- fit$iter
  
  # Eliminate saturated lambda values
  ind <- !is.na(iter)
  lambda <- lambda.seq[ind]
  beta <- beta[, ind, drop = FALSE]
  gamma <- gamma[, ind, drop = FALSE]
  loss <- loss[ind]
  eta <- eta[, ind, drop = FALSE]
  df <- df[ind]
  iter <- iter[ind]
  
  if (iter[1] == max.total.iter){
    stop("Algorithm failed to converge for any values of lambda", call. = FALSE)
  }
  if (sum(iter) == max.total.iter){
    warning("Algorithm failed to converge for all values of lambda", call. = FALSE)
  }
  
  
  # Original scale
  beta <- unorthogonalize(beta, std.Z$std.Z, group)
  rownames(beta) <- colnames(Z)
  if (std.Z$reorder == TRUE){  # original order of beta
    beta <- beta[std.Z$ord.inv, , drop = F]
  } 
  if (standardize == T) {
    unstandardize.para <- unstandardize(beta, gamma, std.Z)
    beta <- unstandardize.para$beta
    gamma <- unstandardize.para$gamma
  }
  
  # Names
  dimnames(beta) <- list(Z.char, round(lambda, digits = 4))
  if (nrow(gamma) == 1 & length(lambda.seq) == 1){
    gamma <- t(gamma)
  }
  dimnames(gamma) <- list(names(n.prov), round(lambda, digits = 4))
  colnames(eta) <- round(lambda, digits = 4)
  
  result <- structure(list(beta = beta, 
                           gamma = gamma,
                           group = factor(initial.group), 
                           lambda = lambda,
                           loss = loss,
                           linear.predictors = eta,
                           df = df, 
                           iter = iter,
                           group.multiplier = group.multiplier),
                      class = "gr_ppLasso")  #define a list for prediction
  if (returnX == TRUE){
    result$returnX <- std.Z
  }
  return(result)
}


cv.grp.lasso <- function(data, Y.char, Z.char, prov.char, group = 1:length(Z.char), ..., nfolds = 10, 
                         seed, fold, trace.cv = FALSE){
  if (missing(prov.char)){ #single intercept
    prov.char <- "intercept"
    data$intercept <- matrix(1, nrow = nrow(data))
  }
  
  # "...": additional arguments to "grp.lasso"
  # "fold": a vector that specifies the fold that observations belongs to
  fit.args <- list(...)
  fit.args$data <- data
  fit.args$Y.char <- Y.char
  fit.args$Z.char <- Z.char
  fit.args$prov.char <- prov.char
  fit.args$group <- group
  fit.args$returnX <- TRUE
  
  fit <- do.call("grp.lasso", fit.args)  #fit the entire model
  
  # get standardized Z
  ZG <- fit$returnX
  Z <- ZG$std.Z
  y <- as.matrix(data[, Y.char])
  returnX <- list(...)$returnX
  if (is.null(returnX) || !returnX){ 
    fit$returnX <- NULL
  } 
  
  #setup folds
  if (!missing(seed)) {
    set.seed(seed)
  }
  
  n.obs <- length(y)
  
  if (missing(fold)){ # make sure each fold contains same proportion of Y = 0 and Y = 1
    ind1 <- which(y == 1)
    ind0 <- which(y == 0)
    n1 <- length(ind1)
    n0 <- length(ind0)
    fold1 <- 1:n1 %% nfolds  # assign observations to different folds
    fold0 <- (n1 + 1:n0) %% nfolds
    fold1[fold1==0] <- nfolds # set fold "0" to fold max
    fold0[fold0==0] <- nfolds
    fold <- integer(n.obs) 
    fold[y == 1] <- sample(fold1)
    fold[y == 0] <- sample(fold0)
  } else {
    nfolds <- max(fold)
  }
  
  # Do Cross-Validation
  E <- Y <- matrix(NA, nrow = length(y), ncol = length(fit$lambda))
  PE <- E
  cv.args <- list(...)
  cv.args$lambda <- fit$lambda
  cv.args$group <- ZG$g
  cv.args$group.multiplier <- ZG$m
  
  for (i in 1:nfolds) {
    if (trace.cv == TRUE){
      cat("Starting CV fold #", i, sep="","...\n")
    }
    
    res <- cvf.grplasso(i, data, Y.char, Z.char, prov.char, fold, cv.args)
    Y[fold == i, 1:res$nl] <- res$yhat
    E[fold == i, 1:res$nl] <- res$loss
    PE[fold == i, 1:res$nl] <- res$pe  # wrong predict class
  }
  
  # Eliminate saturated lambda values, if any
  ind <- which(apply(is.finite(E), 2, all))
  E <- E[, ind, drop=FALSE]
  Y <- Y[, ind]
  lambda <- fit$lambda[ind]
  
  cve <- apply(E, 2, mean)  # mean cross entropy loss within each lambda value
  cvse <- apply(E, 2, sd) / sqrt(n.obs) # standardized loss
  min <- which.min(cve)  #find index of lambda with minimum cve
  pe <- apply(PE[, ind], 2, mean)  # mean predict class error
  
  result <- structure(list(cve = cve, 
                           cvse = cvse, 
                           lambda = lambda, 
                           fit = fit, #model with entire data
                           pe = pe,
                           fold = fold, 
                           min = min, 
                           lambda.min = lambda[min]),
                      class = "cv.gr_ppLasso")
  return(result)
}