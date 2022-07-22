## all provider effect will not be penalized
## group lasso algorithm is implemented for penalizing risk factors coefficient; 
  ### "group information" is specified by group = c(1,1,2,1,3,...), same index implies same group
  ### unpenalized risk factor coefficients can be specified by group 0, (i.e. group = c(0,0,1,1,2,...))
## coefficient estimation has been converted back to original scale

grp.lasso <- function(data, Y.char, Z.char, prov.char, method = c("lasso", "grlasso"), group = 1:ncol(Z), group.multiplier, 
                      standardize = T, lambda, nlambda = 100, lambda.min.ratio = 1e-4, lambda.early.stop = FALSE, nvar.max = p, 
                      group.max = length(unique(group)), stop.dev.ratio = 1e-3, bound = 10.0, backtrack = FALSE, tol = 1e-4, 
                      max.iter = 10000, returnX = FALSE, trace.lambda = FALSE, threads = 1, MM = FALSE, ...){
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
    std.Z <- newZG.Std(data, Z.char, group, group.multiplier)
    Z <- std.Z$std.Z[, ]  # standardized covariate matrix
    group <- std.Z$g  # new group order
    group.multiplier <- std.Z$m # new group multiplier
  } else {
    std.Z <- newZG.Unstd(data, Z.char, group, group.multiplier)
    Z <- std.Z$std.Z[, ] 
    group <- std.Z$g 
    group.multiplier <- std.Z$m 
  }
  
  p <- ncol(Z)
  Y <- newY(data, Y.char)
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
    lambda.fit <- set.lambda(Y, Z, ID, group, n.prov, gamma.prov, beta, group.multiplier,
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
  # first penalized group index (start from 0)
  K0 <- as.integer(if (min(group) == 0) K[1] else 0) 
  # end index of each group (start from 1; include group 0)
  # e.g. If group1 = c(0, 0, 0, 0, 1, 2, 2, 3, 4, 4), then K0 = 4 and K1 = c(4, 5, 7, 8, 10); 
  # If group2 = c(1, 2, 2, 3, 4, 4), then K0 = 0 and K1 = c(0, 1, 3, 4, 6)
  K1 <- as.integer(if (min(group) == 0) cumsum(K) else c(0, cumsum(K)))  
  if (K0 == 0){
    # Since our stop rule is based on beta change, we need at least one beta to be "updated" during each iteration
    # If all beta's are penalized, we set the first group with smallest group size as the first active group.
    initial.active.group <- which(K == min(K))[1] - 1  
  } else {
    # If some beta's are unpenalized, those ones will be "updated" during each iteration. So we don't need to modify initial active group
    initial.active.group <- 0
  }
  
  # main algorithm
  fit <- grp_lasso(Y, Z, n.prov, gamma.prov, beta, K0, K1, lambda.seq, lambda.early.stop, stop.dev.ratio, 
                   group.multiplier, max.iter, tol, nullDev, backtrack, bound, initial.active.group, 
                   nvar.max, group.max, trace.lambda, single.intercept, threads)
  
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
  
  if (iter[1] == max.iter){
    stop("Algorithm failed to converge for any values of lambda", call. = FALSE)
  }
  if (sum(iter) == max.iter){
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


cv.grp.lasso <- function(data, Y.char, Z.char, prov.char, group = 1:ncol(Z), ..., nfolds = 10, 
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
    
    res <- cvf(i, data, Y.char, Z.char, prov.char, fold, cv.args)
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