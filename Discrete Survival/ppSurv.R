pp.Surv <- function(data, Event.char, Z.char, Time.char, standardize = T, lambda, nlambda = 100, lambda.min.ratio = 1e-4, 
                    penalize.x = rep(1, length(Z.char)), penalized.multiplier, lambda.early.stop = FALSE, nvar.max = p, 
                    stop.dev.ratio = 1e-3, bound = 10.0, backtrack = FALSE, tol = 1e-4, max.each.iter = 1e4, 
                    max.total.iter = (max.each.iter * nlambda), actSet = TRUE, actIter = max.each.iter, actVarNum = sum(penalize.x == 1), 
                    actSetRemove = F, returnX = FALSE, trace.lambda = FALSE, threads = 1, MM = FALSE, ...){

  ## data structure:  
  ##   status         Z1         Z2         Z3 Z4 Z5 Z6   time
  ##        0 -0.9144950 -1.9930845 -2.1065963  0  1  1    5.9
  ##        0 -0.8693130  0.4488433  0.7767090  1  1  1    1.3
  ##        1  0.8318571  1.0570449  1.9481646  1  0  0    1.1
  ##        1  0.0822508  1.2584988  2.3215286  0  0  0    1.8
  ##        1  1.1215585 -0.2789035 -0.4760478  1  0  0    0.3
  ##        1 -0.7457385 -0.6656072 -0.6108483  1  0  0    1.1
  ##  Event.char = "status"
  ##  Z.char = c("Z1", "Z2", "Z3", "Z4", "Z5", "Z6")  
  ##  Time.char = "time"
  
  # Convert the observed time in discrete intervals
  # "time" is converted into order index (integers)
  count.gamma <- length(unique(data[, Time.char]))
  timepoint.increase <- sort(unique(data[, Time.char]))
  # new time start from 1, and time points are {1, 2, 3, ...}
  for (i in 1:count.gamma){
    data[, Time.char][which(data[, Time.char] == timepoint.increase[i])] <- i
  }
  max.timepoint <- length(timepoint.increase)  # the number of gamma that we need
  
  pseudo.group <- 1:length(Z.char)    #each variable forms a group
  if (min(penalize.x) == 0){
    pseudo.group[penalize.x == 0] = 0 # set unpenalized variables
  } 
  
  #计算每个时间点的failure人数
  fit <- survival::survfit(Surv(data[, Time.char], data[, Event.char]) ~ 1)
  sum.failure <- fit$n.event
  KM.baseline.hazard <- c(fit$cumhaz[1], fit$cumhaz[2:length(fit$cumhaz)] - fit$cumhaz[1:(length(fit$cumhaz) - 1)])
  
  if (standardize == T){
    std.Z <- newZG.Std.grplasso(data, Z.char, pseudo.group, penalized.multiplier)
    Z <- std.Z$std.Z[, , drop = F]  # standardized covariate matrix
    pseudo.group <- std.Z$g  # new group order
    penalized.multiplier <- std.Z$m # new group multiplier
  } else {
    std.Z <- newZG.Unstd.grplasso(data, Z.char, pseudo.group, penalized.multiplier)
    Z <- std.Z$std.Z[, , drop = F] 
    pseudo.group <- std.Z$g 
    penalized.multiplier <- std.Z$m 
  }
  
  delta.obs <- data[, Event.char]  #each observation's event variable
  time <- data[, Time.char]
  p <- ncol(Z)
  nvar.max <- as.integer(nvar.max)

  # gamma start from KM estimator
  KM.baseline.hazard.add.small <- KM.baseline.hazard
  KM.baseline.hazard.add.small[which(KM.baseline.hazard.add.small == 0)] <- 1e-10
  gamma <- log(KM.baseline.hazard.add.small/(1 - KM.baseline.hazard.add.small))  
  beta <- rep(0, ncol(Z))
  
  if (missing(lambda)) {
    if (nlambda < 2) {
      stop("nlambda must be at least 2", call. = FALSE)
    } else if (nlambda != round(nlambda)){
      stop("nlambda must be a positive integer", call. = FALSE)
    } 
    lambda.fit <- set.lambda.Surv(delta.obs, Z, time, gamma, beta, pseudo.group, 
                                  penalized.multiplier, nlambda = nlambda, 
                                  lambda.min.ratio = lambda.min.ratio)
    lambda.seq <- lambda.fit$lambda.seq
    beta <- lambda.fit$beta
    gamma <- lambda.fit$gamma
  } else {
    nlambda <- length(lambda)  # Note: lambda can be a single value
    lambda.seq <- as.vector(sort(lambda, decreasing = TRUE))
  }
  
  K <- as.integer(table(pseudo.group)) #number of features in each group
  K0 <- as.integer(if (min(pseudo.group) == 0) K[1] else 0) 
  K1 <- as.integer(if (min(pseudo.group) == 0) cumsum(K) else c(0, cumsum(K)))  
  
  initial.active.variable <- -1 ## "-1" means this variable can be any value and we will not use it in our cpp function.
  if (actSet == TRUE){
    if (K0 == 0){ #all variables are penalized
      initial.active.variable <- which(K == min(K))[1] - 1  ## which is actually the first variable
    }
  } else {  ## if we don't use active set method, then the initial active set should contain all penalized variables
    actIter <- max.each.iter
  }
  
  # main algorithm
  fit <- pp_Surv_lasso(delta.obs, max.timepoint, Z, time, gamma, beta, K0, K1, sum.failure, lambda.seq, 
                       penalized.multiplier, max.total.iter, max.each.iter, tol, 
                       backtrack, MM, bound, initial.active.variable, nvar.max, trace.lambda, threads, 
                       actSet, actIter, actVarNum, actSetRemove)
  
  gamma <- fit$gamma
  beta <- fit$beta
  eta <- fit$Eta
  df <- fit$Df
  iter <- fit$iter
  
  # Eliminate saturated lambda values
  ind <- !is.na(iter)
  lambda <- lambda.seq[ind]
  beta <- beta[, ind, drop = FALSE]
  gamma <- gamma[, ind, drop = FALSE]
  eta <- eta[, ind, drop = FALSE]
  df <- df[ind]
  iter <- iter[ind]
  
  #if (iter[1] == max.total.iter){
  #  stop("Algorithm failed to converge for any values of lambda", call. = FALSE)
  #}
  if (sum(iter) == max.total.iter){
    warning("Algorithm failed to converge for all values of lambda", call. = FALSE)
  }
  
  
  # Original scale
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
  
  dimnames(gamma) <- list(paste0("T_", 1:count.gamma), round(lambda, digits = 4))
  colnames(eta) <- round(lambda, digits = 4)
  
  result <- structure(list(beta = beta, 
                           gamma = gamma,
                           lambda = lambda,
                           linear.components = eta,  #eta =  X * beta
                           df = df, 
                           iter = iter,
                           penalized.multiplier = penalized.multiplier,
                           penalize.x = penalize.x),
                      class = "ppSurv")  #define a list for prediction
  if (returnX == TRUE){
    result$returnX <- std.Z
  }
  return(result)
}



cv.pp.Surv <- function(data, Event.char, Z.char, Time.char, penalize.x = rep(1, length(Z.char)), ..., nfolds = 10, 
                       seed, fold, trace.cv = FALSE){
  # "...": additional arguments to "grp.lasso"
  # "fold": a vector that specifies the fold that observations belongs to
  fit.args <- list(...)
  fit.args$data <- data
  fit.args$Event.char <- Event.char
  fit.args$Z.char <- Z.char
  fit.args$Time.char <- Time.char
  fit.args$penalize.x <- penalize.x
  fit.args$returnX <- TRUE
  
  fit <- do.call("pp.Surv", fit.args)  #fit the entire model
  
  # get standardized Z
  ZG <- fit$returnX
  Z <- ZG$std.Z
  delta.obs <- as.matrix(data[, Event.char])
  returnX <- list(...)$returnX
  if (is.null(returnX) || !returnX){ 
    fit$returnX <- NULL
  } 
  
  #setup folds
  if (!missing(seed)) {
    set.seed(seed)
  }
  n.obs <- length(delta.obs)
  original.count.gamma <- length(unique(data[, Time.char])) 
  
  for (s in 1:100){
    if (missing(fold)){ # make sure each fold contains same proportion of censor and failure
      ind1 <- which(delta.obs == 1)
      ind0 <- which(delta.obs == 0)
      n1 <- length(ind1)
      n0 <- length(ind0)
      fold1 <- 1:n1 %% nfolds  # assign observations to different folds
      fold0 <- (n1 + 1:n0) %% nfolds
      fold1[fold1 == 0] <- nfolds # set fold "0" to fold max
      fold0[fold0 == 0] <- nfolds
      fold <- integer(n.obs) 
      fold[delta.obs == 1] <- sample(fold1)
      fold[delta.obs == 0] <- sample(fold0)
    } else {
      nfolds <- max(fold)
    }
    
    #check if some time points have lost within one subset
    len <- rep(NA, nfolds)
    for (i in 1:nfolds){
      len[i] <- length(unique(data[fold != i, Time.char])) != original.count.gamma
    }
    
    if (sum(len) == 0) {
      break
    }# else {
    #  print(paste0("attempt ", s, " failed..."))
    #}
  }
  
  data.small <- cbind(fold, data[, c(Event.char, Time.char)])
  expand.fold <- discSurv::dataLong(dataShort = data.small, timeColumn = Time.char, eventColumn = Event.char, timeAsFactor = TRUE)$fold
  
  # Do Cross-Validation
  E <- Y <- matrix(NA, nrow = sum(data[, Time.char]), ncol = length(fit$lambda)) # stored as expanded matrix
  original.count.gamma <- length(unique(data[, Time.char])) 
  
  cv.args <- list(...)
  cv.args$lambda <- fit$lambda
  cv.args$group <- ZG$g
  cv.args$group.multiplier <- ZG$m
  
  for (i in 1:nfolds) {
    if (trace.cv == TRUE){
      cat("Starting CV fold #", i, sep = "", "...\n")
    }
    res <- cvf.ppSurv(i, data, Event.char, Z.char, Time.char, fold, original.count.gamma, cv.args)
    Y[expand.fold == i, 1:res$nl] <- res$yhat
    E[expand.fold == i, 1:res$nl] <- res$loss
  }
  
  # Eliminate saturated lambda values, if any
  ind <- which(apply(is.finite(E), 2, all))
  E <- E[, ind, drop=FALSE]
  Y <- Y[, ind]
  lambda <- fit$lambda[ind]
  
  cve <- apply(E, 2, mean)  # mean cross entropy loss within each lambda value
  cvse <- apply(E, 2, sd) / sqrt(n.obs) # standardized loss
  min <- which.min(cve)  #find index of lambda with minimum cve
  
  result <- structure(list(cve = cve, 
                           cvse = cvse, 
                           lambda = lambda, 
                           fit = fit, #model with entire data
                           fold = fold, 
                           min = min, 
                           lambda.min = lambda[min]),
                      class = "cv.ppSurv")
  return(result)
}


