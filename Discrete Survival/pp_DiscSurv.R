pp.DiscSurv <- function(data, Event.char, prov.char, Z.char, Time.char, lambda, nlambda = 100, 
                     lambda.min.ratio = 1e-4, penalize.x = rep(1, length(Z.char)), penalized.multiplier, 
                     lambda.early.stop = FALSE, nvar.max = p, stop.dev.ratio = 1e-3, bound = 10.0, backtrack = FALSE, 
                     tol = 1e-4, max.each.iter = 1e4, max.total.iter = (max.each.iter * nlambda), actSet = TRUE, 
                     actIter = max.each.iter, actVarNum = sum(penalize.x == 1), actSetRemove = F, returnX = FALSE, 
                     trace.lambda = FALSE, threads = 1, MM = FALSE, ...){

  ## data structure:  
  ##   Prov.ID        Z_1         Z_2        Z_3        Z_4        Z_5 time status
  ## 1       1  1.0343958 0.003741411  1.4221614 -0.8159581 0.43875657   10      0
  ## 2       1  0.2061218 1.079745297 -0.5828745  0.6271257 0.09106685    5      1
  ## 3       1 -0.6566084 0.153930186 -0.3639364  0.3772989 0.53571441   10      0
  ## 4       1  0.3758861 0.332605759  0.5301672  1.1472704 0.18038985    5      1
  ## 5       1  0.4798284 0.295880942  0.9460193  0.8258463 1.92481395   10      0
  ## 6       1  0.3853669 1.009904976  1.1096166  1.3584786 1.72384165   10      0
  ##  Event.char = "status"
  ##  prov.char = "
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
  
  data <- data[order(factor(data[, prov.char])), ] 
  ID <- as.matrix(data[, prov.char]) # ID vector
  colnames(ID) <- prov.char 
  
  pseudo.group <- 1:length(Z.char)    #each variable forms a group
  if (min(penalize.x) == 0){
    pseudo.group[penalize.x == 0] = 0 # set unpenalized variables
  } 
  
  # failure at each time point
  fit <- survival::survfit(survival::Surv(data[, Time.char], data[, Event.char]) ~ 1)
  sum.failure <- fit$n.event
  KM.baseline.hazard <- c(fit$cumhaz[1], fit$cumhaz[2:length(fit$cumhaz)] - fit$cumhaz[1:(length(fit$cumhaz) - 1)])
  
  failure.each.center <- tapply(data$status, data$Prov.ID, sum)
  
  
  # !!! "standardize = T" may cause problems in transforming gamma and alpha back, current we only consider "standardize = F"
  
  #if (standardize == T){
  #  std.Z <- newZG.Std.grplasso(data, Z.char, pseudo.group, penalized.multiplier)
  #  Z <- std.Z$std.Z[, , drop = F]  # standardized covariate matrix
  #  pseudo.group <- std.Z$g  # new group order
  #  penalized.multiplier <- std.Z$m # new group multiplier
  #} else {
    std.Z <- newZG.Unstd.grplasso(data, Z.char, pseudo.group, penalized.multiplier)
    Z <- std.Z$std.Z[, , drop = F] 
    pseudo.group <- std.Z$g 
    penalized.multiplier <- std.Z$m 
  #}
  
  delta.obs <- data[, Event.char]  #each observation's event variable
  time <- data[, Time.char]
  p <- ncol(Z)
  nvar.max <- as.integer(nvar.max)

  # gamma start from KM estimator
  KM.baseline.hazard.add.small <- KM.baseline.hazard
  KM.baseline.hazard.add.small[which(KM.baseline.hazard.add.small == 0)] <- 1e-10
  gamma <- log(KM.baseline.hazard.add.small/(1 - KM.baseline.hazard.add.small))  
  beta <- rep(0, ncol(Z))
  n.prov <- sapply(split(delta.obs, ID), length) 
  alpha.prov <- rep(log(mean(delta.obs)/(1 - mean(delta.obs))), length(n.prov))
  
  
  
  #####---check-----#####
  if (missing(lambda)) {
    if (nlambda < 2) {
      stop("nlambda must be at least 2", call. = FALSE)
    } else if (nlambda != round(nlambda)){
      stop("nlambda must be a positive integer", call. = FALSE)
    } 
    lambda.fit <- set.lambda.Surv2(delta.obs, Z, time, ID, gamma, beta, alpha.prov, prov.char, 
                                   pseudo.group, penalized.multiplier, nlambda = nlambda,
                                   lambda.min.ratio = lambda.min.ratio)
    lambda.seq <- lambda.fit$lambda.seq
    beta <- lambda.fit$beta
    alpha <- lambda.fit$alpha
    gamma <- lambda.fit$gamma
  } else {
    nlambda <- length(lambda)  # Note: lambda can be a single value
    lambda.seq <- as.vector(sort(lambda, decreasing = TRUE))
  }
  #####---check-----#####
  
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

  fit <- pp_DiscSurv_lasso(delta.obs, max.timepoint, Z, n.prov, time, gamma, beta, alpha.prov, K0, K1, sum.failure, failure.each.center, 
                        lambda.seq, penalized.multiplier, max.total.iter, max.each.iter, tol, backtrack, MM, bound, initial.active.variable, 
                        nvar.max, trace.lambda, threads, actSet, actIter, actVarNum, actSetRemove)
  
  gamma <- fit$gamma
  beta <- fit$beta
  alpha <- fit$alpha
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
  beta <- unorthogonalize(beta, std.Z$std.Z, pseudo.group)
  rownames(beta) <- colnames(Z)
  if (std.Z$reorder == TRUE){  # original order of beta
    beta <- beta[std.Z$ord.inv, , drop = F]
  } 
  #if (standardize == T) {
  #  unstandardize.para <- unstandardize(beta, gamma, std.Z)
  #  beta <- unstandardize.para$beta
  #  gamma <- unstandardize.para$gamma
  #}
  
  # Names
  dimnames(beta) <- list(Z.char, round(lambda, digits = 4))
  if (nrow(gamma) == 1 & length(lambda.seq) == 1){
    gamma <- t(gamma)
    alpha <- t(alpha)
  }
  
  dimnames(gamma) <- list(paste0("T_", 1:count.gamma), round(lambda, digits = 4))
  dimnames(alpha) <- list(names(n.prov), round(lambda, digits = 4))
  
  colnames(eta) <- round(lambda, digits = 4)
  
  result <- structure(list(beta = beta, 
                           alpha = alpha,
                           gamma = gamma,
                           lambda = lambda,
                           linear.components = eta,  #eta = alpha +  Z * beta
                           df = df, 
                           iter = iter,
                           penalized.multiplier = penalized.multiplier,
                           penalize.x = penalize.x),
                      class = "ppDiscSurv")  #define a list for prediction
  if (returnX == TRUE){
    result$returnX <- std.Z
  }
  return(result)
}
