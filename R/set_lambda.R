#####----- pplasso & grplasso -----#####
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

#####----- pp_DiscSurv -----#####
set.lambda.pp_DiscSurv <- function(delta.obs, Z, time, ID, alpha, beta, gamma.prov, prov.char, group,
                                   group.multiplier, nlambda = 100, lambda.min.ratio = 1e-04){
  n <- nrow(Z)
  K <- table(group)
  K1 <- if (min(group) == 0) cumsum(K) else c(0, cumsum(K))
  storage.mode(K1) <- "integer"
  if (K1[1] != 0) {  ## if some beta are unpenalized
    dummy_center <- fastDummies::dummy_cols(ID, select_columns = prov.char, remove_selected_columns = TRUE, 
                                            remove_first_dummy = TRUE)
    new.Z <- cbind(Z[, group == 0, drop = F], dummy_center) #create new dummy variables for center effect (remove 1st center, which is treated as reference center)
    beta.new <- c(beta[1:sum(group == 0)], gamma.prov[2:(ncol(dummy_center) + 1)])
    n.true_beta <- sum(group == 0)
    fit <- pp.DiscSurv.residuals(delta.obs, new.Z, time, alpha, beta.new, n.true_beta)
    r <- fit$residual
    
    beta.initial <- c(fit$beta[1:sum(group == 0)], rep(0, length(beta) - sum(group == 0))) 
    gamma.initial <- c(gamma.prov[1], fit$gamma) # the first center effect is not estimated
    alpha.initial <- fit$alpha  #initial time effect
  } else {  ## only time and center
    new.Z <- fastDummies::dummy_cols(ID, select_columns = prov.char, remove_selected_columns = TRUE, 
                                     remove_first_dummy = TRUE)
    beta.new <- gamma.prov[2:(ncol(new.Z) + 1)]
    n.true_beta <- 0
    fit <- pp.DiscSurv.residuals(delta.obs, new.Z, time, alpha, beta.new, n.true_beta)
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

#####----- DiscSurv -----#####
set.lambda.DiscSurv <- function(delta.obs, Z, time, alpha, beta, group, group.multiplier, 
                                nlambda = 100, lambda.min.ratio = 1e-04){
  n <- nrow(Z)
  K <- table(group)
  K1 <- if (min(group) == 0) cumsum(K) else c(0, cumsum(K))
  storage.mode(K1) <- "integer"
  if (K1[1] != 0) {  ## use Di's code
    n.true_beta <- sum(group == 0)
    fit <- DiscSurv.residuals(delta.obs, Z[, group == 0, drop = F], time, alpha, beta[1:sum(group == 0)])
    r <- fit$residual
    beta.initial <- c(fit$beta[1:sum(group == 0)], rep(0, length(beta) - sum(group == 0)))
    alpha.initial <- fit$alpha
  } else {  ## use KM results
    plogis.alpha <- plogis(alpha)
    r <- rep(0, n)
    for (i in 1:n){
      r[i] <- delta.obs[i] - sum(plogis.alpha[1:time[i]])
    }
    beta.initial <- beta
    alpha.initial <- alpha
  }
  lambda.max <- Z_max_grLasso(Z, r, K1, as.double(group.multiplier))/n
  lambda.seq <- exp(seq(log(lambda.max), log(lambda.min.ratio * lambda.max), length = nlambda))
  lambda.seq[1] <- lambda.seq[1] + 1e-5
  ls <- list(beta = beta.initial, alpha.initial = alpha, lambda.seq = lambda.seq)
  return(ls)
}

#####----- stratified cox -----#####
set.lambda.cox <- function(delta.obs, Z, time, ID, beta, group, group.multiplier, n.each_prov,
                           nlambda = 100, lambda.min.ratio = 1e-04){
  K <- table(group)
  K1 <- if (min(group) == 0) cumsum(K) else c(0, cumsum(K))
  storage.mode(K1) <- "integer"
  n <- sum(n.each_prov)
  if (K1[1] != 0) {  ## some covariates are not penalized (use cox stratified)
    nullFit <- survival::coxph(survival::Surv(time, delta.obs) ~ Z[, group == 0, drop = F] + strata(ID))
    eta <- nullFit$linear.predictors
    beta.initial <- c(nullFit$beta, rep(0, length(beta) - length(nullFit$beta)))
    rsk <- c()
    for (i in unique(ID)){
      rsk <- c(rsk, rev(cumsum(rev(exp(eta[ID == i])))))
    }
    r <- delta.obs - exp(eta) * cumsum(delta.obs / rsk)
  } else { ## all covariates are penalized
    w <- c()
    h <- c()
    for (i in 1:nrow(unique(ID))){
      temp.w <- 1 / (n.each_prov[i] - (1:n.each_prov[i]) + 1)
      w <- c(w, temp.w)
      h <- c(h, cumsum(delta.obs[ID == i] * temp.w))
    }
    
    r <- delta.obs - h
    beta.initial <- beta
  }
  lambda.max <- Z_max_grLasso(Z, r, K1, as.double(group.multiplier))/n
  lambda.seq <- exp(seq(log(lambda.max), log(lambda.min.ratio * lambda.max), length = nlambda))
  lambda.seq[1] <- lambda.seq[1] + 1e-5
  ls <- list(beta = beta.initial, lambda.seq = lambda.seq)
  return(ls)
}