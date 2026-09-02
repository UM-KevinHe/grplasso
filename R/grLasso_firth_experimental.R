## ---------------------------------------------------------------------------
## grpLassoFirth() -- EXPERIMENTAL, NOT EXPORTED, NOT FINISHED
## ---------------------------------------------------------------------------
##
## This was formerly R/grLasso2.R, where it defined a second function also
## called grp.lasso().  R sources R/ alphabetically, so grLasso2.R silently
## overrode grLasso.R and this half-finished version was the one users got --
## renaming either file would have flipped the behaviour of the package with
## no other visible change.  It is renamed here so that grp.lasso() in
## R/grLasso.R is unambiguously the exported implementation.
##
## WHAT IT ADDS
##   Firth's bias correction on the provider effects gamma.  When a provider's
##   outcomes are completely separated (all 0 or all 1) its gamma diverges
##   under maximum likelihood; the Firth penalty is supposed to pull it back to
##   a finite value.  The correction itself lives in C++ (FirthGamma, in
##   src/); this function only threads a `firth` flag through to
##   set.lambda.grplasso() and the compiled grp_lasso().  Everything else is
##   identical to grp.lasso().
##
## WHY IT IS NOT EXPORTED (measured, not assumed)
##   1. It fails on exactly the data it exists for.  With one provider forced
##      to all Y = 1, firth = TRUE aborts with "median(): detected NaN", while
##      firth = FALSE returns gamma = 8.56 (diverging, clipped by `bound`).
##      On ordinary data it runs and does shift the fit (max |dgamma| ~ 0.20,
##      max |dbeta| ~ 0.013), so the flag is live, not inert.
##   2. Nothing downstream supports it.  cv.grp.lasso(), pp.lasso(),
##      Strat.cox(), DiscSurv() and pp.DiscSurv() have no `firth` argument, so
##      a Firth fit cannot be cross-validated or compared within the package.
##
## TO RESUME DEVELOPMENT
##   - Trace the NaN: it surfaces at the clamp gamma = clamp(gamma,
##     median(gamma) +- bound), so a NaN is already present in gamma by then.
##     FirthGamma() is where to look.
##   - Then add `firth` to cv.grp.lasso() and its cvf helper.
##   - Then export it (or fold it back into grp.lasso() as an option), give it
##     its own roxygen block, and add a test for the separated-provider case.
##
## Called as grplasso:::grpLassoFirth(...) meanwhile.  Its return value is a
## "gr_ppLasso" object, so coef(), plot() and predict() already work on it.
## Kept inside R/ deliberately: R CMD check parses it on every build, so it
## cannot rot silently the way code parked outside the package would.
## ---------------------------------------------------------------------------

grpLassoFirth <- function(data, Y.char, Z.char, prov.char, group = 1:length(Z.char), group.multiplier,
                      standardize = TRUE, lambda, nlambda = 100, lambda.min.ratio = 1e-4, lambda.early.stop = FALSE,
                      nvar.max = p, group.max = length(unique(group)), stop.dev.ratio = 1e-3, bound = 10.0,
                      backtrack = FALSE, tol = 1e-4, max.each.iter = 1e4, max.total.iter = (max.each.iter * nlambda),
                      actSet = TRUE, actIter = max.each.iter, actGroupNum = sum(unique(group) != 0), actSetRemove = FALSE,
                      returnX = FALSE, trace.lambda = FALSE, threads = 1, firth = FALSE, ...){
# *** END MODIFIED ***
  #if (!is.null(data$included)){  # data after using preparation function
  #  data <- data[data$included == 1, ]
  #}

  if (missing(prov.char)){ #single intercept
    warning("Provider information not provided. All data is assumed to originate from a single provider!", call. = FALSE)
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
  if (standardize == TRUE){
    std.Z <- newZG.Std.grplasso(data, Z.char, group, group.multiplier)
    Z <- std.Z$std.Z[, , drop = FALSE]  # standardized covariate matrix
    group <- std.Z$g  # new group order
    group.multiplier <- std.Z$m # new group multiplier
  } else {
    std.Z <- newZG.Unstd.grplasso(data, Z.char, group, group.multiplier)
    Z <- std.Z$std.Z[, , drop = FALSE]
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
                                      nlambda = nlambda, lambda.min.ratio = lambda.min.ratio,
                                      firth = firth)  # *** MODIFIED: pass firth ***
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

  # main algorithm (note: "MM" algorithm must be used for group lasso problem)
  # *** MODIFIED: pass firth to grp_lasso ***
  fit <- grp_lasso(Y, Z, n.prov, gamma.prov, beta, K0, K1, lambda.seq, lambda.early.stop, stop.dev.ratio, group.multiplier,
                   max.total.iter, max.each.iter, tol, nullDev, backtrack, bound, initial.active.group, nvar.max, group.max,
                   trace.lambda, single.intercept, threads, actSet, actIter, actGroupNum, actSetRemove,
                   firth)  # *** MODIFIED: added firth ***

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
    beta <- beta[std.Z$ord.inv, , drop = FALSE]
  }
  if (standardize == TRUE) {
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
