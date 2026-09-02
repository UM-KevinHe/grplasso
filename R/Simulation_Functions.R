##################---------LINEAR MODEL----------##################
## linear_data <- sim.linear(m = 100, n.beta = 10, prov.size.mean = 80, rho = 0.7)
#' Simulate clustered data with a continuous outcome
#'
#' Generates a linear model in which each provider contributes its own random
#' intercept and the covariates are correlated with that intercept.
#'
#' @param m number of providers.
#' @param n.beta number of covariates.
#' @param prov.size.mean mean provider size; sizes are Poisson draws floored at 11.
#' @param rho correlation between the covariates and the provider effect, and
#'   between any two covariates.
#'
#' @return A list with the simulated `data`, the `true.beta` and `true.gamma`
#'   used to generate it, and the column names `prov.char`, `Z.char`, `Y.char`.
#'
#' @examples
#' set.seed(1)
#' sim <- sim.linear(m = 10, n.beta = 5, prov.size.mean = 30)
#' str(sim$data)
#' sim$true.beta
#'
#' @export
sim.linear <- function(m, n.beta, prov.size.mean = 80, rho = 0.7) {
  prov.size <- pmax(c(rpois(m, prov.size.mean)), 11)
  N <- sum(prov.size) # total number of discharges
  Z.char <- paste0('z', 1:n.beta)
  sd.err <- 1
  beta <- round(MASS::mvrnorm(n = 1, mu = matrix(runif(1, -1, 1), nrow = n.beta), 
                              Sigma = diag(1, n.beta)), digits = 2)
  gamma <- rnorm(m, 0, sqrt(0.5))
  gamma.rep <- rep(gamma, prov.size)
  prov <- rep(1:m, times = prov.size) 
  rZ <- function(i, rho, n.beta)  
    MASS::mvrnorm(n = prov.size[i], mu = (gamma[i] * rho / sqrt(0.5)) * matrix(1, nrow = n.beta),
                  Sigma = diag(1 - rho, n.beta) + (rho - rho^2) * matrix(1, ncol = n.beta, nrow = n.beta))
  Z <- do.call(rbind, lapply(1:m, FUN = function(i) rZ(i, rho, n.beta))) 
  Y <- gamma.rep + Z %*% beta + rnorm(N, 0, sd.err)
  data <- as.data.frame(cbind(Y, prov, Z))
  colnames(data) <- c("Y", "prov.ID", Z.char)
  result <- list(data = data,
                 true.beta = beta,
                 true.gamma = gamma,
                 prov.char = "prov.ID",
                 Z.char = Z.char,
                 Y.char = "Y")
  return(result)
}


##################---------BINARY GROUP LASSO SETTINGS----------##################
# "ls": list(m, n.beta, n.groups, prop.NonZero.group, prop.outlier, rho)
#                m: number of providers
#           n.beta: number of risk factors
#         n.groups: number of penalized groups
#   prop.NonZero.group: proportion of non-zero coefficients penalized group
#     prop.outlier: proportion of outlying providers
#              rho: corr(Z, gamma) = rho (note: also equals to corr(Z_p1, Z_p2), for all p1, p2 \in {1,2, ..., n.beta} )

#' Simulate clustered data with a binary outcome and grouped covariates
#'
#' Generates the logistic model `logit P(Y = 1) = gamma_k + z' beta`, with one
#' provider effect `gamma_k` per provider and the covariates arranged in groups.
#'
#' @param ls a list with elements `m` (providers), `n.beta` (covariates),
#'   `n.groups` (penalised groups), `prop.NonZero.group` (proportion of groups
#'   with non-zero coefficients), `prop.outlier` (proportion of providers whose
#'   effect is drawn from an outlying distribution) and `rho`.
#' @param prov.size.mean mean provider size.
#' @param unpenalized.beta whether some covariates should be left unpenalised
#'   (assigned to group 0).
#' @param prop.unpenalized.beta the proportion left unpenalised when
#'   `unpenalized.beta = TRUE`.
#'
#' @return A list with `sim.data`, and the `beta`, `gamma` and `group` used to
#'   generate it.
#'
#' @section Note on column order:
#' The covariate COLUMN POSITIONS are shuffled, and the returned `beta` and
#' `group` are aligned with the new positions rather than with the variable
#' names. Take the covariate names from `colnames(sim.data)` rather than
#' constructing them, or `beta` and the fitted coefficients will not correspond.
#'
#' @examples
#' set.seed(2)
#' sim <- Simulation_data_GroupLasso(
#'   list(m = 10, n.beta = 6, n.groups = 3, prop.NonZero.group = 0.5,
#'        prop.outlier = 0, rho = 0.5),
#'   prov.size.mean = 30)
#' # the covariate columns are returned in a shuffled order, so read their
#' # names off the data rather than assuming Z_1 ... Z_p
#' Z.char <- colnames(sim$sim.data)[3:8]
#' fit <- grp.lasso(sim$sim.data, "Y", Z.char, "Prov.ID",
#'                  group = sim$group, nlambda = 10)
#' dim(fit$beta)
#'
#' @export
Simulation_data_GroupLasso <- function(ls, prov.size.mean = 80, unpenalized.beta = FALSE, 
                                       prop.unpenalized.beta = 0.25) {
  m <- ls$m  # number of providers
  n <- ls$n.beta  # number of risk factors
  n.groups <- ls$n.groups  # number of groups
  n.outlier <- floor(ls$prop.outlier * m) # proportion of outlying providers (some gamma's are not generated from normal distribution)
  rho <- ls$rho
  prov.size <- pmax(c(rpois(m, prov.size.mean)), 11) # mean 80 based on real data; 11: provider size cutoff
  gamma <- rnorm(m, log(4/11), 0.4) # log odds; mean and sd based on real data
  gamma[sample.int(m, n.outlier)] <- log(4/11) +
    (2 * rbinom(n.outlier, 1, 0.5) - 1) *
    rnorm(n.outlier, mean = 4 * 0.4, sd = 0.5 * 0.4)  # some gamma's are "outliers"

  n.NonZero.groups <- ceiling(n.groups * ls$prop.NonZero.group)
  beta <- rep(0, n)

  if (unpenalized.beta == FALSE){
    if (n %% n.groups != 0){
      group <- sort(c(rep(1:n.groups, n%/%n.groups), sort(sample(1:n.groups, n%%n.groups))))
    } else {
      group <- sort(c(rep(1:n.groups, n/n.groups)))
    }
    for (i in 1:n.NonZero.groups){
      ind <- which(group == i)
      beta[ind] <- round(MASS::mvrnorm(n = 1, mu = matrix(runif(1, -1, 1), nrow = length(ind)), Sigma = diag(1, length(ind))), digits = 3)
    }
  } else {
    n.penalized <- n - ceiling(n * prop.unpenalized.beta) #number of penalized beta
    if (n.penalized %% n.groups != 0){
      group <- sort(c(rep(0, n - n.penalized), rep(1:(n.groups), n.penalized%/%n.groups), sample(1:n.groups, n.penalized%%(n.groups))))
    } else {
      group <- sort(c(rep(0, n - n.penalized), rep(1:(n.groups), n.penalized%/%n.groups)))
    }
    zero.group.ind <- which(group == 0)
    beta[zero.group.ind] <- round(MASS::mvrnorm(n = 1, mu = matrix(0, nrow = length(zero.group.ind)), Sigma = diag(1, length(zero.group.ind))), digits = 3)
    for (i in 1:n.NonZero.groups){
      ind <- which(group == i)
      beta[ind] <- round(MASS::mvrnorm(n = 1, mu = matrix(runif(1, -1, 1), nrow = length(ind)), Sigma = diag(1, length(ind))), digits = 3)
    }
  }

  Y.char <- 'Y'
  prov.char <- 'Prov.ID'
  Z.char <- paste0('Z_', 1:n)
  
  sim.f <- function(m, prov.size, gamma, beta, Y.char, Z.char, prov.char, rho) {
    N <- sum(prov.size) # total number of discharges
    n.beta <- length(beta)
    gamma.dis <- rep(gamma, prov.size)
    prov <- rep(1:m, prov.size) # provider IDs
    rZ <- function(i, rho, n.beta)  #conditional distribution of Z_ij given gamma_i
      MASS::mvrnorm(n = prov.size[i], mu = ((gamma[i] - log(4/11)) * rho / 0.4) * matrix(1, nrow = n.beta),
                    Sigma = diag(1 - rho, n.beta) + (rho - rho^2) * matrix(1, ncol = n.beta, nrow = n.beta))
    Z <- do.call(rbind, lapply(1:m, FUN = function(i) rZ(i, rho, n.beta)))  # n.beta: length of Z vector
    mu <- plogis(gamma.dis + Z %*% beta)
    Y <- rbinom(N, 1, mu)
    data <- as.data.frame(cbind(Y, prov, Z, mu))
    colnames(data) <- c(Y.char, prov.char, Z.char, "mu")
    return(data)
  }
  
  data <- sim.f(m, prov.size, gamma, beta, Y.char, Z.char, prov.char, rho)
  data <- data[sample(1:nrow(data)),]
  random.ord <- sample(n)
  group <- group[random.ord]
  data[, 3:(2 + n)] <- data[, random.ord + 2]
  colnames(data)[3:(2 + n)] <- colnames(data)[3:(2 + n)][random.ord]
  beta <- beta[random.ord]
  Sim.result <- list(sim.data = data, beta = beta, gamma = gamma, group = group)
  return(Sim.result)
}


##################---------Discrete Survival Model----------##################
# ls <- list(n.center = 1000, n.beta = 5, n.time.point = 10, n.groups = 5, prop.outlier.center = 0, prop.NonZero.group = 1)
# data <- sim.disc(ls, censor_max_t = 10, non_integer_time = FALSE)$data
#' Simulate discrete-time survival data with centre effects
#'
#' Generates a discrete-time logistic hazard model with a separate effect for
#' each centre and a Weibull-shaped baseline hazard.
#'
#' @param ls a list with elements `n.center`, `n.beta`, `n.time.point`,
#'   `n.groups`, `prop.outlier.center` and `prop.NonZero.group`.
#' @param censor_max_t the largest possible censoring time.
#' @param prop.continuous proportion of covariates that are continuous; the rest
#'   are Bernoulli(0.5).
#' @param prov.size.mean mean centre size.
#' @param unpenalized.beta,prop.unpenalized.beta as in
#'   [Simulation_data_GroupLasso()].
#' @param rho correlation between the covariates and the centre effect.
#' @param baseline.Hazard.r,baseline.Hazard.lambda shape and scale of the
#'   Weibull baseline hazard.
#' @param non_integer_time whether to map the integer event times onto
#'   non-integer values.
#'
#' @return A list with `data`, the generating `beta`, `group`, `alpha` and
#'   `gamma`, and the column names.
#'
#' @examples
#' set.seed(3)
#' sim <- sim.disc(
#'   list(n.center = 10, n.beta = 6, n.time.point = 5, n.groups = 3,
#'        prop.outlier.center = 0, prop.NonZero.group = 0.5),
#'   censor_max_t = 5, prov.size.mean = 30)
#' str(sim$data)
#' table(sim$data[[sim$Event.char]])
#'
#' @export
sim.disc <- function(ls, censor_max_t, prop.continuous = 0.8, prov.size.mean = 80, 
                     unpenalized.beta = FALSE, prop.unpenalized.beta = 0, rho = 0.8, 
                     baseline.Hazard.r = 1.2, baseline.Hazard.lambda = 0.01,
                     non_integer_time = TRUE) {
  K <- ls$n.center
  n.groups <- ls$n.groups
  n.outlier.center <- floor(ls$prop.outlier.center * K)
  prov.size <- pmax(c(rpois(K, prov.size.mean)), 11)

  # "alpha" is a vector that contains logit transformation of the baseline hazard at different days
  baseline.hazard <- baseline.Hazard.lambda * baseline.Hazard.r * (baseline.Hazard.lambda * 1:ls$n.time.point)^(baseline.Hazard.r - 1)
  alpha <- log(baseline.hazard)/(1 - baseline.hazard)

  # "gamma" denotes the center effect
  gamma <- rnorm(K, 0, 1)
  gamma[sample.int(K, n.outlier.center)] <- (2 * rbinom(n.outlier.center, 1, 0.5) - 1) * rnorm(n.outlier.center, mean = 4 * 0.4, sd = 0.5 * 0.4)

  n.NonZero.groups <- ceiling(ls$n.groups * ls$prop.NonZero.group)
  beta <- rep(0, ls$n.beta)
  if (unpenalized.beta == FALSE){
    if (ls$n.beta %% ls$n.groups != 0){
      group <- sort(c(rep(1:ls$n.groups, ls$n.beta%/%ls$n.groups), sort(sample(1:ls$n.groups, ls$n.beta%%ls$n.groups))))
    } else {
      group <- sort(c(rep(1:ls$n.groups, ls$n.beta/ls$n.groups)))
    }
    for (i in 1:n.NonZero.groups){
      ind <- which(group == i)
      beta[ind] <- round(MASS::mvrnorm(n = 1, mu = matrix(runif(1, -1, 1), nrow = length(ind)), Sigma = diag(1, length(ind))), digits = 3)
    }
  } else {
    n.penalized <- ls$n.beta  - ceiling(ls$n.beta * prop.unpenalized.beta) #number of penalized beta
    if (n.penalized %% ls$n.groups != 0){
      group <- sort(c(rep(0, ls$n.beta - n.penalized), rep(1:(ls$n.groups), n.penalized%/%ls$n.groups), sample(1:ls$n.groups, n.penalized%%(ls$n.groups))))
    } else {
      group <- sort(c(rep(0, ls$n.beta - n.penalized), rep(1:(ls$n.groups), n.penalized%/%ls$n.groups)))
    }
    zero.group.ind <- which(group == 0)
    beta[zero.group.ind] <- round(MASS::mvrnorm(n = 1, mu = matrix(0, nrow = length(zero.group.ind)), Sigma = diag(1, length(zero.group.ind))), digits = 3)
    for (i in 1:n.NonZero.groups){
      ind <- which(group == i)
      beta[ind] <- round(MASS::mvrnorm(n = 1, mu = matrix(runif(1, -1, 1), nrow = length(ind)), Sigma = diag(1, length(ind))), digits = 3)
    }
  }

  prov.char <- 'Prov.ID'
  Z.char <- paste0('Z', 1:ls$n.beta)

  N <- sum(prov.size) # total number of subjects
  gamma.dis <- rep(gamma, prov.size)
  prov <- rep(1:K, prov.size) # provider IDs
  rZ <- function(i, rho, n.beta) { #conditional distribution of Z_ij given alpha_i
    MASS::mvrnorm(n = prov.size[i], mu = ((gamma[i] - 0) * rho / 1) * matrix(1, nrow = n.beta),
                  Sigma = diag(1 - rho, n.beta) + (rho - rho^2) * matrix(1, ncol = n.beta, nrow = n.beta))
  }  

  n.beta.continuous <- floor(prop.continuous * ls$n.beta) # number of continuous variable
  n.beta.binary <- ls$n.beta - n.beta.continuous

  Z.continuous <- do.call(rbind, lapply(1:K, FUN = function(i) rZ(i, rho, n.beta.continuous)))  # n.beta: length of Z vector
  Z.binary <- matrix(rbinom(N * n.beta.binary, 1, 0.5), N, n.beta.binary) # binary covariate, p = 0.5
  Z <- cbind(Z.continuous, Z.binary)

  day <- 1 #current day
  idx.atrisk <- 1:N  # current at risk people
  days.to.event <- rep(length(alpha), N) # initial days to event: maximum date
  status <- rep(0, N) # 0:censor; 1:failure
  probs <- plogis(alpha[1] + (gamma.dis + as.matrix(Z) %*% beta)) #probability of failure at t = 1
  idx.event <- idx.atrisk[rbinom(length(probs), 1, probs) == 1] # index of people who fail at t = 1
  status[idx.event] <- 1 # set the event of people who fail at t = 1 to "1"
  days.to.event[idx.event] <- day  # if failure, set their observed time
  idx.out <- idx.event  # who is not at risk
  censoring <- runif(N, 1, censor_max_t)  # randomly specified "potential" censoring time of each observations. Follows Uniform distribution
  conTime <- data.frame(time = censoring)  #everyone's censor time
  # "contToDisc": Continuous to Discrete Transformation. This function can transfer continuous "time" to integer (ceiling)
  censoring_time <- as.numeric(discSurv::contToDisc(dataShort = conTime, timeColumn = "time", intervalLimits = 1:(censor_max_t))$timeDisc)
  for (x in tail(alpha,-1)) { # tail(alpha, -1) exclude the first element of time effect
    day <- day + 1 # move to the next day
    idx.atrisk <- c(1:N)[-idx.out]  # remove people already failure
    probs <- plogis(x + (gamma.dis[idx.atrisk] + as.matrix(Z[idx.atrisk, ]) %*% beta)) # x now is the baseline hazard at "day"'s time
    idx.event <- idx.atrisk[rbinom(length(probs), 1, probs) == 1]  #some at risk people have large hazard, such that failure
    status[idx.event] <- 1 # new failure index
    days.to.event[idx.event] <- day #which day failure
    idx.out <- unique(c(idx.out, idx.event))  #why "unique"
  }

  tcens <- as.numeric(censoring < days.to.event) # censoring indicator
  delta <- 1-tcens #failure = 1, censor = 0
  time <- days.to.event * (delta == 1) + censoring_time * (delta == 0)  # final time consider cencor and failure
  delta[-idx.out] <- 0 #those have no event at last date also censor
  data <- as.data.frame(cbind(prov, Z, time, delta))
  colnames(data) <- c(prov.char, Z.char, "time", "status")
  char <- c(prov.char, Z.char, "time", "status")
  
  # change "integer time" into "non-integer"
  if (non_integer_time == TRUE){
    non_int_time <- sort(round(runif(ls$n.time.point, 0, 2 * ls$n.time.point), 2))
    time.ref <- as.data.frame(cbind(sort(unique(data$time)), non_int_time[1:length(unique(data$time))]))
    colnames(time.ref) <- c("time", "non_int_time")
    data <- plyr::join(data, time.ref, by = c("time"))
    data <- data[, !(colnames(data) %in% c("time"))]
    colnames(data)[ncol(data)] <- "time"
  }
  
  result <- list(data = data,
                 beta = beta,
                 group = group,
                 alpha = alpha,
                 gamma = gamma,
                 prov.char = prov.char,
                 Z.char = Z.char,
                 Time.char = "time",
                 Event.char = "status")
  return(result)
}

##################---------Cox model----------##################
#sim <- sim.cox(n.center = 20, n.beta = 5)
#data <- sim$data
#' Simulate right-censored survival data with provider effects
#'
#' Generates survival times from a proportional hazards model in which the
#' hazard of provider k is `gamma_k * exp(z' beta)`, so `gamma_k` acts as a
#' provider-specific baseline hazard.
#'
#' @param n.center number of providers.
#' @param n.beta number of covariates.
#' @param prop.zero_beta the proportion of covariates with a NON-zero
#'   coefficient. The name is historical and does not match the behaviour.
#' @param prov.size.mean mean provider size.
#' @param unpenalized.beta,prop.unpenalized.beta currently unused.
#' @param rho correlation between the covariates and the provider effect.
#'
#' @return A list with `data`, the generating `beta`, and the column names
#'   `prov.char`, `Z.char`, `Time.char`, `Event.char`.
#'
#' @section Known issue:
#' The covariates are centred using the mean and standard deviation of the
#' provider effect in [Simulation_data_GroupLasso()] (`log(4/11)` and `0.4`)
#' rather than those of this function's own `gamma ~ Uniform(0, 2)`. The
#' covariate mean is therefore displaced by roughly `2*rho` to `6*rho`, the
#' linear predictor inherits a large offset, and the hazard either explodes or
#' vanishes: measured event rates range from 0.000 to 0.998 depending on `rho`
#' and `n.beta`. Standardise `gamma` with its own moments (mean 1, standard
#' deviation `sqrt(1/3)`) before relying on this function.
#'
#' @examples
#' set.seed(4)
#' sim <- sim.cox(n.center = 10, n.beta = 6, prov.size.mean = 30)
#' d <- as.data.frame(sim$data)
#' # NOTE the event rate: see the "Known issue" section above.  This example
#' # documents the behaviour as shipped, it does not endorse it.
#' mean(d[[sim$Event.char]])
#' fit <- Strat.cox(d, Event.char = sim$Event.char, Z.char = sim$Z.char,
#'                  Time.char = sim$Time.char, prov.char = sim$prov.char,
#'                  nlambda = 5)
#' dim(fit$beta)
#'
#' @export
sim.cox <- function(n.center, n.beta, prop.zero_beta = 0.1, prov.size.mean = 80,
                    unpenalized.beta = FALSE, prop.unpenalized.beta = 0, rho = 0.8) {
  
  prov.size <- pmax(c(rpois(n.center, prov.size.mean)), 11)
  gamma <- runif(n.center, 0, 2)
  gamma_subject <- rep(gamma, prov.size)
  N <- sum(prov.size)
  
  prov <- rep(1:n.center, prov.size) # provider IDs
  
  beta <- rep(0, n.beta)
  n.NonZero_beta <- ceiling(n.beta * prop.zero_beta)
  beta[1:n.NonZero_beta] <- round(MASS::mvrnorm(n = 1, mu = matrix(runif(1, -1, 1), nrow = n.NonZero_beta),
                                                Sigma = diag(1, n.NonZero_beta)), digits = 2)

  rZ <- function(i, rho, n.beta){
    Z <- MASS::mvrnorm(n = prov.size[i], mu = ((gamma[i] - log(4/11)) * rho / 0.4) * matrix(1, nrow = n.beta),
                       Sigma = diag(1 - rho, n.beta) + (rho - rho^2) * matrix(1, ncol = n.beta, nrow = n.beta))
    return(Z)
  }  
  Z <- do.call(rbind, lapply(1:n.center, FUN = function(i) rZ(i, rho, n.beta)))  
  
  U <- runif(N, 0, 1)
  pre_time <- -log(U) / (gamma_subject * exp(Z %*% beta))
  
  pre_censoring <- runif(N, 1, 10)
  pre_censoring <- pre_censoring * (pre_censoring < 3) + 3 * (pre_censoring >= 3)
  tcens <- (pre_censoring < pre_time) # censoring indicator
  delta <- 1 - tcens
  time <- pre_time * (delta == 1) + pre_censoring * (delta == 0)
  
  data <- cbind(prov, Z, delta, time)
  colnames(data) <- c("Prov.ID", paste0("Z", 1:n.beta), "status", "time")
  data <- data[sample(nrow(data)),]
  
  Sim.result <- list(data = data, 
                     beta = beta,
                     prov.char = "Prov.ID",
                     Z.char = paste0("Z", 1:n.beta),
                     Time.char = "time",
                     Event.char = "status")
  return(Sim.result)
}


