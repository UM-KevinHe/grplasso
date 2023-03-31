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


##################---------GROUP LASSO SETTINGS----------##################

# "ls": list(m, n.beta, n.groups, prop.NonZero.group, prop.outlier, rho)
#                m: number of providers
#           n.beta: number of risk factors
#         n.groups: number of penalized groups
#   prop.NonZero.group: proportion of non-zero coefficients penalized group
#     prop.outlier: proportion of outlying providers
#              rho: corr(Z, gamma) = rho (note: also equals to corr(Z_p1, Z_p2), for all p1, p2 \in {1,2, ..., n.beta} )


Simulation_data_GroupLasso <- function(ls, prov.size.mean = 80, unpenalized.beta = F, prop.unpenalized.beta = 0.25) {
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


##################---------Discrete Survival----------##################

#ls <- list(n.center = 10, n.beta = 10, n.time.point = 10, n.groups = 10, prop.outlier.center = 0, prop.NonZero.group = 1)

sim.disc <- function(ls, censor_max_t, prop.continuous = 0.8,
                     prov.size.mean = 80, unpenalized.beta = F, prop.unpenalized.beta = 0,
                     rho = 0.8, baseline.Hazard.r = 1.2, baseline.Hazard.lambda = 0.01) {
  K <- ls$n.center
  n.groups <- ls$n.groups
  n.outlier.center <- floor(ls$prop.outlier.center * K)
  prov.size <- pmax(c(rpois(K, prov.size.mean)), 11)


  # "gamma" is a vector that contains logit transformation of the baseline hazard at different days
  baseline.hazard <- baseline.Hazard.lambda * baseline.Hazard.r * (baseline.Hazard.lambda * 1:ls$n.time.point)^(baseline.Hazard.r - 1)
  gamma <- log(baseline.hazard)/(1 - baseline.hazard)


  # "alpha" denotes the center effect
  alpha <- rnorm(K, 0, 1)
  alpha[sample.int(K, n.outlier.center)] <- (2 * rbinom(n.outlier.center, 1, 0.5) - 1) * rnorm(n.outlier.center, mean = 4 * 0.4, sd = 0.5 * 0.4)

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
  alpha.dis <- rep(alpha, prov.size)
  prov <- rep(1:K, prov.size) # provider IDs
  rZ <- function(i, rho, n.beta)  #conditional distribution of Z_ij given gamma_i
    MASS::mvrnorm(n = prov.size[i], mu = ((alpha[i] - 0) * rho / 1) * matrix(1, nrow = n.beta),
                  Sigma = diag(1 - rho, n.beta) + (rho - rho^2) * matrix(1, ncol = n.beta, nrow = n.beta))

  n.beta.continuous <- floor(prop.continuous * ls$n.beta) # number of continuous variable
  n.beta.binary <- ls$n.beta - n.beta.continuous

  Z.continuous <- do.call(rbind, lapply(1:K, FUN = function(i) rZ(i, rho, n.beta.continuous)))  # n.beta: length of Z vector
  Z.binary <- matrix(rbinom(N * n.beta.binary, 1, 0.5), N, n.beta.binary) # binary covariate, p = 0.5
  Z <- cbind(Z.continuous, Z.binary)

  day <- 1 #current day
  idx.atrisk <- 1:N  # current at risk people
  days.to.event <- rep(length(gamma), N) # initial days to event: maximum date
  status <- rep(0, N) # 0:censor; 1:failure
  probs <- plogis(gamma[1] + (alpha.dis + as.matrix(Z) %*% beta)) #probability of failure at t = 1
  idx.event <- idx.atrisk[rbinom(length(probs), 1, probs) == 1] # index of people who fail at t = 1
  status[idx.event] <- 1 # set the event of people who fail at t = 1 to "1"
  days.to.event[idx.event] <- day  # if failure, set their observed time
  idx.out <- idx.event  # who is not at risk
  censoring <- runif(N, 1, censor_max_t)  # randomly specified "potential" censoring time of each observations. Follows Unif(1, 100)
  conTime <- data.frame(time = censoring)  #everyone's censor time
  # "contToDisc": Continuous to Discrete Transformation. This function can transfer continuous "time" to integer (ceiling)
  censoring_time <- as.numeric(discSurv::contToDisc(dataShort = conTime, timeColumn = "time", intervalLimits = 1:(censor_max_t))$timeDisc)#3
  for (x in tail(gamma,-1)) { # tail(gamma, -1) exclude the first element of time effect
    day <- day + 1 # move to the next day
    idx.atrisk <- c(1:N)[-idx.out]  # remove people already failure
    probs <- plogis(x + (alpha.dis[idx.atrisk] + as.matrix(Z[idx.atrisk, ]) %*% beta)) # x now is the baseline hazard at "day"'s time
    idx.event <- idx.atrisk[rbinom(length(probs), 1, probs) == 1]  #some at risk people have large hazard, such that failure
    status[idx.event] <- 1 # new failure index
    days.to.event[idx.event] <- day #which day failure
    idx.out <- unique(c(idx.out, idx.event))  #why "unique"
  }

  tcens <- as.numeric(censoring<days.to.event) # censoring indicator
  delta <- 1-tcens #failure = 1, censor = 0
  time <- days.to.event * (delta == 1) + censoring_time * (delta == 0)  # final time consider cencor and failure
  delta[-idx.out] <- 0 #those have no event at last date also censor
  data <- as.data.frame(cbind(prov, Z, time, delta))
  colnames(data) <- c(prov.char, Z.char, "time", "status")
  char <- c(prov.char, Z.char, "time", "status")
  result <- list(data = data,
                 beta = beta,
                 group = group,
                 gamma =gamma,
                 alpha = alpha,
                 prov.char = prov.char,
                 Z.char = Z.char,
                 Time.char = "time",
                 Event.char = "status")
  return(result)
}
