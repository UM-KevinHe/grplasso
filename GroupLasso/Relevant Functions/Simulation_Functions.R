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


Simulation_data_GroupLasso <- function(ls, unpenalized.beta = F, prop.unpenalized.beta = 0.25) { 
  m <- ls$m  # number of providers
  n <- ls$n.beta  # number of risk factors
  n.groups <- ls$n.groups  # number of groups
  n.outlier <- floor(ls$prop.outlier * m) # proportion of outlying providers (some gamma's are not generated from normal distribution)
  rho <- ls$rho 
  prov.size <- pmax(c(rpois(m, 5000)), 11) # mean 80 based on real data; 11: provider size cutoff
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

