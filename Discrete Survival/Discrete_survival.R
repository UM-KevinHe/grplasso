AR1 <- function(tau, m) {
  if(m==1) {R <- 1}
  if(m > 1) {
    R <- diag(1, m)
    for(i in 1:(m-1)) {
      for(j in (i+1):m) {
        R[i,j] <- R[j,i] <- tau^(abs(i-j))
      }
    }
  }
  return(R)
}

simu_z <- function(n, size.groups)
{
  Sigma_z1 <- diag(size.groups) # =p
  Corr1 <- AR1(0.5, size.groups) #correlation structure 0.5 (AR structure)
  diag(Corr1) <- 1
  Sigma_z1 <- Corr1
  pre_z <- rmvnorm(n, mean = rep(0, size.groups), sigma = Sigma_z1) #mean = 0, cov matrix follows AR structure; can be check by cor(pre_Z)
  return(pre_z)
}

##  > Sigma_z1
##       [,1]  [,2] [,3]
##  [1,] 1.00  0.5  0.25
##  [2,] 0.50  1.0  0.50
##  [3,] 0.25  0.5  1.00


sim.disc <- function(beta, day_effect, n.obs, Z.char, censor_max_t, prop.continuous = 0.8) {
  # "day_effect" is a vector that contains logit transformation of the baseline hazard at different days (gamma)
  N <- n.obs 
  p_1 <- floor(prop.continuous * length(beta)) # number of continuous variable
  p_2 <- length(beta) - p_1
  Z1 <- as.matrix(simu_z(N, p_1)) #continuous covariate, AR correlation structure
  Z2 <- matrix(rbinom(N * p_2, 1, 0.5), N, p_2) # binary covariate, p = 0.5
  Z <- cbind(Z1, Z2) # covariate matrix
  day <- 1 #current day
  idx.atrisk <- 1:N  # current at risk people
  days.to.event <- rep(length(day_effect), N) # initial days to event: maximum date
  status <- rep(0, N) # 0:cencor; 1:failure
  probs <- plogis(day_effect[1] + (as.matrix(Z) %*% beta)) #probability of failure at t = 1
  idx.event <- idx.atrisk[rbinom(length(probs), 1, probs) == 1] # index of people who fail at t = 1
  status[idx.event] <- 1 # set the event of people who fail at t = 1 to "1"
  days.to.event[idx.event] <- day  # if failure, set their observed time
  idx.out <- idx.event  # who is not at risk
  censoring <- runif(N, 1, censor_max_t)  # randomly specified "potential" censoring time of each observations. Follows Unif(1, 100)
  conTime <- data.frame(time = censoring)  #everyone's censor time
  # "contToDisc": Continuous to Discrete Transformation. This function can transfer continuous "time" to integer (ceiling)
  censoring_time <- as.numeric(contToDisc(dataShort = conTime, timeColumn = "time", intervalLimits = 1:(censor_max_t))$timeDisc)#3
  for (x in tail(day_effect,-1)) { # tail(day_effect, -1) exclude the first element of day_effect
    day <- day + 1 # move to the next day
    idx.atrisk <- c(1:N)[-idx.out]  # remove people already failure   
    probs <- plogis(x + (as.matrix(Z[idx.atrisk, ]) %*% beta)) # x now is the baseline hazard at "day"'s time
    idx.event <- idx.atrisk[rbinom(length(probs), 1, probs) == 1]  #some at risk people have large hazard, such that failure
    status[idx.event] <- 1 # new failure index
    days.to.event[idx.event] <- day #which day failure
    idx.out <- unique(c(idx.out, idx.event))  #why "unique"
  }
  
  tcens <- as.numeric(censoring<days.to.event) # censoring indicator
  delta <- 1-tcens #failure = 1, censor = 0
  time <- days.to.event * (delta == 1) + censoring_time * (delta == 0)  # final time consider cencor and failure
  delta[-idx.out] <- 0 #those have no event at last date also censor
  data <- as.data.frame(cbind(delta, Z, time))
  colnames(data) <- c("status", Z.char, "time")
  return(data)
}


# main function
discreSurv_logit <- function(t, X, ind, tol = 1e-4, max_iter = 1000){
  maxt <- max(t)
  c <- ncol(X)
  r <- nrow(X)
  od <- order(t, decreasing = T)
  t <- t[od]
  ind <- ind[od]
  X <- X[od, ]
  ind <- as.matrix(ind)
  beta_t <- as.matrix(rep(0, maxt))
  beta_v <- as.matrix(rep(0, c))
  t <- as.matrix(t)
  X <- as.matrix(X)
  return (NR_residuals(t, X, ind, beta_t, beta_v, tol, max_iter))
}




