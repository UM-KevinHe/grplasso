library(MASS)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(grid)
library(dplyr)

setwd("/home/ybshao/GroupLasso/Simulations/Figures&Tables/FE_vs_Pooled/GLM_Binary/FDP/Table")

# ls <- list(m = 40, poisson.mean = 20, rho = 0.6, n.rep = 50, n.beta = 100, prop.NonZero.beta = 0.05, beta = beta)
sim.FDP <- function(ls){
  print(paste0("Iteration for beta = ", ls$beta, " starts.."))
  start.time <- Sys.time()
  sd.gamma <- 1
  m <- ls$m
  sd.z <- 1
  true.beta <- c(rep(ls$beta, ls$n.beta * ls$prop.NonZero.beta), 
                 rep(0, ls$n.beta * (1 - ls$prop.NonZero.beta)))  # 10% beta's are non-zero
  loop_length <- ls$n.rep
  rho <- ls$rho
  sim.continuous <- function() {
    data.proveff3 <- function(m, prov.size, gamma, sd.gamma, sd.z,
                              beta, Y.char, Z.char, prov.char, rho) {
      N <- sum(prov.size) # total number of discharges
      n.beta <- length(beta)
      gamma.rep <- rep(gamma, prov.size)
      prov <- rep(paste0("prov", sprintf(paste0("%0",nchar(m),"d"),1:m)), prov.size) # provider IDs
      KW2013 <- function(i, rho, n.beta){
        MASS::mvrnorm(n = prov.size[i], mu = (sd.z * gamma[i] * rho / sd.gamma) * matrix(1, nrow = n.beta),
                      Sigma = sd.z^2 * (diag(1 - rho, n.beta) + (rho - rho^2) * matrix(1, ncol = n.beta, nrow = n.beta)))
      }
      Z <- do.call(rbind, lapply(1:m, function(i) KW2013(i, rho, n.beta)))
      Y <- rbinom(N, 1, plogis(gamma.rep + Z %*% beta))
      data <- data.frame(Y, prov, Z, stringsAsFactors = F)
      colnames(data) <- c(Y.char, prov.char, Z.char)
      return(data)
    }
    prov.size <- pmax(c(rpois(m, ls$poisson.mean)), 20)
    gamma <- rnorm(m, 0, sd.gamma)
    Y.char <- 'Y'
    prov.char <- 'unit'
    Z.char <- paste0('z', 1:length(true.beta))
    tb <- data.proveff3(m, prov.size, gamma, sd.gamma, sd.z, true.beta, Y.char, Z.char, prov.char, rho)
    
    true.non.zero.beta.index <- 1:(ls$n.beta * ls$prop.NonZero.beta)
    true.zero.beta.index <- (ls$n.beta * ls$prop.NonZero.beta + 1):ls$n.beta
    
    # (1) Fixed effect model 
    cv.FE.lasso <- cv.pp.lasso(tb, Y.char, Z.char, prov.char, nfolds = 5) 
    cv.FE.fit <- cv.FE.lasso$fit
    FE.beta.under.best.lambda <- cv.FE.fit$beta[, cv.FE.lasso$min]  #beta estimate under best lambda (by 10-fold cross validation)
    FE.estimated.non.zero.beta.index <- which(FE.beta.under.best.lambda != 0)
    FE.estimated.zero.beta.index <- which(FE.beta.under.best.lambda == 0)
    
    FE.FD <- sum(!(FE.estimated.non.zero.beta.index %in% true.non.zero.beta.index))
    FE.FN <- sum(!(FE.estimated.zero.beta.index %in% true.zero.beta.index))
    if (length(FE.estimated.non.zero.beta.index) == 0) {
      FE.FDP <- 0
    } else {
      FE.FDP <- FE.FD/length(FE.estimated.non.zero.beta.index)
    }
    FE.FNP <- FE.FN/length(FE.estimated.zero.beta.index)
    
    
    # (2) Pooled model 
    cv.pooled.lasso <- cv.pp.lasso(tb, Y.char, Z.char, nfolds = 5)
    cv.pooled.fit <- cv.pooled.lasso$fit
    pooled.beta.under.best.lambda <- cv.pooled.fit$beta[, cv.pooled.lasso$min]  #beta estimate under best lambda (by 10-fold cross validation)
    pooled.estimated.non.zero.beta.index <- which(pooled.beta.under.best.lambda != 0)
    pooled.estimated.zero.beta.index <- which(pooled.beta.under.best.lambda == 0)
    
    pooled.FD <- sum(!(pooled.estimated.non.zero.beta.index %in% true.non.zero.beta.index))
    pooled.FN <- sum(!(pooled.estimated.zero.beta.index %in% true.zero.beta.index))
    if (length(pooled.estimated.non.zero.beta.index) == 0){
      pooled.FDP <- 0
    } else {
      pooled.FDP <- pooled.FD/length(pooled.estimated.non.zero.beta.index)
    }
    pooled.FNP <- pooled.FN/length(pooled.estimated.zero.beta.index)
    
    res <- c(FE.FD, FE.FN, FE.FDP, FE.FNP, pooled.FD, pooled.FN, pooled.FDP, pooled.FNP) 
    return(res)
  }
  
  `%dopar%` <- foreach::`%dopar%`
  cl <- parallel::makeCluster(5)
  doParallel::registerDoParallel(cl)
  system.time(
    mat <- foreach::foreach(i = 1:loop_length, .combine = cbind, .packages = c("TmpLasso", "Matrix")) %dopar% {
      sim.continuous()})
  parallel::stopCluster(cl)
  
  #get statistics after loops end
  Mean.Performance <- rowMeans(mat, na.rm = T)
  Sd.Performance <- apply(mat, 1, sd, na.rm = T)
  
  rst <- c(ls$m, ls$poisson.mean, rho, ls$beta, Mean.Performance, Sd.Performance) 
  names(rst) <- c("num.units", "possion.mean", "rho", "current.beta",
                  "FE.FD.Mean", "FE.FN.Mean", "FE.FDP.Mean", "FE.FNP.Mean", 
                  "pooled.FD.Mean", "pooled.FN.Mean", "pooled.FDP.Mean", "pooled.FNP.Mean",
                  "FE.FD.Sd", "FE.FN.Sd", "FE.FDP.Sd", "FE.FNP.Sd", 
                  "pooled.FD.Sd", "pooled.FN.Sd", "pooled.FDP.Sd", "pooled.FNP.Sd")
  end.time <- Sys.time()
  process.time <- difftime(end.time, start.time, units = 'mins')
  print(paste0("Processing time: ", round(process.time, digits = 4), " mins"))
  return(rst)
}

sim.FDP.res <- NULL
for (beta in c(2)){
  ls <- list(m = 50, poisson.mean = 25, rho = 0.6, n.rep = 1000, n.beta = 100, prop.NonZero.beta = 0.05, beta = beta)
  sim.FDP.res <- rbind(sim.FDP.res, sim.FDP(ls))
}
sim.FDP.res <- round(as.data.frame(sim.FDP.res), digits = 3)
save(sim.FDP.res, file = paste0("table_2_", Sys.Date(), ".RData"))


sim.FDP.res.cb <- NULL
sim.FDP.res.sub <- sim.FDP.res
sim.FDP.res.cb <- rbind(sim.FDP.res.cb, sim.FDP.res.sub)

sim.FDP.res <- sim.FDP.res.cb
sim.FDP.df <- cbind(sim.FDP.res$current.beta,
                    paste0(sim.FDP.res$FE.FD.Mean, " (", sim.FDP.res$FE.FD.Sd, ") "),
                    paste0(sim.FDP.res$FE.FN.Mean, " (", sim.FDP.res$FE.FN.Sd, ") "),
                    paste0(sim.FDP.res$pooled.FD.Mean, " (", sim.FDP.res$pooled.FD.Sd, ") "),
                    paste0(sim.FDP.res$pooled.FN.Mean, " (", sim.FDP.res$pooled.FN.Sd, ") "),
                    paste0(sim.FDP.res$FE.FDP.Mean, " (", sim.FDP.res$FE.FDP.Sd, ") "),
                    paste0(sim.FDP.res$FE.FNP.Mean, " (", sim.FDP.res$FE.FNP.Sd, ") "),
                    paste0(sim.FDP.res$pooled.FDP.Mean, " (", sim.FDP.res$pooled.FDP.Sd, ") "),
                    paste0(sim.FDP.res$pooled.FNP.Mean, " (", sim.FDP.res$pooled.FNP.Sd, ") "))
colnames(sim.FDP.df) <- c("beta", "FD.FE", "FN.FE", "FD.Pooled", "FN.Pooled", 
                          "FDP.FE", "FNP.FE", "FDP.Pooled", "FNP.Pooled")
save(sim.FDP.df, file = paste0("FE_vs_Pooled_FDP_table_NEW_", Sys.Date(), ".RData"))
write.xlsx(sim.FDP.df, file = paste0("FE_vs_Pooled_FDP_table_", Sys.Date(), ".xlsx"), sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE)

