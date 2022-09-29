setwd("/home/ybshao/GroupLasso/Simulations/Figures&Tables/Accuracy_RME&RMSE")
library(MASS)
library(Matrix)
library(grpreg)
library(fastDummies)
library(foreach)
library(doParallel)
library(glmnet)
library(TmpLasso)
library(parallel)

######----------- GroupLASSO Part ------------######
grouplasso.RMSE <- function(ls, beta.sequence, ind, n.rep){
  beta <- beta.sequence[ind]
  data.loop <- 1:n.rep
  
  multiResultClass <- function(RMSE = NULL, RME = NULL){
    result <- list(
      RMSE = RMSE,
      RME = RME
    )
    class(result) <- append(class(result), "multiResultClass")
    return(result)
  }
  
  cl <- makeCluster(5)
  registerDoParallel(cl) 
  Model.Comparison <- 
    foreach (i = data.loop, .packages = c("grpreg", "fastDummies", "RcppArmadillo", "MASS", "Matrix", "TmpLasso")) %dopar% {
      sd.gamma <- 1
      m <- 50
      poisson.mean <- 100
      sd.z <- 1
      true.beta <- c(rep(beta, ls$n.beta * ls$prop.NonZero.beta), 
                     rep(0, ls$n.beta * (1 - ls$prop.NonZero.beta)))  # 10/50 beta's are non-zero
      rho <- ls$rho
      group <- rep(1:10, each = 5)
      
      data.proveff3 <- function(m, prov.size, gamma, sd.gamma, sd.z,
                                beta, Y.char, Z.char, prov.char, rho) {
        N <- sum(prov.size) # total number of discharges
        n.beta <- length(beta)
        gamma.rep <- rep(gamma, prov.size)
        prov <- rep(1:m, prov.size) # provider IDs
        KW2013 <- function(i, rho, n.beta){
          MASS::mvrnorm(n = prov.size[i], mu = (sd.z * gamma[i] * rho / sd.gamma) * matrix(1, nrow = n.beta),
                        Sigma = sd.z^2 * (diag(1 - rho, n.beta) + (rho - rho^2) * matrix(1, ncol = n.beta, nrow = n.beta)))
        }
        Z <- do.call(rbind, lapply(1:m, function(i) KW2013(i, rho, n.beta)))
        Y <- rbinom(N, 1, plogis(gamma.rep + Z %*% beta))
        mu <- plogis(gamma.rep + Z %*% beta)
        data <- data.frame(Y, prov, Z, mu, stringsAsFactors = F)
        colnames(data) <- c(Y.char, prov.char, Z.char, "true.mu")
        return(data)
      }
      prov.size <- pmax(c(rpois(m, poisson.mean)), 20)
      gamma <- rnorm(m, 0, sd.gamma)
      Y.char <- 'Y'
      prov.char <- 'unit'
      Z.char <- paste0('z', 1:length(true.beta))
      sim_data <- data.proveff3(m, prov.size, gamma, sd.gamma, sd.z, true.beta, Y.char, Z.char, prov.char, rho)
      true.mu <- sim_data$true.mu
      
      # results from grpLasso
      cv.model_grp_lasso <- cv.grp.lasso(sim_data, Y.char, Z.char, prov.char, group = group, trace.lambda = F,
                                         nfolds = 10, trace.cv = F, MM = F)
      cv_BestModel_grp_lasso <- cv.model_grp_lasso$fit
      best.beta.grp_lasso <- cv_BestModel_grp_lasso$beta[, cv.model_grp_lasso$min]
      best.eta.grp_lasso <- cv_BestModel_grp_lasso$linear.predictors[, cv.model_grp_lasso$min]
      
      RMSE.grp_lasso <- sqrt(sum((best.beta.grp_lasso - true.beta) ^ 2) / ls$n.beta) #RMSE
      RME.grp_lasso <- sqrt(sum((plogis(best.eta.grp_lasso) - true.mu) ^ 2) / nrow(sim_data)) #RME
      
      
      
      ## results from grpreg
      dummy_data <- dummy_cols(sim_data, select_columns = prov.char, remove_selected_columns = TRUE, 
                               remove_first_dummy = TRUE)
      ID.char <- rep(NA, m - 1)
      for (i in 1:(m - 1)){
        ID.char[i] <- paste0("unit_", i + 1)
      }
      cv.model_grpreg <- cv.grpreg(dummy_data[,c(Z.char, ID.char)], dummy_data[,Y.char], family = "binomial",
                                   penalty = "grLasso", group = c(group, rep(0, length(ID.char))), alpha = 1, 
                                   nfolds = 10, trace.cv = F)
      cv_BestModel_grpreg <- cv.model_grpreg$fit
      best.beta.grpreg <- cv_BestModel_grpreg$beta[2:(1 + ls$n.beta), cv.model_grpreg$min]
      best.eta.grpreg <- cv_BestModel_grpreg$linear.predictors[, cv.model_grpreg$min]
      
      RMSE.grpreg <- sqrt(sum((best.beta.grpreg - true.beta)^2) / ls$n.beta)  #RMSE
      RME.grpreg <- sqrt(sum((plogis(best.eta.grpreg) - true.mu)^2)/nrow(sim_data))  #RME
      
      result <- multiResultClass()
      result$RMSE <- matrix(c(RMSE.grp_lasso, RMSE.grpreg), nrow = 2)
      result$RME <- matrix(c(RME.grp_lasso, RME.grpreg), nrow = 2)
      return(result)
    }
  stopCluster(cl)
  
  n.data.loop <- length(data.loop)
  
  RME.Mean <- as.data.frame(matrix(rep(0, length(beta.sequence) * 2), nrow = 2))
  colnames(RME.Mean) <- paste0("beta = ", beta.sequence)
  rownames(RME.Mean) <- c("grplasso", "grpreg")
  
  RME.sd <- as.data.frame(matrix(rep(0, length(beta.sequence) * 2), nrow = 2))
  colnames(RME.sd) <- paste0("beta = ", beta.sequence)
  rownames(RME.sd) <- c("grplasso", "grpreg")
  
  RMSE.Mean <- as.data.frame(matrix(rep(0, length(beta.sequence) * 2), nrow = 2))
  colnames(RMSE.Mean) <- paste0("beta = ", beta.sequence)
  rownames(RMSE.Mean) <- c("grplasso", "grpreg")
  
  RMSE.sd <- as.data.frame(matrix(rep(0, length(beta.sequence) * 2), nrow = 2))
  colnames(RMSE.sd) <- paste0("beta = ", beta.sequence)
  rownames(RMSE.sd) <- c("grplasso", "grpreg")
  
  RMSE <- matrix(rep(0, 2 * n.data.loop), nrow = 2)
  rownames(RMSE) <- c("grLasso", "grpreg")
  colnames(RMSE) <- paste0("Data_", data.loop)
  for (i in data.loop){
    RMSE[, i] <- Model.Comparison[[i]]$RMSE
  }
  RMSE.Mean[, ind] <- round(apply(RMSE, 1, mean), digits = 3)
  RMSE.sd[, ind] <- round(apply(RMSE, 1, sd), digits = 3)
  RMSE.df <- RMSE.Mean
  RMSE.df[, ind] <- paste0(RMSE.Mean[, ind], " (", RMSE.sd[, ind], ") ")
  
  
  RME <- matrix(rep(0, 2 * n.data.loop), nrow = 2)
  rownames(RME) <- c("grLasso", "grpreg")
  colnames(RME) <- paste0("Data_", data.loop)
  for (i in data.loop){
    RME[, i] <- Model.Comparison[[i]]$RME
  }
  RME.Mean[, ind] <- round(apply(RME, 1, mean), digits = 3)
  RME.sd[, ind] <- round(apply(RME, 1, sd), digits = 3)
  RME.df <- RME.Mean
  RME.df[, ind] <- paste0(RME.Mean[, ind], " (", RME.sd[, ind], ") ")
  
  return(list(RMSE.df = RMSE.df,
              RME.df = RME.df))
}


ls <- list(rho = 0.6, n.beta = 50, prop.NonZero.beta = 0.2)
beta.sequence <- c(0.5, 1, 2)

### simulation

grplasso.beta05 <- grouplasso.RMSE(ls, beta.sequence, ind = 1, n.rep = 100)
save(grplasso.beta05, file = paste0("grplasso_RME_RMSE_beta05_", Sys.Date(), ".RData"))

grplasso.beta1 <- grouplasso.RMSE(ls, beta.sequence, ind = 2, n.rep = 100)
save(grplasso.beta1, file = paste0("grplasso_RME_RMSE_beta1_", Sys.Date(), ".RData"))

grplasso.beta2 <- grouplasso.RMSE(ls, beta.sequence, ind = 3, n.rep = 100)
save(grplasso.beta2, file = paste0("grplasso_RME_RMSE_beta2_", Sys.Date(), ".RData"))


grplasso.RMSE <- cbind(grplasso.beta05$RMSE.df[, 1, drop = F], grplasso.beta1$RMSE.df[, 2, drop = F], grplasso.beta2$RMSE.df[, 3, drop = F])
grplasso.RME <- cbind(grplasso.beta05$RME.df[, 1, drop = F], grplasso.beta1$RME.df[, 2, drop = F], grplasso.beta2$RME.df[, 3, drop = F])

save(grplasso.RMSE, grplasso.RME, file = paste0("grplasso_RME_RMSE_combine_", Sys.Date(), ".RData"))

######----------- LASSO Part ------------######
lasso.RMSE <- function(ls, beta.sequence, ind, n.rep){
  beta <- beta.sequence[ind]
  data.loop <- 1:n.rep
  
  multiResultClass <- function(RMSE = NULL, RME = NULL){
    result <- list(
      RMSE = RMSE,
      RME = RME
    )
    class(result) <- append(class(result), "multiResultClass")
    return(result)
  }
  
  cl <- makeCluster(5)
  registerDoParallel(cl) 
  Model.Comparison <- 
    foreach (i = data.loop, .packages = c("grpreg", "glmnet", "fastDummies", "RcppArmadillo", "MASS", "Matrix", "TmpLasso")) %dopar% {
      sd.gamma <- 1
      m <- 50
      poisson.mean <- 100
      sd.z <- 1
      true.beta <- c(rep(beta, ls$n.beta * ls$prop.NonZero.beta), 
                     rep(0, ls$n.beta * (1 - ls$prop.NonZero.beta)))  # 10/50 beta's are non-zero
      rho <- ls$rho
      
      data.proveff3 <- function(m, prov.size, gamma, sd.gamma, sd.z,
                                beta, Y.char, Z.char, prov.char, rho) {
        N <- sum(prov.size) # total number of discharges
        n.beta <- length(beta)
        gamma.rep <- rep(gamma, prov.size)
        prov <- rep(1:m, prov.size) # provider IDs
        KW2013 <- function(i, rho, n.beta){
          MASS::mvrnorm(n = prov.size[i], mu = (sd.z * gamma[i] * rho / sd.gamma) * matrix(1, nrow = n.beta),
                        Sigma = sd.z^2 * (diag(1 - rho, n.beta) + (rho - rho^2) * matrix(1, ncol = n.beta, nrow = n.beta)))
        }
        Z <- do.call(rbind, lapply(1:m, function(i) KW2013(i, rho, n.beta)))
        Y <- rbinom(N, 1, plogis(gamma.rep + Z %*% beta))
        mu <- plogis(gamma.rep + Z %*% beta)
        data <- data.frame(Y, prov, Z, mu, stringsAsFactors = F)
        colnames(data) <- c(Y.char, prov.char, Z.char, "true.mu")
        return(data)
      }
      prov.size <- pmax(c(rpois(m, poisson.mean)), 20)
      gamma <- rnorm(m, 0, sd.gamma)
      Y.char <- 'Y'
      prov.char <- 'unit'
      Z.char <- paste0('z', 1:length(true.beta))
      sim_data <- data.proveff3(m, prov.size, gamma, sd.gamma, sd.z, true.beta, Y.char, Z.char, prov.char, rho)
      true.mu <- sim_data$true.mu
      
      # results from grpLasso
      cv.model_pp_lasso <- cv.pp.lasso(sim_data, Y.char, Z.char, prov.char, trace.lambda = F,
                                       nfolds = 10, trace.cv = F, MM = F)
      cv_BestModel_pp_lasso <- cv.model_pp_lasso$fit
      best.beta.pp_lasso <- cv_BestModel_pp_lasso$beta[, cv.model_pp_lasso$min]
      best.eta.pp_lasso <- cv_BestModel_pp_lasso$linear.predictors[, cv.model_pp_lasso$min]

      RMSE.pp_lasso <- sqrt(sum((best.beta.pp_lasso - true.beta) ^ 2) / ls$n.beta) #RMSE
      RME.pp_lasso <- sqrt(sum((plogis(best.eta.pp_lasso) - true.mu) ^ 2) / nrow(sim_data)) #RME
      
      dummy_data <- dummy_cols(sim_data, select_columns = prov.char, remove_selected_columns = TRUE, 
                               remove_first_dummy = TRUE)
      ID.char <- rep(NA, m - 1)
      for (i in 1:(m - 1)){
        ID.char[i] <- paste0("unit_", i + 1)
      }
      
      ## results from grpreg
      cv.model_grpreg <- cv.grpreg(dummy_data[,c(Z.char, ID.char)], dummy_data[,Y.char], family = "binomial",
                                   penalty = "grLasso", group = c(1:length(Z.char), rep(0, length(ID.char))), alpha = 1, 
                                   nfolds = 10, trace.cv = F)
      cv_BestModel_grpreg <- cv.model_grpreg$fit
      best.beta.grpreg <- cv_BestModel_grpreg$beta[2:(1 + ls$n.beta), cv.model_grpreg$min]
      best.eta.grpreg <- cv_BestModel_grpreg$linear.predictors[, cv.model_grpreg$min]
      RMSE.grpreg <- sqrt(sum((best.beta.grpreg - true.beta) ^ 2) / ls$n.beta)  #RMSE
      RME.grpreg <- sqrt(sum((plogis(best.eta.grpreg) - true.mu) ^ 2) / nrow(sim_data))  #RME
      
      
      
      ## results from glmnet     
      cv.model_glmnet <- cv.glmnet(as.matrix(dummy_data[ ,c(Z.char, ID.char)]), dummy_data[,Y.char], family = "binomial", 
                                   penalty.factor = c(rep(1, length(Z.char)), rep(0, length(ID.char))), 
                                   alpha = 1, nfolds = 10)
      
      best.beta.glmnet <- coef(cv.model_glmnet, s = "lambda.min")[2:(1 + ls$n.beta), ]
      best.eta.glmnet <- predict(cv.model_glmnet, newx = as.matrix(dummy_data[ ,c(Z.char, ID.char)]), s = "lambda.min")
      RMSE.glmnet <- sqrt(sum((best.beta.glmnet - true.beta) ^ 2) / ls$n.beta)  #RMSE
      RME.glmnet <- sqrt(sum((plogis(best.eta.glmnet) - true.mu) ^ 2) / nrow(sim_data))  #RME
      
      
      result <- multiResultClass()
      result$RMSE <- matrix(c(RMSE.pp_lasso, RMSE.grpreg, RMSE.glmnet), nrow = 3)
      result$RME <- matrix(c(RME.pp_lasso, RME.grpreg, RME.glmnet), nrow = 3)
      return(result)
    }
  stopCluster(cl)
  
  n.data.loop <- length(data.loop)
  
  RME.Mean <- as.data.frame(matrix(rep(0, length(beta.sequence) * 3), nrow = 3))
  colnames(RME.Mean) <- paste0("beta = ", beta.sequence)
  rownames(RME.Mean) <- c("pplasso", "grpreg", "glmnet")
  
  RME.sd <- as.data.frame(matrix(rep(0, length(beta.sequence) * 3), nrow = 3))
  colnames(RME.sd) <- paste0("beta = ", beta.sequence)
  rownames(RME.sd) <- c("pplasso", "grpreg", "glmnet")
  
  RMSE.Mean <- as.data.frame(matrix(rep(0, length(beta.sequence) * 3), nrow = 3))
  colnames(RMSE.Mean) <- paste0("beta = ", beta.sequence)
  rownames(RMSE.Mean) <- c("pplasso", "grpreg", "glmnet")
  
  RMSE.sd <- as.data.frame(matrix(rep(0, length(beta.sequence) * 3), nrow = 3))
  colnames(RMSE.sd) <- paste0("beta = ", beta.sequence)
  rownames(RMSE.sd) <- c("pplasso", "grpreg", "glmnet")
  
  RMSE <- matrix(rep(0, 3 * n.data.loop), nrow = 3)
  rownames(RMSE) <- c("pplasso", "grpreg", "glmnet")
  colnames(RMSE) <- paste0("Data_", data.loop)
  for (i in data.loop){
    RMSE[, i] <- Model.Comparison[[i]]$RMSE
  }
  RMSE.Mean[, ind] <- round(apply(RMSE, 1, mean), digits = 3)
  RMSE.sd[, ind] <- round(apply(RMSE, 1, sd), digits = 3)
  RMSE.df <- RMSE.Mean
  RMSE.df[, ind] <- paste0(RMSE.Mean[, ind], " (", RMSE.sd[, ind], ") ")
  
  
  RME <- matrix(rep(0, 3 * n.data.loop), nrow = 3)
  rownames(RME) <- c("pplasso", "grpreg", "glmnet")
  colnames(RME) <- paste0("Data_", data.loop)
  for (i in data.loop){
    RME[, i] <- Model.Comparison[[i]]$RME
  }
  RME.Mean[, ind] <- round(apply(RME, 1, mean), digits = 3)
  RME.sd[, ind] <- round(apply(RME, 1, sd), digits = 3)
  RME.df <- RME.Mean
  RME.df[, ind] <- paste0(RME.Mean[, ind], " (", RME.sd[, ind], ") ")
  
  return(list(RMSE.df = RMSE.df,
              RME.df = RME.df))
}


ls <- list(rho = 0.6, n.beta = 50, prop.NonZero.beta = 0.2)
beta.sequence <- c(0.5, 1, 2)

### simulation

pplasso.beta05 <- lasso.RMSE(ls, beta.sequence, ind = 1, n.rep = 100)
save(pplasso.beta05, file = paste0("pplasso_RME_RMSE_beta05_", Sys.Date(), ".RData"))

pplasso.beta1 <- lasso.RMSE(ls, beta.sequence, ind = 2, n.rep = 100)
save(pplasso.beta1, file = paste0("pplasso_RME_RMSE_beta1_", Sys.Date(), ".RData"))

pplasso.beta2 <- lasso.RMSE(ls, beta.sequence, ind = 3, n.rep = 100)
save(pplasso.beta2, file = paste0("pplasso_RME_RMSE_beta2_", Sys.Date(), ".RData"))


pplasso.RMSE <- cbind(pplasso.beta05$RMSE.df[, 1, drop = F], pplasso.beta1$RMSE.df[, 2, drop = F], pplasso.beta2$RMSE.df[, 3, drop = F])
pplasso.RME <- cbind(pplasso.beta05$RME.df[, 1, drop = F], pplasso.beta1$RME.df[, 2, drop = F], pplasso.beta2$RME.df[, 3, drop = F])

save(pplasso.RMSE, pplasso.RME, file = paste0("pplasso_RME_RMSE_combine_", Sys.Date(), ".RData"))






