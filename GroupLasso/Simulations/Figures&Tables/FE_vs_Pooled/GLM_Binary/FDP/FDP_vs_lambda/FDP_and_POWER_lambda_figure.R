library(MASS)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(grid)
library(dplyr)

setwd("/home/ybshao/GroupLasso/Simulations/Figures&Tables/FE_vs_Pooled/GLM_Binary/FDP/FDP_vs_lambda")

sim4 <- function(ls){
  print(paste0("Iteration for lambda = ", ls$lambda, " starts.."))
  start.time <- Sys.time()
  sd.gamma <- 1
  m <- ls$m
  sd.z <- 1
  true.beta <- c(rep(ls$beta, ls$n.beta * ls$prop.NonZero.beta), 
                 rep(0, ls$n.beta * (1 - ls$prop.NonZero.beta)))  # 5% beta's are non-zero
  loop_length <- ls$n.rep
  rho <- ls$rho
  sim.continuous <- function() {
    data.proveff4 <- function(m, prov.size, gamma, sd.gamma, sd.z,
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
    lambda <- ls$lambda
    prov.size <- pmax(c(rpois(m, ls$poisson.mean)), 20)
    gamma <- rnorm(m, 0, sd.gamma)
    Y.char <- 'Y'
    prov.char <- 'unit'
    Z.char <- paste0('z', 1:length(true.beta))
    tb <- data.proveff4(m, prov.size, gamma, sd.gamma, sd.z, true.beta, Y.char, Z.char, prov.char, rho)
    
    true.non.zero.beta.index <- 1:(ls$n.beta * ls$prop.NonZero.beta)
    true.zero.beta.index <- (ls$n.beta * ls$prop.NonZero.beta + 1):ls$n.beta
    
    # (1) Fixed effect model 
    FE.lasso.lambda <- pp.lasso(tb, Y.char, Z.char, prov.char, lambda = lambda, MM = T)
    FE.beta <- FE.lasso.lambda$beta
    FE.estimated.non.zero.beta.index <- which(FE.beta != 0)
    FE.estimated.zero.beta.index <- which(FE.beta == 0)

    if (length(FE.estimated.non.zero.beta.index) == 0) {
      FE.FDP <- 0
    } else {
      FE.FDP <- sum(!(FE.estimated.non.zero.beta.index %in% true.non.zero.beta.index))/length(FE.estimated.non.zero.beta.index)
    }
    FE.FNP <- sum(!(FE.estimated.zero.beta.index %in% true.zero.beta.index))/length(FE.estimated.zero.beta.index)

    
    # (2) Pooled model 
    pooled.lasso.lambda <- pp.lasso(tb, Y.char, Z.char, lambda = lambda, MM = T)
    pooled.beta <- pooled.lasso.lambda$beta
    pooled.estimated.non.zero.beta.index <- which(pooled.beta != 0)
    pooled.estimated.zero.beta.index <- which(pooled.beta == 0)
    
    if (length(pooled.estimated.non.zero.beta.index) == 0){
      pooled.FDP <- 0
    } else {
      pooled.FDP <- sum(!(pooled.estimated.non.zero.beta.index %in% true.non.zero.beta.index))/length(pooled.estimated.non.zero.beta.index)
    }
    pooled.FNP <- sum(!(pooled.estimated.zero.beta.index %in% true.zero.beta.index))/length(pooled.estimated.zero.beta.index)

    res <- c(FE.FDP, FE.FNP, pooled.FDP, pooled.FNP) 
    return(res)
  }
  
  `%dopar%` <- foreach::`%dopar%`
  cl <- parallel::makeCluster(5)
  doParallel::registerDoParallel(cl)
  system.time(
    mat <- foreach::foreach(i = 1:loop_length, .combine = cbind, .packages = c("TmpLasso", "Matrix")) %dopar% {
      sim.continuous()})
  parallel::stopCluster(cl)
  
  Mean.Performance <- rowMeans(mat)
  Sd.Performance <- apply(mat, 1, sd)
  
  rst <- c(ls$m, ls$poisson.mean, rho, ls$lambda, Mean.Performance, Sd.Performance) 
  
  names(rst) <- c("num.units", "possion.mean", "rho", "lambda",
                  "FE.FDP.Mean", "FE.FNP.Mean", "pooled.FDP.Mean", "pooled.FNP.Mean",
                  "FE.FDP.Sd", "FE.FNP.Sd", "pooled.FDP.Sd", "pooled.FNP.Sd")
  end.time <- Sys.time()
  process.time <- difftime(end.time, start.time, units = 'mins')
  print(paste0("Processing time: ", round(process.time, digits = 4), " mins"))
  return(rst)
}

shared.lambda.seq <- exp(seq(log(1e-01), log(1e-03), length = 100))  #fixed lambda sequence, from 0.1 to 0.001

### (1) beta = 0.5
sim4.res.beta1 <- NULL
for (lambda in shared.lambda.seq){
  ls <- list(m = 200, poisson.mean = 100, rho = 0.6, n.rep = 200, n.beta = 100, prop.NonZero.beta = 0.1, beta = 0.5, lambda = lambda)
  sim4.res.beta1 <- rbind(sim4.res.beta1, sim4(ls))
}
sim4.res.beta1 <- as.data.frame(sim4.res.beta1)

## FDP
FE.FDP.df.beta1 <- cbind(log(sim4.res.beta1["lambda"]), sim4.res.beta1["FE.FDP.Mean"], sim4.res.beta1["FE.FDP.Mean"] - sim4.res.beta1["FE.FDP.Sd"], 
                         sim4.res.beta1["FE.FDP.Mean"] + sim4.res.beta1["FE.FDP.Sd"])
colnames(FE.FDP.df.beta1) <- NA
pooled.FDP.df.beta1 <- cbind(log(sim4.res.beta1["lambda"]), sim4.res.beta1["pooled.FDP.Mean"], sim4.res.beta1["pooled.FDP.Mean"] - sim4.res.beta1["pooled.FDP.Sd"], 
                             sim4.res.beta1["pooled.FDP.Mean"] + sim4.res.beta1["pooled.FDP.Sd"])
colnames(pooled.FDP.df.beta1) <- NA
FDP.plot.df.beta1 <- rbind(FE.FDP.df.beta1, pooled.FDP.df.beta1)
colnames(FDP.plot.df.beta1) <- c("lambda", "FD.Mean", "FD.lower", "FD.upper")
FDP.plot.df.beta1$model <- rep(c("FE", "Pooled"), each = length(shared.lambda.seq))


FDP.plot.beta1 <- ggplot(FDP.plot.df.beta1, aes(lambda, group = factor(model))) +
  geom_ribbon(aes(ymin = FD.lower, ymax = FD.upper, fill = factor(model)), alpha = 0.5) +
  geom_line(aes(y = FD.Mean, color = factor(model), linetype = factor(model)), size = 1)  + 
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black")) + 
  theme(axis.title = element_text(size = 24, family = "serif", face = "bold"),
        legend.text = element_text(size = 14, family = "serif")) + 
  theme(axis.text = element_text(face = "bold.italic", size = 17, family = "serif")) + 
  theme(legend.position = c(0.78, 0.22),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank()) + 
  theme(legend.key.height = unit(0.6, 'cm'),
        legend.key.width = unit(1.2, 'cm'),
        legend.spacing.y = unit(0, 'cm')) +
  labs(title = "", 
       x = bquote(log(lambda)),
       y = bquote("FDP")) +
  scale_linetype_manual(values = c("solid", "dotdash"), name = "", labels = c("Fixed Effect", "Pooled")) + 
  scale_color_manual(values = c("blue", "red"), name = "", labels = c("Fixed Effect", "Pooled")) + 
  scale_fill_manual(values = c("grey85", "grey90"), name = "", labels = c("Fixed Effect", "Pooled")) +
  scale_x_continuous(trans = scales::reverse_trans(), breaks = round(seq(round(max(log(shared.lambda.seq)), 0), round(min(log(shared.lambda.seq)), 0), by = - 1), 1)) +
  ylim(-0.05, 0.9)

## FNP
FE.FNP.df.beta1 <- cbind(log(sim4.res.beta1["lambda"]), sim4.res.beta1["FE.FNP.Mean"], sim4.res.beta1["FE.FNP.Mean"] - sim4.res.beta1["FE.FNP.Sd"], 
                         sim4.res.beta1["FE.FNP.Mean"] + sim4.res.beta1["FE.FNP.Sd"])
colnames(FE.FNP.df.beta1) <- NA
pooled.FNP.df.beta1 <- cbind(log(sim4.res.beta1["lambda"]), sim4.res.beta1["pooled.FNP.Mean"], sim4.res.beta1["pooled.FNP.Mean"] - sim4.res.beta1["pooled.FNP.Sd"], 
                             sim4.res.beta1["pooled.FNP.Mean"] + sim4.res.beta1["pooled.FNP.Sd"])
colnames(pooled.FNP.df.beta1) <- NA
FNP.plot.df.beta1 <- rbind(FE.FNP.df.beta1, pooled.FNP.df.beta1)
colnames(FNP.plot.df.beta1) <- c("lambda", "FN.Mean", "FN.lower", "FN.upper")
FNP.plot.df.beta1$model <- rep(c("FE", "Pooled"), each = length(shared.lambda.seq))


FNP.plot.beta1 <- ggplot(FNP.plot.df.beta1, aes(lambda, group = factor(model))) +
  geom_ribbon(aes(ymin = FN.lower, ymax = FN.upper, fill = factor(model)), alpha = 0.5) +
  geom_line(aes(y = FN.Mean, color = factor(model), linetype = factor(model)), size = 1)  + 
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black")) + 
  theme(axis.title = element_text(size = 24, family = "serif", face = "bold"),
        legend.text = element_text(size = 14, family = "serif")) + 
  theme(axis.text = element_text(face = "bold.italic", size = 17, family = "serif")) + 
  theme(legend.position = c(0.78, 0.22),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank()) + 
  theme(legend.key.height = unit(0.6, 'cm'),
        legend.key.width = unit(1.2, 'cm'),
        legend.spacing.y = unit(0, 'cm')) +
  labs(title = "", 
       x = bquote(log(lambda)),
       y = bquote("FNP")) +
  scale_linetype_manual(values = c("solid", "dotdash"), name = "", labels = c("Fixed Effect", "Pooled")) + 
  scale_color_manual(values = c("blue", "red"), name = "", labels = c("Fixed Effect", "Pooled")) + 
  scale_fill_manual(values = c("grey85", "grey90"), name = "", labels = c("Fixed Effect", "Pooled")) +
  scale_x_continuous(trans = scales::reverse_trans(), breaks = round(seq(round(max(log(shared.lambda.seq)), 0), round(min(log(shared.lambda.seq)), 0), by = - 1), 1)) +
  ylim(-0.005, 0.1)



### (2) beta = 1
sim4.res.beta2 <- NULL
for (lambda in shared.lambda.seq){
  ls <- list(m = 200, poisson.mean = 100, rho = 0.6, n.rep = 500, n.beta = 100, prop.NonZero.beta = 0.1, beta = 1, lambda = lambda)
  sim4.res.beta2 <- rbind(sim4.res.beta2, sim4(ls))
}
sim4.res.beta2 <- as.data.frame(sim4.res.beta2)


##FDP
FE.FDP.df.beta2 <- cbind(log(sim4.res.beta2["lambda"]), sim4.res.beta2["FE.FDP.Mean"], sim4.res.beta2["FE.FDP.Mean"] - sim4.res.beta2["FE.FDP.Sd"], 
                         sim4.res.beta2["FE.FDP.Mean"] + sim4.res.beta2["FE.FDP.Sd"])
colnames(FE.FDP.df.beta2) <- NA
pooled.FDP.df.beta2 <- cbind(log(sim4.res.beta2["lambda"]), sim4.res.beta2["pooled.FDP.Mean"], sim4.res.beta2["pooled.FDP.Mean"] - sim4.res.beta2["pooled.FDP.Sd"], 
                             sim4.res.beta2["pooled.FDP.Mean"] + sim4.res.beta2["pooled.FDP.Sd"])
colnames(pooled.FDP.df.beta2) <- NA
FDP.plot.df.beta2 <- rbind(FE.FDP.df.beta2, pooled.FDP.df.beta2)
colnames(FDP.plot.df.beta2) <- c("lambda", "FD.Mean", "FD.lower", "FD.upper")
FDP.plot.df.beta2$model <- rep(c("FE", "Pooled"), each = length(shared.lambda.seq))


FDP.plot.beta2 <- ggplot(FDP.plot.df.beta2, aes(lambda, group = factor(model))) +
  geom_ribbon(aes(ymin = FD.lower, ymax = FD.upper, fill = factor(model)), alpha = 0.5) +
  geom_line(aes(y = FD.Mean, color = factor(model), linetype = factor(model)), size = 1)  + 
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black")) + 
  theme(axis.title = element_text(size = 24, family = "serif", face = "bold"),
        legend.text = element_text(size = 14, family = "serif")) + 
  theme(axis.text = element_text(face = "bold.italic", size = 17, family = "serif")) + 
  theme(legend.position = c(0.78, 0.22),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank()) + 
  theme(legend.key.height = unit(0.6, 'cm'),
        legend.key.width = unit(1.2, 'cm'),
        legend.spacing.y = unit(0, 'cm')) +
  labs(title = "", 
       x = bquote(log(lambda)),
       y = bquote("FDP")) +
  scale_linetype_manual(values = c("solid", "dotdash"), name = "", labels = c("Fixed Effect", "Pooled")) + 
  scale_color_manual(values = c("blue", "red"), name = "", labels = c("Fixed Effect", "Pooled")) + 
  scale_fill_manual(values = c("grey85", "grey90"), name = "", labels = c("Fixed Effect", "Pooled")) +
  scale_x_continuous(trans = scales::reverse_trans(), breaks = round(seq(round(max(log(shared.lambda.seq)), 0), round(min(log(shared.lambda.seq)), 0), by = - 1), 1)) +
  ylim(-0.05, 0.9)

## FNP
FE.FNP.df.beta2 <- cbind(log(sim4.res.beta2["lambda"]), sim4.res.beta2["FE.FNP.Mean"], sim4.res.beta2["FE.FNP.Mean"] - sim4.res.beta2["FE.FNP.Sd"], 
                         sim4.res.beta2["FE.FNP.Mean"] + sim4.res.beta2["FE.FNP.Sd"])
colnames(FE.FNP.df.beta2) <- NA
pooled.FNP.df.beta2 <- cbind(log(sim4.res.beta2["lambda"]), sim4.res.beta2["pooled.FNP.Mean"], sim4.res.beta2["pooled.FNP.Mean"] - sim4.res.beta2["pooled.FNP.Sd"], 
                             sim4.res.beta2["pooled.FNP.Mean"] + sim4.res.beta2["pooled.FNP.Sd"])
colnames(pooled.FNP.df.beta2) <- NA
FNP.plot.df.beta2 <- rbind(FE.FNP.df.beta2, pooled.FNP.df.beta2)
colnames(FNP.plot.df.beta2) <- c("lambda", "FN.Mean", "FN.lower", "FN.upper")
FNP.plot.df.beta2$model <- rep(c("FE", "Pooled"), each = length(shared.lambda.seq))


FNP.plot.beta2 <- ggplot(FNP.plot.df.beta2, aes(lambda, group = factor(model))) +
  geom_ribbon(aes(ymin = FN.lower, ymax = FN.upper, fill = factor(model)), alpha = 0.5) +
  geom_line(aes(y = FN.Mean, color = factor(model), linetype = factor(model)), size = 1)  + 
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black")) + 
  theme(axis.title = element_text(size = 24, family = "serif", face = "bold"),
        legend.text = element_text(size = 14, family = "serif")) + 
  theme(axis.text = element_text(face = "bold.italic", size = 17, family = "serif")) + 
  theme(legend.position = c(0.78, 0.22),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank()) + 
  theme(legend.key.height = unit(0.6, 'cm'),
        legend.key.width = unit(1.2, 'cm'),
        legend.spacing.y = unit(0, 'cm')) +
  labs(title = "", 
       x = bquote(log(lambda)),
       y = bquote("FNP")) +
  scale_linetype_manual(values = c("solid", "dotdash"), name = "", labels = c("Fixed Effect", "Pooled")) + 
  scale_color_manual(values = c("blue", "red"), name = "", labels = c("Fixed Effect", "Pooled")) + 
  scale_fill_manual(values = c("grey85", "grey90"), name = "", labels = c("Fixed Effect", "Pooled")) +
  scale_x_continuous(trans = scales::reverse_trans(), breaks = round(seq(round(max(log(shared.lambda.seq)), 0), round(min(log(shared.lambda.seq)), 0), by = - 1), 1)) + 
  ylim(-0.005, 0.1)



plt.beta1 <- annotate_figure(ggarrange(FDP.plot.beta1,
                                       FNP.plot.beta1,
                                       nrow = 2, ncol = 1, common.legend = F, 
                                       labels = c("A.1", "A.2")))

plt.beta2 <- annotate_figure(ggarrange(FDP.plot.beta2,
                                       FNP.plot.beta2,
                                       nrow = 2, ncol = 1, common.legend = F, 
                                       labels = c("B.1", "B.2")))

save(sim4.res.beta1, sim4.res.beta2,
     FDP.plot.df.beta1, FDP.plot.beta1, 
     FNP.plot.df.beta1, FNP.plot.beta1,
     FDP.plot.df.beta2, FDP.plot.beta2, 
     FNP.plot.df.beta2, FNP.plot.beta2, 
     plt.beta1, plt.beta2,
     file = paste0("FE_vs_Pooled_lambda_vs_FDP_figure_", Sys.Date(), ".RData"))
