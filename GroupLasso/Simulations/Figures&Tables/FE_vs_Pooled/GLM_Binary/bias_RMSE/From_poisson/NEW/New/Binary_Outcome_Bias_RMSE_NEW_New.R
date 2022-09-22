library(MASS)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(grid)
library(dplyr)
setwd("/home/ybshao/GroupLasso/Simulations/Figures&Tables/FE_vs_Pooled/GLM_Binary/bias_RMSE/From_poisson/NEW/New")

######----------- Simulation Functions ------------######
# ls <- list(m = 200, poisson.mean = 100, rho = rho, n.rep = 10000, sd.gamma = 1, sd.z = 0.3)
sim <- function(ls){
  print(paste0("Iteration for rho = ", ls$rho, " starts.."))
  start.time <- Sys.time()
  m <- ls$m
  sd.z <- ls$sd.z
  sd.gamma <- ls$sd.gamma
  beta <- 1
  beta_true <- beta
  loop_length <- ls$n.rep
  rho <- ls$rho
  
  sim.continuous <- function() {
    data.proveff2 <- function(m, prov.size, gamma, sd.gamma, sd.z,
                              beta, Y.char, Z.char, prov.char, rho) {
      N <- sum(prov.size) # total number of discharges
      gamma.rep <- rep(gamma, times = prov.size)
      prov <- rep(paste0("prov", sprintf(paste0("%0", nchar(m), "d"), 1:m)),
                  times = prov.size) # provider IDs
      KW2013 <- function(i, rho){
        rnorm(n = prov.size[i], sd.z * gamma[i] * rho / sd.gamma, sd.z * sqrt(1 - rho ^ 2))
      }
      Z <- do.call(c, lapply(1:m, function(i) KW2013(i, rho)))
      Y <- rbinom(N, 1, plogis(gamma.rep + Z * beta))
      data <- data.frame(Y, prov, Z, stringsAsFactors = F)
      colnames(data) <- c(Y.char, prov.char, Z.char)
      return(data)
    }
    prov.size <- pmax(c(rpois(m, ls$poisson.mean)), 20)
    gamma <- rnorm(m, 0, sd.gamma)
    Y.char <- 'Y'
    prov.char <- 'unit'
    Z.char <- "X"
    tb <- data.proveff2(m, prov.size, gamma, sd.gamma, sd.z, beta,
                        Y.char, Z.char, prov.char, rho)
    
    #fit two model and obtain estimators
    pool <- pp.lasso(tb, Y.char, Z.char, lambda = 0, MM = T)
    beta_pool <- pool$beta[1]
    
    fixed <- pp.lasso(tb, Y.char, Z.char, prov.char, lambda = 0, MM = T)
    beta_fe <- fixed$beta[1]
    
    c(beta_pool, beta_fe) 
  }
  
  `%dopar%` <- foreach::`%dopar%`
  cl <- parallel::makeCluster(4)
  doParallel::registerDoParallel(cl)
  system.time(
    mat <- foreach::foreach(i = 1:loop_length, .combine = cbind, .packages = c("TmpLasso", "Matrix")) %dopar% {
      sim.continuous()})
  parallel::stopCluster(cl)

  #get statistics after loops end
  bias.beta <- rowMeans(mat - beta_true)
  ESD.beta <- apply(mat, 1, sd)
  RMSE.beta <- sqrt(bias.beta ^ 2 + ESD.beta ^ 2)
  sd.z.gamma <- paste0("sd.gamma = ", sd.gamma, "; sd.z = ", sd.z)
  rst <- c(m, sd.z.gamma, rho, bias.beta, RMSE.beta) 
  names(rst) <- c("num.units", "sd.z.gamma", "rho", 
                  "bias.beta.pooled", "bias.beta.FE",
                  "RMSE.beta.pooled", "RMSE.beta.FE")
  end.time <- Sys.time()
  process.time <- difftime(end.time, start.time, units = 'mins')
  print(paste0("Processing time: ", round(process.time, digits = 4), " mins (current rho = ", rho, ")"))
  return(rst)
}

######----------- Simulations: m = 50 ------------######
# 1. 50 units & poisson(100) & sd.gamma = 2 & sd.z = 0.3;
sim.res <- NULL
for (rho in seq(-8, 8, 1)/10){
  ls <- list(m = 50, poisson.mean = 100, rho = rho, n.rep = 10000, sd.gamma = 2, sd.z = 0.3)
  sim.res <- rbind(sim.res, sim(ls))
}
sim.res.m50.2.03 <- as.data.frame(sim.res)

# 2. 50 units & poisson(100) & sd.gamma = 1 & sd.z = 0.3;
sim.res <- NULL
for (rho in seq(-8, 8, 1)/10){
  ls <- list(m = 50, poisson.mean = 100, rho = rho, n.rep = 10000, sd.gamma = 1, sd.z = 0.3)
  sim.res <- rbind(sim.res, sim(ls))
}
sim.res.m50.1.03 <- as.data.frame(sim.res)

# 3. 50 units & poisson(100) & sd.gamma = 1 & sd.z = 1;
sim.res <- NULL
for (rho in seq(-8, 8, 1)/10){
  ls <- list(m = 50, poisson.mean = 100, rho = rho, n.rep = 10000, sd.gamma = 2, sd.z = 1)
  sim.res <- rbind(sim.res, sim(ls))
}
sim.res.m50.2.1 <- as.data.frame(sim.res)

Binary.combine.m50 <- rbind(sim.res.m50.2.03, sim.res.m50.1.03, sim.res.m50.2.1, stringsAsFactors = F)

save(Binary.combine.m50, sim.res.m50.2.03, sim.res.m50.1.03, sim.res.m50.2.1, 
     file = paste0("Binary_m50_Bias_RMSE_", Sys.Date(), ".RData"))


######----------- Simulations: m = 100 ------------######
# 1. 100 units & poisson(100) & sd.gamma = 2 & sd.z = 0.3;
sim.res <- NULL
for (rho in seq(-8, 8, 1)/10){
  ls <- list(m = 100, poisson.mean = 100, rho = rho, n.rep = 10000, sd.gamma = 2, sd.z = 0.3)
  sim.res <- rbind(sim.res, sim(ls))
}
sim.res.m100.2.03 <- as.data.frame(sim.res)

# 2. 100 units & poisson(100) & sd.gamma = 1 & sd.z = 0.3;
sim.res <- NULL
for (rho in seq(-8, 8, 1)/10){
  ls <- list(m = 100, poisson.mean = 100, rho = rho, n.rep = 10000, sd.gamma = 1, sd.z = 0.3)
  sim.res <- rbind(sim.res, sim(ls))
}
sim.res.m100.1.03 <- as.data.frame(sim.res)

# 3. 100 units & poisson(100) & sd.gamma = 2 & sd.z = 1;
sim.res <- NULL
for (rho in seq(-8, 8, 1)/10){
  ls <- list(m = 100, poisson.mean = 100, rho = rho, n.rep = 10000, sd.gamma = 2, sd.z = 1)
  sim.res <- rbind(sim.res, sim(ls))
}
sim.res.m100.2.1 <- as.data.frame(sim.res)

Binary.combine.m100 <- rbind(sim.res.m100.2.03, sim.res.m100.1.03, sim.res.m100.2.1, stringsAsFactors = F)

save(Binary.combine.m100, sim.res.m100.2.03, sim.res.m100.1.03, sim.res.m100.2.1, 
     file = paste0("Binary_m100_Bias_RMSE_", Sys.Date(), ".RData"))

######----------- Simulations: m = 200 ------------######
# 1. 200 units & poisson(100) & sd.gamma = 2 & sd.z = 0.3;
sim.res <- NULL
for (rho in seq(-8, 8, 1)/10){
  ls <- list(m = 200, poisson.mean = 100, rho = rho, n.rep = 10000, sd.gamma = 2, sd.z = 0.3)
  sim.res <- rbind(sim.res, sim(ls))
}
sim.res.m200.2.03 <- as.data.frame(sim.res)

# 2. 200 units & poisson(100) & sd.gamma = 1 & sd.z = 0.3;
sim.res <- NULL
for (rho in seq(-8, 8, 1)/10){
  ls <- list(m = 200, poisson.mean = 100, rho = rho, n.rep = 10000, sd.gamma = 1, sd.z = 0.3)
  sim.res <- rbind(sim.res, sim(ls))
}
sim.res.m200.1.03 <- as.data.frame(sim.res)

# 3. 200 units & poisson(100) & sd.gamma = 2 & sd.z = 1;
sim.res <- NULL
for (rho in seq(-8, 8, 1)/10){
  ls <- list(m = 200, poisson.mean = 100, rho = rho, n.rep = 10000, sd.gamma = 2, sd.z = 1)
  sim.res <- rbind(sim.res, sim(ls))
}
sim.res.m200.2.1 <- as.data.frame(sim.res)

Binary.combine.m200 <- rbind(sim.res.m200.2.03, sim.res.m200.1.03, sim.res.m200.2.1, stringsAsFactors = F)

save(Binary.combine.m200, sim.res.m200.2.03, sim.res.m200.1.03, sim.res.m200.2.1, 
     file = paste0("Binary_m200_Bias_RMSE_", Sys.Date(), ".RData"))



######----------- Corresponding Plot Functions ------------######
plot.function.bias <- function(res.df, label){
  bias.df <- res.df[, c("sd.z.gamma", "rho", "bias.beta.pooled")]
  bias.df$rho <- as.numeric(as.character(bias.df$rho))
  bias.df$bias.beta.pooled <- as.numeric(as.character(bias.df$bias.beta.pooled))

  Bias.plot <- ggplot(bias.df, aes(rho, group = factor(sd.z.gamma))) +
    geom_abline(color = "black", linetype = "solid", size = 1.5, alpha = 0.8, intercept = 0, slope = 0) +
    geom_line(aes(y = bias.beta.pooled, color = factor(sd.z.gamma), 
                  linetype = factor(sd.z.gamma)), size = 1)  + 
    theme(panel.grid = element_blank(), panel.background = element_blank(),
          axis.line = element_line(colour = "black")) + 
    theme(axis.title = element_text(size = 24, family = "serif", face = "bold"),
          legend.text = element_text(size = 14, family = "serif")) + 
    theme(axis.text = element_text(face = "bold.italic", size = 17, family = "serif")) + 
    theme(legend.position = c(0.28, 0.78),
          legend.background = element_rect(fill = "white", color = "black"),
          legend.title = element_blank()) + 
    theme(legend.key.height = unit(0.6, 'cm'),
          legend.key.width = unit(1.2, 'cm'),
          legend.spacing.y = unit(0, 'cm')) +
    labs(title = "", 
         x = bquote(rho),
         y = "Bias") +
    scale_x_continuous(breaks = seq(-8, 8, 2)/10) + 
    scale_linetype_manual(values = c("dotted", "dashed", "dotdash"), name = "", 
                          labels = c(bquote(sigma[gamma] ~ " = 2, " ~ sigma[X] ~ " = 0.3"), 
                                     bquote(sigma[gamma] ~ " = 1, " ~ sigma[X] ~ " = 0.3"), 
                                     bquote(sigma[gamma] ~ " = 2, " ~ sigma[X] ~ " = 1"))) + 
    scale_color_manual(values = c("red", "blue", "darkorange"), name = "", 
                       labels = c(bquote(sigma[gamma] ~ " = 2, " ~ sigma[X] ~ " = 0.3"), 
                                  bquote(sigma[gamma] ~ " = 1, " ~ sigma[X] ~ " = 0.3"), 
                                  bquote(sigma[gamma] ~ " = 2, " ~ sigma[X] ~ " = 1"))) +
    ylim(-4.5, 4.5)
  Bias.plot <- annotate_figure(ggarrange(Bias.plot,
                                         nrow = 1, ncol = 1,
                                         labels = label))
  return(Bias.plot)
}

plot.bias.binary.m50 <- plot.function.bias(Binary.combine.m50, "A")
plot.bias.binary.m100 <- plot.function.bias(Binary.combine.m100, "B")
plot.bias.binary.m200 <- plot.function.bias(Binary.combine.m200, "C")

save(plot.bias.binary.m50,
     plot.bias.binary.m100,
     plot.bias.binary.m200,
     file = paste0("Binary_bias_", Sys.Date(), ".RData"))
