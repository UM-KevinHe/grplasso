library(MASS)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(grid)
library(dplyr)
setwd("/home/ybshao/GroupLasso/Simulations/Figures&Tables/FE_vs_Pooled/GLM_Binary/bias_RMSE/From_poisson/NEW")

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
  rst <- c(ls$m, ls$poisson.mean, rho, bias.beta, RMSE.beta) 
  names(rst) <- c("num.units", "poisson.mean", "rho", 
                  "bias.beta.pooled", "bias.beta.FE",
                  "RMSE.beta.pooled", "RMSE.beta.FE")
  end.time <- Sys.time()
  process.time <- difftime(end.time, start.time, units = 'mins')
  print(paste0("Processing time: ", round(process.time, digits = 4), " mins (current rho = ", rho, ")"))
  return(rst)
}

######----------- Corresponding Plot Functions ------------######
plot.function.bias.RMSE <- function(res.df, label){
  bias.df <- res.df[, c("rho", "bias.beta.pooled", "bias.beta.FE")]
  bias.df <- as.data.frame(bias.df) %>% dplyr::select(rho, starts_with("bias.beta")) %>% 
    pivot_longer(cols = starts_with("bias.beta"), names_to = "model.type", 
                 names_prefix = "bias.beta.", values_to = "bias")
  Bias.plot <- ggplot(bias.df, aes(rho, group = factor(model.type))) +
    geom_line(aes(y = bias, color = factor(model.type), 
                  linetype = factor(model.type)), size = 1)  + 
    theme(panel.grid = element_blank(), panel.background = element_blank(),
          axis.line = element_line(colour = "black")) + 
    theme(axis.title = element_text(size = 24, family = "serif", face = "bold"),
          legend.text = element_text(size = 14, family = "serif")) + 
    theme(axis.text = element_text(face = "bold.italic", size = 17, family = "serif")) + 
    theme(legend.position = c(0.85, 0.88),
          legend.background = element_rect(fill = "white", color = "black"),
          legend.title = element_blank()) + 
    theme(legend.key.height = unit(0.6, 'cm'),
          legend.key.width = unit(1.2, 'cm'),
          legend.spacing.y = unit(0, 'cm')) +
    labs(title = "", 
         x = bquote(rho),
         y = "Bias") +
    scale_x_continuous(breaks = seq(-8, 8, 2)/10) + 
    scale_linetype_manual(values = c("solid", "dotdash"), name = "", 
                          labels = c("FE", "Pooled")) + 
    scale_color_manual(values = c("red", "blue"), name = "", 
                       labels = c("FE", "Pooled")) + ylim(-4, 4)
  
  
  RMSE.df <- res.df[, c("rho", "RMSE.beta.pooled", "RMSE.beta.FE")]
  RMSE.df <- as.data.frame(RMSE.df) %>% dplyr::select(rho, starts_with("RMSE.beta")) %>% 
    pivot_longer(cols = starts_with("RMSE.beta"), names_to = "model.type", 
                 names_prefix = "RMSE.beta.", values_to = "RMSE")
  RMSE.plot <- ggplot(RMSE.df, aes(rho, group = factor(model.type))) +
    geom_line(aes(y = RMSE, color = factor(model.type), linetype = factor(model.type)), size = 0.8)  + 
    theme(panel.grid = element_blank(), panel.background = element_blank(),
          axis.line = element_line(colour = "black")) + 
    theme(axis.title = element_text(size = 24, family = "serif", face = "bold"),
          legend.text = element_text(size = 14, family = "serif")) + 
    theme(axis.text = element_text(face = "bold.italic", size = 17, family = "serif")) + 
    theme(legend.position = c(0.85, 0.88),
          legend.background = element_rect(fill = "white", color = "black"),
          legend.title = element_blank()) + 
    theme(legend.key.height = unit(0.6, 'cm'),
          legend.key.width = unit(1.2, 'cm'),
          legend.spacing.y = unit(0, 'cm')) +
    labs(title = "", 
         x = bquote(rho),
         y = "RMSE") +
    scale_x_continuous(breaks = seq(-8, 8, 2)/10) + 
    scale_linetype_manual(values = c("solid", "dotdash"), name = "", 
                          labels = c("FE", "Pooled")) + 
    scale_color_manual(values = c("red", "blue"), name = "", 
                       labels = c("FE", "Pooled")) + ylim(0, 4.2)
  
  
  
  return.plot <- annotate_figure(ggarrange(Bias.plot, RMSE.plot, 
                                           nrow = 1, ncol = 2,
                                           labels = label))
  
  return(return.plot)
}

plot.function.bias.RMSE.New.m100.1.03 <- plot.function.bias.RMSE(sim1.res.m100, c("A1", "A2"))
plot.function.bias.RMSE.New.m100.2.03 <- plot.function.bias.RMSE(sim3.res.m100, c("B1", "B2"))
plot.function.bias.RMSE.New.m200.1.03 <- plot.function.bias.RMSE(sim1.res.m200, c("C1", "C2"))
plot.function.bias.RMSE.New.m200.2.03 <- plot.function.bias.RMSE(sim3.res.m200, c("D1", "D2"))


save(sim1.res.m100, plot.function.bias.RMSE.New.m100.1.03, 
     sim3.res.m100, plot.function.bias.RMSE.New.m100.2.03,
     sim1.res.m200, plot.function.bias.RMSE.New.m200.1.03,
     sim3.res.m200, plot.function.bias.RMSE.New.m200.2.03,
     file = paste0("Binary_bias.RMSE_NEW_", Sys.Date(), ".RData"))

