library(MASS)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(grid)
library(dplyr)
setwd("/home/ybshao/GroupLasso/Simulations/Figures&Tables/FE_vs_Pooled/Continuous/NEW")

######----------- Simulation Functions (NEW) ------------######
sim.new <- function(ls){
  print(paste0("Iteration for rho = ", ls$rho, " starts.."))
  start.time <- Sys.time()
  m <- ls$m
  sd.z <- ls$sd.z
  sd.gamma <- ls$sd.gamma
  sd.err <- ls$sd.err 
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
      Y <- gamma.rep + Z * beta + rnorm(N, 0, sd.err)
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
    #(1) Pooled model
    pooled <- lm(Y ~ X, data = tb)
    beta_pooled <- coef(pooled)[Z.char]
    
    #(2) Fixed effect model
    fixed <- glm(as.formula(paste(Y.char, "~ 0 + ", prov.char, "+", 
                                  paste0(Z.char, collapse = "+"))), data = tb)
    beta_fe <- coef(fixed)[Z.char]
    c(beta_pooled, beta_fe) 
  }
  
  `%dopar%` <- foreach::`%dopar%`
  cl <- parallel::makeCluster(4)
  doParallel::registerDoParallel(cl)
  system.time(
    mat <- foreach::foreach(i = 1:loop_length, .combine = cbind, .packages = c("Matrix")) %dopar% {
      sim.continuous()})
  parallel::stopCluster(cl)
  
  #get statistics after loops end
  bias.beta <- rowMeans(mat - beta_true)
  bias.beta.pooled <- rowMeans(mat - beta_true)[1]
  ESD.beta <- apply(mat, 1, sd)
  RMSE.beta <- sqrt(bias.beta ^ 2 + ESD.beta ^ 2)
  sd.z.gamma <- paste0("sd.gamma = ", sd.gamma, "; sd.z = ", sd.z)
  rst <- c(m, sd.z.gamma, rho, bias.beta.pooled, RMSE.beta) 
  names(rst) <- c("num.units", "sd.z.gamma",
                  "rho", "bias.beta.pooled",
                  "RMSE.beta.pooled", "RMSE.beta.FE")
  end.time <- Sys.time()
  process.time <- difftime(end.time, start.time, units = 'mins')
  print(paste0("Processing time: ", round(process.time, digits = 4), " mins (current rho = ", rho, ")"))
  return(rst)
}

## 1. 50 units 
# (1) sd.gamma = 1 & sd.z = 0.3
sim.res <- NULL
for (rho in seq(-8, 8, 1)/10){
  ls <- list(m = 50, poisson.mean = 100, rho = rho, n.rep = 1000, sd.gamma = 1, sd.z = 0.3, sd.err = 4)
  sim.res <- rbind(sim.res, sim.new(ls))
}
sim1.res <- as.data.frame(sim.res)

# (2) sd.gamma = 1 & sd.z = 1
sim.res <- NULL
for (rho in seq(-8, 8, 1)/10){
  ls <- list(m = 50, poisson.mean = 100, rho = rho, n.rep = 1000, sd.gamma = 1, sd.z = 1, sd.err = 4)
  sim.res <- rbind(sim.res, sim.new(ls))
}
sim2.res <- as.data.frame(sim.res)

# (3) sd.gamma = 2 & sd.z = 0.3
sim.res <- NULL
for (rho in seq(-8, 8, 1)/10){
  ls <- list(m = 50, poisson.mean = 100, rho = rho, n.rep = 1000, sd.gamma = 2, sd.z = 0.3, sd.err = 4)
  sim.res <- rbind(sim.res, sim.new(ls))
}
sim3.res <- as.data.frame(sim.res)


# (4) sd.gamma = 2 & sd.z = 1
sim.res <- NULL
for (rho in seq(-8, 8, 1)/10){
  ls <- list(m = 50, poisson.mean = 100, rho = rho, n.rep = 1000, sd.gamma = 2, sd.z = 1, sd.err = 4)
  sim.res <- rbind(sim.res, sim.new(ls))
}
sim4.res <- as.data.frame(sim.res)

combine.sim.m50 <- rbind(sim1.res, sim2.res, sim3.res, sim4.res, stringsAsFactors = F)

save(combine.sim.m50, sim1.res, sim2.res, sim3.res, sim4.res, 
     file = paste0("FE_vs_Pooled_Continuous_m50_", Sys.Date(), ".RData"))


## 2. 100 units 
# (1) sd.gamma = 1 & sd.z = 0.3
sim.res <- NULL
for (rho in seq(-8, 8, 1)/10){
  ls <- list(m = 100, poisson.mean = 100, rho = rho, n.rep = 1000, sd.gamma = 1, sd.z = 0.3, sd.err = 4)
  sim.res <- rbind(sim.res, sim.new(ls))
}
sim1.res <- as.data.frame(sim.res)

# (2) sd.gamma = 1 & sd.z = 1
sim.res <- NULL
for (rho in seq(-8, 8, 1)/10){
  ls <- list(m = 100, poisson.mean = 100, rho = rho, n.rep = 1000, sd.gamma = 1, sd.z = 1, sd.err = 4)
  sim.res <- rbind(sim.res, sim.new(ls))
}
sim2.res <- as.data.frame(sim.res)

# (3) sd.gamma = 2 & sd.z = 0.3
sim.res <- NULL
for (rho in seq(-8, 8, 1)/10){
  ls <- list(m = 100, poisson.mean = 100, rho = rho, n.rep = 1000, sd.gamma = 2, sd.z = 0.3, sd.err = 4)
  sim.res <- rbind(sim.res, sim.new(ls))
}
sim3.res <- as.data.frame(sim.res)


# (4) sd.gamma = 2 & sd.z = 1
sim.res <- NULL
for (rho in seq(-8, 8, 1)/10){
  ls <- list(m = 100, poisson.mean = 100, rho = rho, n.rep = 1000, sd.gamma = 2, sd.z = 1, sd.err = 4)
  sim.res <- rbind(sim.res, sim.new(ls))
}
sim4.res <- as.data.frame(sim.res)

combine.sim.m100 <- rbind(sim1.res, sim2.res, sim3.res, sim4.res, stringsAsFactors = F)

save(combine.sim.m100, sim1.res, sim2.res, sim3.res, sim4.res, 
     file = paste0("FE_vs_Pooled_Continuous_m100_", Sys.Date(), ".RData"))

## 3. 200 units 
# (1) sd.gamma = 1 & sd.z = 0.3
sim.res <- NULL
for (rho in seq(-8, 8, 1)/10){
  ls <- list(m = 200, poisson.mean = 100, rho = rho, n.rep = 1000, sd.gamma = 1, sd.z = 0.3, sd.err = 4)
  sim.res <- rbind(sim.res, sim.new(ls))
}
sim1.res <- as.data.frame(sim.res)

# (2) sd.gamma = 1 & sd.z = 1
sim.res <- NULL
for (rho in seq(-8, 8, 1)/10){
  ls <- list(m = 200, poisson.mean = 100, rho = rho, n.rep = 1000, sd.gamma = 1, sd.z = 1, sd.err = 4)
  sim.res <- rbind(sim.res, sim.new(ls))
}
sim2.res <- as.data.frame(sim.res)

# (3) sd.gamma = 2 & sd.z = 0.3
sim.res <- NULL
for (rho in seq(-8, 8, 1)/10){
  ls <- list(m = 200, poisson.mean = 100, rho = rho, n.rep = 1000, sd.gamma = 2, sd.z = 0.3, sd.err = 4)
  sim.res <- rbind(sim.res, sim.new(ls))
}
sim3.res <- as.data.frame(sim.res)


# (4) sd.gamma = 2 & sd.z = 1
sim.res <- NULL
for (rho in seq(-8, 8, 1)/10){
  ls <- list(m = 200, poisson.mean = 100, rho = rho, n.rep = 1000, sd.gamma = 2, sd.z = 1, sd.err = 4)
  sim.res <- rbind(sim.res, sim.new(ls))
}
sim4.res <- as.data.frame(sim.res)

combine.sim.m200 <- rbind(sim1.res, sim2.res, sim3.res, sim4.res, stringsAsFactors = F)

save(combine.sim.m50, combine.sim.m100, combine.sim.m200, 
     file = paste0("Continuous_RMSE_Bias_", Sys.Date(), ".RData"))




######----------- New plot: theoretical FE bias & empirical pooled bias  ------------######
plot.function.bias <- function(res.df, label){
  bias.df <- res.df[, c("sd.z.gamma", "rho", "bias.beta.pooled")]
  bias.df$rho <- as.numeric(as.character(res.df$rho))
  bias.df$bias.beta.pooled <- as.numeric(as.character(res.df$bias.beta.pooled))
  
  Bias.plot <- ggplot(bias.df, aes(rho, group = factor(sd.z.gamma))) +
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
    scale_linetype_manual(values = c("solid", "dashed", "dotdash", "dotted"), name = "", 
                          labels = c(bquote(sigma[gamma] ~ " = 1, " ~ sigma[X] ~ " = 0.3"), 
                                     bquote(sigma[gamma] ~ " = 1, " ~ sigma[X] ~ " = 1"), 
                                     bquote(sigma[gamma] ~ " = 2, " ~ sigma[X] ~ " = 0.3"), 
                                     bquote(sigma[gamma] ~ " = 2, " ~ sigma[X] ~ " = 1"))) + 
    scale_color_manual(values = c("red", "blue", "black", "dark grey"), name = "", 
                       labels = c(bquote(sigma[gamma] ~ " = 1, " ~ sigma[X] ~ " = 0.3"), 
                                  bquote(sigma[gamma] ~ " = 1, " ~ sigma[X] ~ " = 1"), 
                                  bquote(sigma[gamma] ~ " = 2, " ~ sigma[X] ~ " = 0.3"), 
                                  bquote(sigma[gamma] ~ " = 2, " ~ sigma[X] ~ " = 1"))) 
  Bias.plot <- annotate_figure(ggarrange(Bias.plot,
                                         nrow = 1, ncol = 1,
                                         labels = label))
  return(Bias.plot)
}

plot.bias.m50 <- plot.function.bias(combine.sim.m50, "A")
plot.bias.m100 <- plot.function.bias(combine.sim.m100, "B")
plot.bias.m200 <- plot.function.bias(combine.sim.m200, "C")

save(plot.bias.m50, plot.bias.m100, plot.bias.m200, 
     file = paste0("Pooled_bias_figure", Sys.Date(), ".RData"))


plot.function.RMSE <- function(res.df, sd.gamma, sd.z, label){
  sd.z.gamma <- paste0("sd.gamma = ", sd.gamma, "; sd.z = ", sd.z)
  res.df <- res.df[res.df$sd.z.gamma == sd.z.gamma,]
  RMSE.df <- as.data.frame(res.df) %>% dplyr::select(rho, starts_with("RMSE.beta")) %>% 
    pivot_longer(cols = starts_with("RMSE.beta"), names_to = "model.type", 
                 names_prefix = "RMSE.beta.", values_to = "RMSE")
  
  RMSE.df$rho <- as.numeric(as.character(RMSE.df$rho))
  RMSE.df$RMSE <- as.numeric(as.character(RMSE.df$RMSE))
  
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
                       labels = c("FE", "Pooled")) + ylim(0, 6)
  
  RMSE.plot <- annotate_figure(ggarrange(RMSE.plot,
                                         nrow = 1, ncol = 1,
                                         labels = label))
  return(RMSE.plot)
}

plot.RMSE.1.03 <- plot.function.RMSE(combine.sim.m200, 1, 0.3, "A")
plot.RMSE.1.1 <- plot.function.RMSE(combine.sim.m200, 1, 1, "B")
plot.RMSE.2.03 <- plot.function.RMSE(combine.sim.m200, 2, 0.3, "C")
plot.RMSE.2.1 <- plot.function.RMSE(combine.sim.m200, 2, 1, "D")

save(plot.RMSE.1.03, plot.RMSE.1.1, plot.RMSE.2.03, plot.RMSE.2.1,
     file = paste0("Pooled_RMSE_figure_", Sys.Date(), ".RData"))

