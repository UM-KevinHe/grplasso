library(MASS)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(grid)
library(dplyr)

setwd("/home/ybshao/GroupLasso/Simulations/Figures&Tables/FE_vs_Pooled/GLM_Binary/bias_RMSE/From_poisson")

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
plot.function.bias <- function(res.df, sd.gamma, sd.z){
  bias.df <- as.data.frame(res.df) %>% select(rho, starts_with("bias.beta")) %>% 
    pivot_longer(cols = starts_with("bias.beta"), names_to = "model.type", 
                 names_prefix = "bias.beta.", values_to = "bias")
  
  Bias.plot <- ggplot(bias.df, aes(rho, group = factor(model.type))) +
    geom_line(aes(y = bias, color = factor(model.type), linetype = factor(model.type)), size = 0.8)  + 
    theme(panel.grid = element_blank(), panel.background = element_blank(),
          axis.line = element_line(colour = "black")) + 
    theme(axis.title = element_text(size = 13, family = "serif"),
          legend.text = element_text(size = 10, family = "serif")) + 
    theme(axis.text = element_text(face = "italic", size = 11, family = "serif")) + 
    theme(legend.position = c(0.22, 0.82),
          legend.background = element_rect(fill = "white", color = "black"),
          legend.title = element_blank()) + 
    theme(legend.key.height = unit(0.6, 'cm'),
          legend.key.width = unit(1.2, 'cm'),
          legend.spacing.y = unit(0, 'cm')) +
    labs(title = "", 
         x = "",
         y = "",
         caption = bquote("(" ~ sigma[gamma] ~ " = " ~ .(sd.gamma) ~"; " 
                          ~ sigma[Z] ~ " = " ~ .(sd.z) ~ ")")) +
    scale_x_continuous(breaks = seq(-7, 7, 2)/10) + 
    scale_linetype_manual(values = c("solid", "dotdash"), name = "", labels = c("Fixed Effects", "Pooled")) + 
    scale_color_manual(values = c("red", "blue"), name = "", labels = c("Fixed Effects", "Pooled")) 
  
  return(Bias.plot)
}

plot.function.RMSE <- function(res.df, sd.gamma, sd.z){
  RMSE.df <- as.data.frame(res.df) %>% dplyr::select(rho, starts_with("RMSE.beta")) %>% 
    pivot_longer(cols = starts_with("RMSE.beta"), names_to = "model.type", 
                 names_prefix = "RMSE.beta.", values_to = "RMSE")
  
  RMSE.plot <- ggplot(RMSE.df, aes(rho, group = factor(model.type))) +
    geom_line(aes(y = RMSE, color = factor(model.type), linetype = factor(model.type)), size = 0.8)  + 
    theme(panel.grid = element_blank(), panel.background = element_blank(),
          axis.line = element_line(colour = "black")) + 
    theme(axis.title = element_text(size = 13, family = "serif"),
          legend.text = element_text(size = 10, family = "serif")) + 
    theme(axis.text = element_text(face = "italic", size = 11, family = "serif")) + 
    theme(legend.position = c(0.22, 0.82),
          legend.background = element_rect(fill = "white", color = "black"),
          legend.title = element_blank()) + 
    theme(legend.key.height = unit(0.6, 'cm'),
          legend.key.width = unit(1.2, 'cm'),
          legend.spacing.y = unit(0, 'cm')) +
    labs(title = "", 
         x = "",
         y = "",
         caption = bquote("(" ~ sigma[gamma] ~ " = " ~ .(sd.gamma) ~"; " 
                          ~ sigma[Z] ~ " = " ~ .(sd.z) ~ ")")) +
    scale_x_continuous(breaks = seq(-7, 7, 2)/10) + 
    scale_linetype_manual(values = c("solid", "dotdash"), name = "", 
                          labels = c("Fixed Effects", "Pooled")) + 
    scale_color_manual(values = c("red", "blue"), name = "", 
                       labels = c("Fixed Effects", "Pooled"))
  return(RMSE.plot)
}




######----------- Simulations: m = 200 ------------######
# 1. 200 units & poisson(100) & sd.gamma = 1 & sd.z = 0.3;
sim.res <- NULL
for (rho in seq(-7, 7, 1)/10){
  ls <- list(m = 200, poisson.mean = 100, rho = rho, n.rep = 10000, sd.gamma = 1, sd.z = 0.3)
  sim.res <- rbind(sim.res, sim(ls))
}
sim1.res.m200 <- as.data.frame(sim.res)

# 2. 200 units & poisson(100) & sd.gamma = 1 & sd.z = 1;
sim.res <- NULL
for (rho in seq(-7, 7, 1)/10){
  ls <- list(m = 200, poisson.mean = 100, rho = rho, n.rep = 10000, sd.gamma = 1, sd.z = 1)
  sim.res <- rbind(sim.res, sim(ls))
}
sim2.res.m200 <- as.data.frame(sim.res)

# 3. 200 units & poisson(100) & sd.gamma = 2 & sd.z = 0.3;
sim.res <- NULL
for (rho in seq(-7, 7, 1)/10){
  ls <- list(m = 200, poisson.mean = 100, rho = rho, n.rep = 10000, sd.gamma = 2, sd.z = 0.3)
  sim.res <- rbind(sim.res, sim(ls))
}
sim3.res.m200 <- as.data.frame(sim.res)

plot1.bias <- plot.function.bias(sim1.res.m200, 1, 0.3)
plot2.bias <- plot.function.bias(sim2.res.m200, 1, 1)
plot3.bias <- plot.function.bias(sim3.res.m200, 2, 0.3)

plot1.RMSE <- plot.function.RMSE(sim1.res.m200, 1, 0.3)
plot2.RMSE <- plot.function.RMSE(sim2.res.m200, 1, 1)
plot3.RMSE <- plot.function.RMSE(sim3.res.m200, 2, 0.3)

extract_legend <- function(myPlot) {
  tmp <- ggplot_gtable(ggplot_build(myPlot))
  tmp2 <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[tmp2]]
  return(legend)
}

shared_legend <- extract_legend(plot1.RMSE)

plt.m200 <- annotate_figure(ggarrange(plot1.bias + theme(legend.position = "none"), plot1.RMSE + theme(legend.position = "none"),
                                      plot2.bias + theme(legend.position = "none"), plot2.RMSE + theme(legend.position = "none"),
                                      plot3.bias + theme(legend.position = "none"), plot3.RMSE + theme(legend.position = "none"),
                                      nrow = 3, ncol = 2, common.legend = T, 
                                      labels = c("A.1", "A.2", "B.1", "B.2", "C.1", "C.2")),
                            bottom = text_grob(bquote(rho), size = 14),
                            left = text_grob(bquote("Empirical Bias of" ~ beta), face = "bold", 
                                             family = "serif", size = 14, rot = 90),
                            right = text_grob(bquote("Empirical RMSE of" ~ beta), face = "bold", 
                                              family = "serif", size = 14, rot = 90))
# width = 700; Height = 900


save(sim1.res.m200, sim2.res.m200, sim3.res.m200, plt.m200,
     file = paste0("FE_vs_Pooled_m200_Binary_", Sys.Date(), ".RData"))

######----------- Simulations: m = 100 ------------######
# 1. 100 units & poisson(100) & sd.gamma = 1 & sd.z = 0.3;
sim.res <- NULL
for (rho in seq(-7, 7, 1)/10){
  ls <- list(m = 100, poisson.mean = 100, rho = rho, n.rep = 10000, sd.gamma = 1, sd.z = 0.3)
  sim.res <- rbind(sim.res, sim(ls))
}
sim1.res.m100 <- as.data.frame(sim.res)

# 2. 100 units & poisson(100) & sd.gamma = 1 & sd.z = 1;
sim.res <- NULL
for (rho in seq(-7, 7, 1)/10){
  ls <- list(m = 200, poisson.mean = 100, rho = rho, n.rep = 10000, sd.gamma = 1, sd.z = 1)
  sim.res <- rbind(sim.res, sim(ls))
}
sim2.res.m100 <- as.data.frame(sim.res)

# 3. 200 units & poisson(100) & sd.gamma = 2 & sd.z = 0.3;
sim.res <- NULL
for (rho in seq(-7, 7, 1)/10){
  ls <- list(m = 200, poisson.mean = 100, rho = rho, n.rep = 10000, sd.gamma = 2, sd.z = 0.3)
  sim.res <- rbind(sim.res, sim(ls))
}
sim3.res.m100 <- as.data.frame(sim.res)

plot1.bias <- plot.function.bias(sim1.res.m100, 1, 0.3)
plot2.bias <- plot.function.bias(sim2.res.m100, 1, 1)
plot3.bias <- plot.function.bias(sim3.res.m100, 2, 0.3)

plot1.RMSE <- plot.function.RMSE(sim1.res.m100, 1, 0.3)
plot2.RMSE <- plot.function.RMSE(sim2.res.m100, 1, 1)
plot3.RMSE <- plot.function.RMSE(sim3.res.m100, 2, 0.3)

extract_legend <- function(myPlot) {
  tmp <- ggplot_gtable(ggplot_build(myPlot))
  tmp2 <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[tmp2]]
  return(legend)
}

shared_legend <- extract_legend(plot1.RMSE)

plt.m100 <- annotate_figure(ggarrange(plot1.bias + theme(legend.position = "none"), plot1.RMSE + theme(legend.position = "none"),
                                      plot2.bias + theme(legend.position = "none"), plot2.RMSE + theme(legend.position = "none"),
                                      plot3.bias + theme(legend.position = "none"), plot3.RMSE + theme(legend.position = "none"),
                                      nrow = 3, ncol = 2, common.legend = T, 
                                      labels = c("A.1", "A.2", "B.1", "B.2", "C.1", "C.2")),
                            bottom = text_grob(bquote(rho), size = 14),
                            left = text_grob(bquote("Empirical Bias of" ~ beta), face = "bold", 
                                             family = "serif", size = 14, rot = 90),
                            right = text_grob(bquote("Empirical RMSE of" ~ beta), face = "bold", 
                                              family = "serif", size = 14, rot = 90))

save(sim1.res.m100, sim2.res.m100, sim3.res.m100, plt.m100,
     file = paste0("FE_vs_Pooled_m100_Binary_", Sys.Date(), ".RData"))
