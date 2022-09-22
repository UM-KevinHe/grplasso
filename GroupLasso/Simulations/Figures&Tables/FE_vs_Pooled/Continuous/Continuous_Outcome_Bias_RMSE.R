library(MASS)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(grid)
library(dplyr)
setwd("/home/ybshao/GroupLasso/Simulations/Figures&Tables/FE_vs_Pooled/Continuous")

######----------- Simulation Functions ------------######
# ls <- list(m = 200, poisson.mean = 100, rho = rho, n.rep = 10000, sd.gamma = 1, sd.z = 0.3, sd.err = 4)
sim5 <- function(ls){
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

# rho range from -0.8 to +0.8
# 1. 200 units & poisson(100) & sd.gamma = 1 & sd.z = 0.3
sim.res <- NULL
for (rho in seq(-8, 8, 1)/10){
  ls <- list(m = 200, poisson.mean = 100, rho = rho, n.rep = 10000, sd.gamma = 1, sd.z = 0.3, sd.err = 4)
  sim.res <- rbind(sim.res, sim5(ls))
}
sim1.res <- as.data.frame(sim.res)

# 2. 200 units & poisson(100) & sd.gamma = 1 & sd.z = 1
sim.res <- NULL
for (rho in seq(-8, 8, 1)/10){
  ls <- list(m = 200, poisson.mean = 100, rho = rho, n.rep = 10000, sd.gamma = 1, sd.z = 1, sd.err = 4)
  sim.res <- rbind(sim.res, sim5(ls))
}
sim2.res <- as.data.frame(sim.res)

# 3. 200 units & poisson(100) & sd.gamma = 2 & sd.z = 0.3
sim.res <- NULL
for (rho in seq(-8, 8, 1)/10){
  ls <- list(m = 200, poisson.mean = 100, rho = rho, n.rep = 10000, sd.gamma = 2, sd.z = 0.3, sd.err = 4)
  sim.res <- rbind(sim.res, sim5(ls))
}
sim3.res <- as.data.frame(sim.res)


# 4. 200 units & poisson(100) & sd.gamma = 2 & sd.z = 1
sim.res <- NULL
for (rho in seq(-8, 8, 1)/10){
  ls <- list(m = 200, poisson.mean = 100, rho = rho, n.rep = 10000, sd.gamma = 2, sd.z = 1, sd.err = 4)
  sim.res <- rbind(sim.res, sim5(ls))
}
sim4.res <- as.data.frame(sim.res)

save(sim1.res, sim2.res, sim3.res, sim4.res, 
     file = paste0("FE_vs_Pooled_poisson100_Continuous_", Sys.Date(), ".RData"))


######----------- Sim5-Corresponding Plot Functions ------------######
plot.function.bias <- function(res.df, sd.gamma, sd.z){
  res.df$bias.beta.pooled.Theore <- seq(-8, 8, 1)/10 * sd.gamma / sd.z
  res.df$bias.beta.FE.Theore <- 0
  bias.df <- as.data.frame(res.df) %>% dplyr::select(rho, starts_with("bias.beta")) %>% 
    pivot_longer(cols = starts_with("bias.beta"), names_to = "model.type", 
                 names_prefix = "bias.beta.", values_to = "bias")
  
  Bias.plot <- ggplot(bias.df, aes(rho, group = factor(model.type))) +
    geom_line(aes(y = bias, color = factor(model.type), linetype = factor(model.type), size = factor(model.type)))  + 
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
         caption = bquote("(" ~ sigma[epsilon] ~ " = 4; " ~ sigma[gamma] ~ " = " ~ .(sd.gamma) ~"; " 
                          ~ sigma[Z] ~ " = " ~ .(sd.z) ~ ")")) +
    scale_x_continuous(breaks = seq(-8, 8, 2)/10) + 
    scale_linetype_manual(values = c("solid", "solid", "dotdash", "dotdash"), name = "", 
                          labels = c("FE: Empirical", "FE: Theoretical", "Pooled: Empirical", "Pooled: Theoretical")) + 
    scale_color_manual(values = c("red", "blue", "red", "blue"), name = "", 
                       labels = c("FE: Empirical", "FE: Theoretical", "Pooled: Empirical", "Pooled: Theoretical")) +
    scale_size_manual(values = c(1, 0.5, 1, 0.5), name = "", 
                      labels = c("FE: Empirical", "FE: Theoretical", "Pooled: Empirical", "Pooled: Theoretical"))
  
  return(Bias.plot)
}

plot1.bias <- plot.function.bias(sim1.res, 1, 0.3)
plot2.bias <- plot.function.bias(sim2.res, 1, 1)
plot3.bias <- plot.function.bias(sim3.res, 2, 0.3)
plot4.bias <- plot.function.bias(sim4.res, 2, 1)

extract_legend <- function(myPlot) {
  tmp <- ggplot_gtable(ggplot_build(myPlot))
  tmp2 <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[tmp2]]
  return(legend)
}

shared_legend <- extract_legend(plot1.bias)

plts.Bias <- annotate_figure(ggarrange(plot1.bias + theme(legend.position = "none"),
                                       plot2.bias + theme(legend.position = "none"),
                                       plot3.bias + theme(legend.position = "none"),
                                       plot4.bias + theme(legend.position = "none"),
                                       nrow = 2, ncol = 2, common.legend = T, 
                                       labels = c("A", "B", "C", "D")),
                             bottom = text_grob(bquote(rho), size = 14),
                             left = text_grob(bquote("Bias of" ~ beta), face = "bold", 
                                              family = "serif", size = 14, rot = 90))


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
         caption = bquote("(" ~ sigma[epsilon] ~ " = 4; " ~ sigma[gamma] ~ " = " ~ .(sd.gamma) ~"; " 
                          ~ sigma[Z] ~ " = " ~ .(sd.z) ~ ")")) +
    scale_x_continuous(breaks = seq(-8, 8, 2)/10) + 
    scale_linetype_manual(values = c("solid", "dotdash"), name = "", 
                          labels = c("FE", "Pooled")) + 
    scale_color_manual(values = c("red", "blue"), name = "", 
                       labels = c("FE", "Pooled"))
  return(RMSE.plot)
}

plot1.RMSE <- plot.function.RMSE(sim1.res, 1, 0.3)
plot2.RMSE <- plot.function.RMSE(sim2.res, 1, 1)
plot3.RMSE <- plot.function.RMSE(sim3.res, 2, 0.3)
plot4.RMSE <- plot.function.RMSE(sim4.res, 2, 1)

extract_legend <- function(myPlot) {
  tmp <- ggplot_gtable(ggplot_build(myPlot))
  tmp2 <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[tmp2]]
  return(legend)
}

shared_legend <- extract_legend(plot1.RMSE)

plt.RMSE <- annotate_figure(ggarrange(plot1.RMSE + theme(legend.position = "none"),
                                      plot2.RMSE + theme(legend.position = "none"),
                                      plot3.RMSE + theme(legend.position = "none"),
                                      plot4.RMSE + theme(legend.position = "none"),
                                      nrow = 2, ncol = 2, common.legend = T, 
                                      labels = c("A", "B", "C", "D")),
                            bottom = text_grob(bquote(rho), size = 14),
                            left = text_grob(bquote("Empirical RMSE of" ~ beta), face = "bold", 
                                             family = "serif", size = 14, rot = 90))
