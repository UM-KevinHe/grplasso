library(MASS)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(grid)
library(dplyr)
setwd("/home/ybshao/GroupLasso/Simulations/Figures&Tables/FE_vs_Pooled/Continuous/NEW/New")

######----------- New plot: theoretical FE bias & empirical pooled bias  ------------######
plot.function.bias <- function(res.df, label, remove = "sd.gamma = 1; sd.z = 1"){
  bias.df <- res.df[, c("sd.z.gamma", "rho", "bias.beta.pooled")]
  bias.df <- bias.df[which(bias.df$sd.z.gamma != remove), ]
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
                          labels = c(bquote(sigma[gamma] ~ " = 1, " ~ sigma[X] ~ " = 0.3"), 
                                     #bquote(sigma[gamma] ~ " = 1, " ~ sigma[X] ~ " = 1"), 
                                     bquote(sigma[gamma] ~ " = 2, " ~ sigma[X] ~ " = 0.3"), 
                                     bquote(sigma[gamma] ~ " = 2, " ~ sigma[X] ~ " = 1"))) + 
    scale_color_manual(values = c("red", "blue", "darkorange"), name = "", 
                       labels = c(bquote(sigma[gamma] ~ " = 1, " ~ sigma[X] ~ " = 0.3"), 
                                  #bquote(sigma[gamma] ~ " = 1, " ~ sigma[X] ~ " = 1"), 
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

plot.RMSE.2.03 <- plot.function.RMSE(combine.sim.m200, 2, 0.3, "A")
plot.RMSE.1.03 <- plot.function.RMSE(combine.sim.m200, 1, 0.3, "B")
plot.RMSE.2.1 <- plot.function.RMSE(combine.sim.m200, 2, 1, "C")

save(plot.RMSE.1.03, plot.RMSE.2.03, plot.RMSE.2.1,
     file = paste0("Pooled_RMSE_figure_", Sys.Date(), ".RData"))