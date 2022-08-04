######----------- Pooled model v.s. Fixed effect model ------------######
setwd("/home/ybshao/grplasso/GroupLasso/Simulations/FE_vs_Pooled")
library(MASS)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(grid)

######----------- Simulation Functions ------------######
sim <- function(ls){
  start.time <- Sys.time()
  sd.gamma<-1
  m <- ls$m
  sd.z <- 0.3
  beta <- 1
  beta_true <- beta
  loop_length <- ls$n.rep
  rho <- ls$rho
  sim.continuous <- function() {
    data.proveff <- function(m, prov.size, gamma, sd.gamma, sd.z,
                             beta, Y.char, Z.char, prov.char, rho) {
      N <- sum(prov.size) # total number of discharges
      gamma.rep <- rep(gamma, times = prov.size)
      prov <- rep(paste0("prov", sprintf(paste0("%0",nchar(m),"d"),1:m)),
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
    prov.size <- rep(ls$n.obs, m)
    gamma <- rnorm(m,0,sd.gamma)
    Y.char <- 'Y'
    prov.char <- 'unit'
    Z.char <- "X"
    tb <- data.proveff(m, prov.size, gamma, sd.gamma, sd.z, beta, 
                       Y.char, Z.char, prov.char, rho)
    
    #fit two model and obtain estimators
    pool <- glm(as.formula(paste(Y.char, "~", paste0(Z.char, collapse = "+"))), data = tb, family = "binomial") 
    beta_pool <- unname(coef(pool)[Z.char])
    names(beta_pool) <- NULL
    
    fixed <- glm(as.formula(paste(Y.char, "~ 0 + ", prov.char, "+", paste0(Z.char, collapse = "+"))), data = tb, family = "binomial")
    beta_fe <- unname(coef(fixed)[Z.char])
    
    c(beta_pool, beta_fe) 
  }
  
  `%dopar%` <- foreach::`%dopar%`
  cl <- parallel::makeCluster(4)
  doParallel::registerDoParallel(cl)
  system.time(
    mat <- foreach::foreach(i = 1:loop_length, .combine = cbind) %dopar% {
      sim.continuous()})
  parallel::stopCluster(cl)
  
  #get statistics after loops end
  bias.beta <- rowMeans(mat - beta_true)
  ESD.beta <- apply(mat, 1, sd)
  RMSE.beta <- sqrt(bias.beta ^ 2 + ESD.beta ^ 2)
  Absbias.beta <- rowMeans(abs(mat - beta_true))
  rst <- c(ls$m, ls$n.obs, rho, bias.beta, ESD.beta, RMSE.beta, Absbias.beta) 
  names(rst) <- c("num.units", "num.obs", "rho", 
                  "bias.beta.pooled", "bias.beta.FE",
                  "ESD.beta.pooled", "ESD.beta.FE",
                  "RMSE.beta.pooled", "RMSE.beta.FE",
                  "Absbias.beta.pooled", "Absbias.beta.FE")
  end.time <- Sys.time()
  process.time <- difftime(end.time, start.time, units = 'mins')
  print(paste0("Processing time: ", round(process.time, digits = 4), " mins (current rho = ", rho, ")"))
  return(rst)
}

plot.function <- function(res.df){
  bias.df <- as.data.frame(res.df) %>% select(rho, starts_with("bias.beta")) %>% 
    pivot_longer(cols = starts_with("bias.beta"), names_to = "model.type", 
                 names_prefix = "bias.beta.", values_to = "bias")
  Bias.plot <- ggplot(bias.df, aes(rho, group = factor(model.type))) +
    geom_line(aes(y = bias, color = factor(model.type), linetype = factor(model.type)), size = 0.8)  + 
    theme(panel.grid = element_blank(), panel.background = element_blank(),
          axis.line = element_line(colour = "black")) + 
    theme(axis.title = element_text(size = 13, family = "serif"),
          legend.text = element_text(size = 12, family = "serif")) + 
    theme(axis.text = element_text(face = "italic", size = 9, family = "serif")) + 
    theme(legend.direction = "horizontal", legend.spacing.x = unit(0.5, 'cm')) + 
    theme(legend.key.height = unit(0.5, 'cm'),
          legend.key.width = unit(1, 'cm')) +
    labs(title = "", 
         x = "", 
         y = bquote("Bias of" ~ beta)) +
    scale_x_continuous(breaks = seq(-8, 8, 2)/10) + 
    scale_linetype_manual(values = c("solid", "dotdash"), name = "", labels = c("Fixed Effects", "Pooled")) + 
    scale_color_manual(values = c("red", "blue"), name = "", labels = c("Fixed Effects", "Pooled")) 
  
  ESD.df <- as.data.frame(res.df) %>% select(rho, starts_with("ESD.beta")) %>% 
    pivot_longer(cols = starts_with("ESD.beta"), names_to = "model.type", 
                 names_prefix = "ESD.beta.", values_to = "ESD")
  ESD.plot <- ggplot(ESD.df, aes(rho, group = factor(model.type))) +
    geom_line(aes(y = ESD, color = factor(model.type), linetype = factor(model.type)), size = 0.8)  + 
    theme(panel.grid = element_blank(), panel.background = element_blank(),
          axis.line = element_line(colour = "black")) + 
    theme(axis.title = element_text(size = 13, family = "serif"),
          legend.text = element_text(size = 12, family = "serif")) + 
    theme(axis.text = element_text(face = "italic", size = 9, family = "serif")) + 
    theme(legend.direction = "horizontal", legend.spacing.x = unit(0.5, 'cm')) + 
    theme(legend.key.height = unit(0.5, 'cm'),
          legend.key.width = unit(1, 'cm')) +
    labs(title = "", 
         x = "", 
         y = bquote("ESD of" ~ beta)) +
    scale_x_continuous(breaks = seq(-8, 8, 2)/10) + 
    scale_linetype_manual(values = c("solid", "dotdash"), name = "", labels = c("Fixed Effects", "Pooled")) + 
    scale_color_manual(values = c("red", "blue"), name = "", labels = c("Fixed Effects", "Pooled")) 
  
  RMSE.df <- as.data.frame(res.df) %>% select(rho, starts_with("RMSE.beta")) %>% 
    pivot_longer(cols = starts_with("RMSE.beta"), names_to = "model.type", 
                 names_prefix = "RMSE.beta.", values_to = "RMSE")
  RMSE.plot <- ggplot(RMSE.df, aes(rho, group = factor(model.type))) +
    geom_line(aes(y = RMSE, color = factor(model.type), linetype = factor(model.type)), size = 0.8)  + 
    theme(panel.grid = element_blank(), panel.background = element_blank(),
          axis.line = element_line(colour = "black")) + 
    theme(axis.title = element_text(size = 13, family = "serif"),
          legend.text = element_text(size = 12, family = "serif")) + 
    theme(axis.text = element_text(face = "italic", size = 9, family = "serif")) + 
    theme(legend.direction = "horizontal", legend.spacing.x = unit(0.5, 'cm')) + 
    theme(legend.key.height = unit(0.5, 'cm'),
          legend.key.width = unit(1, 'cm')) +
    labs(title = "", 
         x = "", 
         y = bquote("RMSE of" ~ beta)) +
    scale_x_continuous(breaks = seq(-8, 8, 2)/10) + 
    scale_linetype_manual(values = c("solid", "dotdash"), name = "", labels = c("Fixed Effects", "Pooled")) + 
    scale_color_manual(values = c("red", "blue"), name = "", labels = c("Fixed Effects", "Pooled")) 
  
  Absbias.df <- as.data.frame(res.df) %>% select(rho, starts_with("Absbias.beta")) %>% 
    pivot_longer(cols = starts_with("Absbias.beta"), names_to = "model.type", 
                 names_prefix = "Absbias.beta.", values_to = "Absbias")
  Absbias.plot <- ggplot(Absbias.df, aes(rho, group = factor(model.type))) +
    geom_line(aes(y = Absbias, color = factor(model.type), linetype = factor(model.type)), size = 0.8)  + 
    theme(panel.grid = element_blank(), panel.background = element_blank(),
          axis.line = element_line(colour = "black")) + 
    theme(axis.title = element_text(size = 13, family = "serif"),
          legend.text = element_text(size = 12, family = "serif")) + 
    theme(axis.text = element_text(face = "italic", size = 9, family = "serif")) + 
    theme(legend.direction = "horizontal", legend.spacing.x = unit(0.5, 'cm')) + 
    theme(legend.key.height = unit(0.5, 'cm'),
          legend.key.width = unit(1, 'cm')) +
    labs(title = "", 
         x = "", 
         y = bquote("Absolute bias of" ~ beta)) +
    scale_x_continuous(breaks = seq(-8, 8, 2)/10) + 
    scale_linetype_manual(values = c("solid", "dotdash"), name = "", labels = c("Fixed Effects", "Pooled")) + 
    scale_color_manual(values = c("red", "blue"), name = "", labels = c("Fixed Effects", "Pooled")) 
  
  
  extract_legend <- function(myPlot) {
    tmp <- ggplot_gtable(ggplot_build(myPlot))
    tmp2 <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[tmp2]]
    return(legend)
  }
  shared_legend <- extract_legend(Bias.plot)
  
  title <- paste0(res.df[1,1], " units, with ", res.df[1,2], " observations per unit \n")
  TitleText <- textGrob(title, gp = gpar(fontfamily = "serif", fontsize = 16, fontface = "bold"))
  plts <- annotate_figure(ggarrange(Bias.plot + theme(legend.position = "none"),
                                    ESD.plot + theme(legend.position = "none"),
                                    RMSE.plot + theme(legend.position = "none"),
                                    Absbias.plot + theme(legend.position = "none"),
                                    nrow = 2, ncol = 2, common.legend = T),
                          bottom = text_grob(bquote(rho), size = 14),
                          top = TitleText)
  return(plts)
}


# rho range from -0.6 to +0.6
# provider counts: 20 & 50 & 100
# obs within each provider: 50 & 100

# 1. 20 units & 50 obs
res <- NULL
for (rho in seq(-6, 6, 1)/10){
  ls <- list(m = 20, n.obs = 50, rho = rho, n.rep = 10000)
  res <- rbind(res, sim(ls))
}
res1 <- as.data.frame(res)
figure1 <- plot.function(res1)

# 2. 50 units & 50 obs
res <- NULL
for (rho in seq(-6, 6, 1)/10){
  ls <- list(m = 50, n.obs = 50, rho = rho, n.rep = 10000)
  res <- rbind(res, sim(ls))
}
res2 <- as.data.frame(res)
figure2 <- plot.function(res2)

# 3. 100 units & 50 obs
res <- NULL
for (rho in seq(-6, 6, 1)/10){
  ls <- list(m = 100, n.obs = 50, rho = rho, n.rep = 10000)
  res <- rbind(res, sim(ls))
}
res3 <- as.data.frame(res)
figure3 <- plot.function(res3)

# 4. 20 units & 100 obs
res <- NULL
for (rho in seq(-6, 6, 1)/10){
  ls <- list(m = 20, n.obs = 100, rho = rho, n.rep = 10000)
  res <- rbind(res, sim(ls))
}
res4 <- as.data.frame(res)
figure4 <- plot.function(res4)

# 5. 50 units & 100 obs
res <- NULL
for (rho in seq(-6, 6, 1)/10){
  ls <- list(m = 50, n.obs = 100, rho = rho, n.rep = 10000)
  res <- rbind(res, sim(ls))
}
res5 <- as.data.frame(res)
figure5 <- plot.function(res5)

# 6. 100 units & 100 obs
res <- NULL
for (rho in seq(-6, 6, 1)/10){
  ls <- list(m = 100, n.obs = 100, rho = rho, n.rep = 10000)
  res <- rbind(res, sim(ls))
}
res6 <- as.data.frame(res)
figure6 <- plot.function(res6)


save(res1, figure1, 
     res2, figure2, 
     res3, figure3, 
     res4, figure4, 
     res5, figure5, 
     res6, figure6, 
     file = paste0("FE_vs_Pooled_", Sys.Date(), ".RData"))



