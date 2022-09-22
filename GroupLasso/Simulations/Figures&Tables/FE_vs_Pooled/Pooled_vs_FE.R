######----------- Pooled model v.s. Fixed effect model ------------######
library(MASS)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(grid)
library(dplyr)

######----------- Corresponding Plot Functions ------------######
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
    theme(legend.key = element_blank(),
          legend.key.height = unit(0.5, 'cm'),
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
    theme(legend.key = element_blank(),
          legend.key.height = unit(0.5, 'cm'),
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
    theme(legend.key = element_blank(),
          legend.key.height = unit(0.5, 'cm'),
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
    theme(legend.key = element_blank(),
          legend.key.height = unit(0.5, 'cm'),
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
  
  #title <- paste0(res.df[1, 1], " units, with ", res.df[1,2], " observations per unit \n")
  title <- paste0(res.df[1, 1], " units, with number of observations derived from Poisson(", res.df[1,2],")")
  TitleText <- textGrob(title, gp = gpar(fontfamily = "serif", fontsize = 12, fontface = "bold"))
  plts <- annotate_figure(ggarrange(Bias.plot + theme(legend.position = "none"),
                                    ESD.plot + theme(legend.position = "none"),
                                    RMSE.plot + theme(legend.position = "none"),
                                    Absbias.plot + theme(legend.position = "none"),
                                    nrow = 2, ncol = 2, common.legend = T),
                          bottom = text_grob(bquote(rho), size = 14),
                          top = TitleText)
  return(plts)
}

######----------- Simulation Functions ------------######

######----------- Sim1: fix obs within each unit ------------######
#ls <- list(m = 20, n.obs = 50, rho = rho, n.rep = 10000)
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
    pool <- grp.lasso(tb, Y.char, Z.char, lambda = 0)
    beta_pool <- coef(pool)$beta
    
    fixed <- grp.lasso(tb, Y.char, Z.char, prov.char, lambda = 0)
    beta_fe <- coef(fixed)$beta
    
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
     file = paste0("FE_vs_Pooled_fixed_", Sys.Date(), ".RData"))

######----------- Sim2: obs within each unit derived from a distribution ------------######
# ls <- list(m = 200, poisson.mean = 100, rho = rho, n.rep = 10000, sd.gamma = 1, sd.z = 0.3)
setwd("/home/ybshao/GroupLasso/Simulations/Figures&Tables/FE_vs_Pooled/GLM_Binary/bias_RMSE/From_fixed_number")
sim2 <- function(ls){
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


# rho range from -0.6 to +0.6
# 1. 200 units & poisson(100)
sim2.res <- NULL
for (rho in seq(-6, 6, 1)/10){
  ls <- list(m = 200, poisson.mean = 100, rho = rho, n.rep = 10000, sd.gamma = 1, sd.z = 0.3)
  sim2.res <- rbind(sim2.res, sim2(ls))
}
sim2.res1 <- as.data.frame(sim2.res)
figure.sim2.res1 <- plot.function(sim2.res1)

save(sim2.res1, figure.sim2.res1, 
     file = paste0("FE_vs_Pooled_poisson100_", Sys.Date(), ".RData"))

# 2. 200 units & poisson(200)
sim2.res <- NULL
for (rho in seq(-6, 6, 1)/10){
  ls <- list(m = 200, poisson.mean = 200, rho = rho, n.rep = 10000)
  sim2.res <- rbind(sim2.res, sim2(ls))
}
sim2.res2 <- as.data.frame(sim2.res)
figure.sim2.res2 <- plot.function(sim2.res2)

save(sim2.res2, figure.sim2.res2, 
     file = paste0("FE_vs_Pooled_poisson200_", Sys.Date(), ".RData"))

######----------- Sim3: simulate FDP comparison between pooled model and fixed effect model------------######
# ls <- list(m = 200, poisson.mean = 100, rho = 0.6, n.rep = 4, n.beta = 100, prop.NonZero.beta = 0.05, beta = 1)
sim3 <- function(ls){
  print(paste0("Iteration for beta = ", ls$beta, " starts.."))
  start.time <- Sys.time()
  sd.gamma <- 1
  m <- ls$m
  sd.z <- 1
  true.beta <- c(rep(ls$beta, ls$n.beta * ls$prop.NonZero.beta), 
                 rep(0, ls$n.beta * (1 - ls$prop.NonZero.beta)))  # 5% beta's are non-zero
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
    # cor(rep(gamma, times = prov.size), tb$z1)
    # cor(tb$z1, tb$z2)
    
    
    true.non.zero.beta.index <- 1:(ls$n.beta * ls$prop.NonZero.beta)
    true.zero.beta.index <- (ls$n.beta * ls$prop.NonZero.beta + 1):ls$n.beta
    # (1) Fixed effect model 
    cv.FE.lasso <- cv.grp.lasso(tb, Y.char, Z.char, prov.char) # by default, each covariate forms a group (simple lasso problem)
    cv.FE.fit <- cv.FE.lasso$fit
    FE.beta.under.best.lambda <- cv.FE.fit$beta[, cv.FE.lasso$min]  #beta estimate under best lambda (by 10-fold cross validation)
    FE.estimated.non.zero.beta.index <- which(FE.beta.under.best.lambda != 0)
    FE.estimated.zero.beta.index <- which(FE.beta.under.best.lambda == 0)

    FE.FD <- sum(!(FE.estimated.non.zero.beta.index %in% true.non.zero.beta.index))
    FE.FN <- sum(!(FE.estimated.zero.beta.index %in% true.zero.beta.index))
    FE.FDP <- FE.FD/length(FE.estimated.non.zero.beta.index)
    FE.power <- sum(true.non.zero.beta.index %in% FE.estimated.non.zero.beta.index)/(ls$n.beta * ls$prop.NonZero.beta)
    FE.FD.add.FN <- FE.FD + FE.FN
    
    # (2) Pooled model 
    cv.pooled.lasso <- cv.grp.lasso(tb, Y.char, Z.char) # treat as "single unit"
    cv.pooled.fit <- cv.pooled.lasso$fit
    pooled.beta.under.best.lambda <- cv.pooled.fit$beta[, cv.pooled.lasso$min]  #beta estimate under best lambda (by 10-fold cross validation)
    pooled.estimated.non.zero.beta.index <- which(pooled.beta.under.best.lambda != 0)
    pooled.estimated.zero.beta.index <- which(pooled.beta.under.best.lambda == 0)
    
    pooled.FD <- sum(!(pooled.estimated.non.zero.beta.index %in% true.non.zero.beta.index))
    pooled.FN <- sum(!(pooled.estimated.zero.beta.index %in% true.zero.beta.index))
    pooled.FDP <- pooled.FD/length(pooled.estimated.non.zero.beta.index)
    pooled.power <- sum(true.non.zero.beta.index %in% pooled.estimated.non.zero.beta.index)/(ls$n.beta * ls$prop.NonZero.beta)
    pooled.FD.add.FN <- pooled.FD + pooled.FN
    
    res <- c(FE.FD, FE.FN, FE.FDP, FE.power, FE.FD.add.FN, pooled.FD, pooled.FN, pooled.FDP, pooled.power, pooled.FD.add.FN) 
    return(res)
  }
  
  `%dopar%` <- foreach::`%dopar%`
  cl <- parallel::makeCluster(4)
  doParallel::registerDoParallel(cl)
  system.time(
    mat <- foreach::foreach(i = 1:loop_length, .combine = cbind, .packages = c("TmpGrlasso", "Matrix")) %dopar% {
      sim.continuous()})
  parallel::stopCluster(cl)
  
  #get statistics after loops end
  Mean.Performance <- rowMeans(mat)
  Sd.Performance <- apply(mat, 1, sd)

  rst <- c(ls$m, ls$poisson.mean, rho, ls$beta, Mean.Performance, Sd.Performance) 
  names(rst) <- c("num.units", "possion.mean", "rho", "current.beta",
                  "FE.FD.Mean", "FE.FN.Mean", "FE.FDP.Mean", "FE.power.Mean", "FE.FD.add.FN.Mean", 
                  "pooled.FD.Mean", "pooled.FN.Mean", "pooled.FDP.Mean", "pooled.power.Mean", "pooled.FD.add.FN.Mean",
                  "FE.FD.Sd", "FE.FN.Sd", "FE.FDP.Sd", "FE.power.Sd", "FE.FD.add.FN.Sd", 
                  "pooled.FD.Sd", "pooled.FN.Sd", "pooled.FDP.Sd", "pooled.power.Sd", "pooled.FD.add.FN.Sd")
  end.time <- Sys.time()
  process.time <- difftime(end.time, start.time, units = 'mins')
  print(paste0("Processing time: ", round(process.time, digits = 4), " mins"))
  return(rst)
}

sim3.res <- NULL
for (beta in seq(3, 10, 1)/10){
  ls <- list(m = 200, poisson.mean = 100, rho = 0.6, n.rep = 200, n.beta = 100, prop.NonZero.beta = 0.05, beta = beta)
  sim3.res <- rbind(sim3.res, sim3(ls))
}

save(sim3.res, file = paste0("/home/ybshao/GroupLasso/Simulations/FE_vs_Pooled/FE_vs_Pooled_FDP_table_", Sys.Date(), ".RData"))
write.csv(sim3.res, file = paste0("/home/ybshao/GroupLasso/Simulations/FE_vs_Pooled/FE_vs_Pooled_FDP_table_", Sys.Date(), ".csv"), quote = F, row.names = T)

######----------- Sim4: simulate lambda v.s FDP -----------######
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
    FE.lasso.lambda <- grp.lasso(tb, Y.char, Z.char, prov.char, lambda = lambda)
    FE.beta <- FE.lasso.lambda$beta
    FE.estimated.non.zero.beta.index <- which(FE.beta != 0)
    FE.estimated.zero.beta.index <- which(FE.beta == 0)
  
    FE.FDP <- sum(!(FE.estimated.non.zero.beta.index %in% true.non.zero.beta.index))/length(FE.estimated.non.zero.beta.index)
    FE.power <- sum(true.non.zero.beta.index %in% FE.estimated.non.zero.beta.index)/(ls$n.beta * ls$prop.NonZero.beta)
    
    # (2) Pooled model 
    pooled.lasso.lambda <- grp.lasso(tb, Y.char, Z.char, lambda = lambda)
    pooled.beta <- pooled.lasso.lambda$beta
    pooled.estimated.non.zero.beta.index <- which(pooled.beta != 0)
    pooled.estimated.zero.beta.index <- which(pooled.beta == 0)
    
    pooled.FDP <- sum(!(pooled.estimated.non.zero.beta.index %in% true.non.zero.beta.index))/length(pooled.estimated.non.zero.beta.index)
    pooled.power <- sum(true.non.zero.beta.index %in% pooled.estimated.non.zero.beta.index)/(ls$n.beta * ls$prop.NonZero.beta)
    
    res <- c(FE.FDP, FE.power, pooled.FDP, pooled.power) 
    return(res)
  }
  

  `%dopar%` <- foreach::`%dopar%`
  cl <- parallel::makeCluster(4)
  doParallel::registerDoParallel(cl)
  system.time(
    mat <- foreach::foreach(i = 1:loop_length, .combine = cbind, .packages = c("TmpGrlasso", "Matrix")) %dopar% {
      sim.continuous()})
  parallel::stopCluster(cl)
  
  Mean.Performance <- rowMeans(mat)
  Sd.Performance <- apply(mat, 1, sd)
  
  rst <- c(ls$m, ls$poisson.mean, rho, ls$lambda, Mean.Performance, Sd.Performance) 
  
  names(rst) <- c("num.units", "possion.mean", "rho", "lambda",
                  "FE.FDP.Mean", "FE.power.Mean", "pooled.FDP.Mean", "pooled.power.Mean",
                  "FE.FDP.Sd", "FE.power.Sd", "pooled.FDP.Sd", "pooled.power.Sd")
  end.time <- Sys.time()
  process.time <- difftime(end.time, start.time, units = 'mins')
  print(paste0("Processing time: ", round(process.time, digits = 4), " mins"))
  return(rst)
}

shared.lambda.seq <- exp(seq(log(1e-01), log(1e-04), length = 100))  #fixed lambda sequence, from 0.1 to 0.0001

sim4.res <- NULL
for (lambda in shared.lambda.seq){
  ls <- list(m = 200, poisson.mean = 100, rho = 0.6, n.rep = 200, n.beta = 100, prop.NonZero.beta = 0.05, beta = 1, lambda = lambda)
  sim4.res <- rbind(sim4.res, sim4(ls))
}
sim4.res <- as.data.frame(sim4.res)

FE.FDP.df <- cbind(log(sim4.res["lambda"]), sim4.res["FE.FDP.Mean"], sim4.res["FE.FDP.Mean"] - sim4.res["FE.FDP.Sd"], 
                   sim4.res["FE.FDP.Mean"] + sim4.res["FE.FDP.Sd"])
colnames(FE.FDP.df) <- NA
pooled.FDP.df <- cbind(log(sim4.res["lambda"]), sim4.res["pooled.FDP.Mean"], sim4.res["pooled.FDP.Mean"] - sim4.res["pooled.FDP.Sd"], 
                       sim4.res["pooled.FDP.Mean"] + sim4.res["pooled.FDP.Sd"])
colnames(pooled.FDP.df) <- NA
FDP.plot.df <- rbind(FE.FDP.df, pooled.FDP.df)
colnames(FDP.plot.df) <- c("lambda", "FD.Mean", "FD.lower", "FD.upper")
FDP.plot.df$model <- rep(c("FE", "Pooled"), each = length(shared.lambda.seq))

FDP.plot <- ggplot(FDP.plot.df, aes(lambda, group = factor(model))) +
  geom_ribbon(aes(ymin = FD.lower, ymax = FD.upper, fill = factor(model)), alpha = 0.5) +
  geom_line(aes(y = FD.Mean, color = factor(model), linetype = factor(model)), size = 1)  + 
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(size = 13, face = "bold", family = "serif"),
        axis.title = element_text(size = 13, family = "serif"),
        plot.caption = element_text(size = 10, face = "italic", family = "serif"),
        legend.text = element_text(size = 13, family = "serif", face = "bold")) + 
  theme(axis.text = element_text(face = "italic", size = 10, family = "serif")) + 
  theme(legend.position = c(0.80, 0.18)) + 
  theme(legend.key.height= unit(0.5, 'cm'),
        legend.key.width= unit(1.5, 'cm')) +
  labs(title = "Porportion of False Discoveries (FDP) of Fixed Effect Model and Pooled Model", 
       x = expression(log(lambda)), 
       y = "FDP") +
  scale_linetype_manual(values = c("solid", "dotdash"), name = "", labels = c("Fixed Effect", "Pooled")) + 
  scale_color_manual(values = c("blue", "red"), name = "", labels = c("Fixed Effect", "Pooled")) + 
  scale_fill_manual(values = c("grey85", "grey90"), name = "", labels = c("Fixed Effect", "Pooled")) +
  scale_x_continuous(trans = scales::reverse_trans(), breaks = round(seq(round(max(log(shared.lambda.seq)), 0), round(min(log(shared.lambda.seq)), 0), by = - 1), 1))

save(FDP.plot.df, FDP.plot,
     file = paste0("/home/ybshao/GroupLasso/Simulations/FE_vs_Pooled/FE_vs_Pooled_lambda_vs_FDP_figure_", Sys.Date(), ".RData"))



######----------- Sim5: Continuous Outcome; From Poisson(100); low dimensional -----------######
# ls <- list(m = 200, poisson.mean = 100, rho = rho, n.rep = 10000, sd.gamma = 1, sd.z = 0.3, sd.err = 4)
setwd("/home/ybshao/GroupLasso/Simulations/Figures&Tables/FE_vs_Pooled/Continuous")
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