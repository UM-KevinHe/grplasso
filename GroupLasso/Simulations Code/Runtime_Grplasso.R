library(MASS)
library(RcppArmadillo)
library(Rcpp)
library(Matrix)
library(grpreg)
library(fastDummies)
library(foreach)
library(doParallel)
library(glmnet)

setwd("/home/ybshao/GroupLasso/Simulations/Figures&Tables/Runtime/GroupLasso")

# n.beta = 50, n.groups = 4, n.provider range from 50 to 500
multiResultClass <- function(Runtime = NULL){#, RME = NULL, RMSE = NULL, cross_entropy = NULL, wrong_prediction_rate = NULL){
  result <- list(Runtime = Runtime)
  ## Set the name for the class
  class(result) <- append(class(result), "multiResultClass")
  return(result)
}


m.sequence <- seq(50, 400, 50)

Runtime.Mean <- matrix(rep(0, length(m.sequence) * 2), nrow = 2)
colnames(Runtime.Mean) <- paste0("n.prov = ", m.sequence)
rownames(Runtime.Mean) <- c("grLasso", "grpreg")

Runtime.sd <- matrix(rep(0, length(m.sequence) * 2), nrow = 2)
colnames(Runtime.sd) <- paste0("n.prov = ", m.sequence)
rownames(Runtime.sd) <- c("grLasso", "grpreg")

data.loop <- 1:10
ind <- 0

for (j in m.sequence){ #outer loop for
  ind <- ind + 1
  cl <- makeCluster(5)
  registerDoParallel(cl) 
  Model.Comparison <- 
    foreach (i = data.loop, .packages = c("grpreg", "fastDummies", "RcppArmadillo", "MASS", "Matrix", "TmpLasso")) %dopar% {
      Y.char <- 'Y'
      prov.char <- 'Prov.ID'
      sim.parameters.GrLasso <- list(m = j, n.beta = 50, n.groups = 10, prop.NonZero.group = 0.5, 
                                     prop.outlier = 0.05, rho = 0.7)
      Sim_GrLasso <- Simulation_data_GroupLasso(sim.parameters.GrLasso, unpenalized.beta = F)
      data_GrLasso <- Sim_GrLasso$sim.data
      group <- Sim_GrLasso$group
      
      true.beta <- Sim_GrLasso$beta  #true beta's for computing RMSE
      Z.char <- paste0('Z_', 1:sim.parameters.GrLasso$n.beta)
      data_prep <- fe.data.prep(data_GrLasso, Y.char, Z.char, prov.char, cutoff = 0, check = FALSE)
      data_prep <- data_prep[data_prep$included == 1, ]
      #true.mu <- data_prep$mu #true mu's for computing RME
      
      # results from grLasso
      start <- Sys.time()
      cv.model_grp_lasso <- cv.grp.lasso(data_prep, Y.char, Z.char, prov.char, group = group, trace.lambda = F,
                                         nfolds = 10, trace.cv = F)
      end <- Sys.time()
      cv.process.time1 <- difftime(end, start, units = 'mins')  #runtime
      
      ## results from grpreg
      start <- Sys.time()
      # dummy data for grpreg
      dummy_data <- dummy_cols(data_prep, select_columns = prov.char, remove_selected_columns = TRUE, 
                               remove_first_dummy = TRUE)
      ID.char <- rep(NA, sim.parameters.GrLasso$m - 1)
      for (i in 1:(sim.parameters.GrLasso$m - 1)){
        ID.char[i] <- paste0("Prov.ID_", i + 1)
      }
      cv.model_grpreg <- cv.grpreg(dummy_data[,c(Z.char, ID.char)], dummy_data[,Y.char], family = "binomial",
                                   penalty = "grLasso", group = c(group, rep(0, length(ID.char))), alpha = 1, 
                                   nfolds = 10, trace.cv = F)
      end <- Sys.time()
      cv.process.time2 <-  difftime(end, start, units = 'mins') #runtime
      result <- multiResultClass()
      result$Runtime <- round(matrix(c(cv.process.time1, cv.process.time2), nrow = 2), digits = 3)
      return(result)
    }
  stopCluster(cl)
  
  n.data.loop <- length(data.loop)
  
  Runtime <- matrix(rep(0, 2 * n.data.loop), nrow = 2)
  rownames(Runtime) <- c("grLasso", "grpreg")
  colnames(Runtime) <- paste0("Data_", data.loop)
  
  for (i in data.loop){
    Runtime[, i] <- Model.Comparison[[i]]$Runtime
  }
  
  Runtime.Mean[, ind] <- round(apply(Runtime, 1, mean), digits = 3)
  Runtime.sd[, ind] <- round(apply(Runtime, 1, sd), digits = 3)
}



#### Figure1: Runtime
library(reshape2)
n.prov <- m.sequence
Runtime.Mean.lower <- Runtime.Mean - Runtime.sd
Runtime.Mean.upper <- Runtime.Mean + Runtime.sd


Runtime.figure.grplasso.df <- data.frame("model" = melt(Runtime.Mean)$Var1,
                                         "n.prov" = rep(n.prov, each = 2),
                                         "Mean" = melt(Runtime.Mean)$value,
                                         "lower" = melt(Runtime.Mean.lower)$value,
                                         "upper" = melt(Runtime.Mean.upper)$value)

Runtime.plot.grplasso <- ggplot(Runtime.figure.grplasso.df, aes(n.prov, group = factor(model))) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(model)), alpha = 0.5) +
  geom_line(aes(y = Mean, color = factor(model), linetype = factor(model)), size = 1)  + 
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(size = 13, face = "bold", family = "serif"),
        axis.title = element_text(size = 14, family = "serif", face = "bold"),
        plot.caption = element_text(size = 10, face = "italic", family = "serif"),
        legend.text = element_text(size = 13, family = "serif")) + 
  theme(axis.text = element_text(face = "italic", size = 11, family = "serif")) + 
  theme(legend.position = c(0.22, 0.82),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank()) + 
  theme(legend.key.height = unit(0.6, 'cm'),
        legend.key.width = unit(1.2, 'cm'),
        legend.spacing.y = unit(0, 'cm')) +
  labs(title = "", 
       x = "Number of providers", 
       y = "Time to convergence (Minutes)") +
  scale_x_continuous(breaks = seq(50, 400, 50)) + 
  scale_linetype_manual(values = c("solid", "dotdash"), name = "", labels = c("grlasso", "grpreg")) + 
  scale_color_manual(values = c("blue", "red"), name = "", labels = c("grlasso", "grpreg")) + 
  scale_fill_manual(values = c("grey85", "grey90"), name = "", labels = c("grlasso", "grpreg"))

Runtime.plot.grplasso <- annotate_figure(ggarrange(Runtime.plot.grplasso, 
                                                   nrow = 1, ncol = 1,
                                                   labels = c("B")))

save(Runtime.figure.grplasso.df, Runtime.plot.grplasso, 
     file = paste0("Runtime_grplasso_", Sys.Date(), ".RData"))

