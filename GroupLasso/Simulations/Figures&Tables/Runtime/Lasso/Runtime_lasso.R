library(MASS)
library(RcppArmadillo)
library(Rcpp)
library(Matrix)
library(grpreg)
library(fastDummies)
library(foreach)
library(doParallel)
library(glmnet)

setwd("/home/ybshao/GroupLasso/Simulations/Figures&Tables/Runtime/Lasso")

multiResultClass2 <- function(Runtime = NULL, total_iter = NULL){
  result <- list(
    Runtime = Runtime,
    total_iter = total_iter
  )
  ## Set the name for the class
  class(result) <- append(class(result), "multiResultClass2")
  return(result)
}


m.sequence <- seq(50, 400, 50)

Runtime.Mean <- matrix(rep(0, length(m.sequence) * 3), nrow = 3)
colnames(Runtime.Mean) <- paste0("n.prov = ", m.sequence)
rownames(Runtime.Mean) <- c("pplasso", "grpreg", "glmnet")

Runtime.sd <- matrix(rep(0, length(m.sequence) * 3), nrow = 3)
colnames(Runtime.sd) <- paste0("n.prov = ", m.sequence)
rownames(Runtime.sd) <- c("pplasso", "grpreg", "glmnet")

iter.Mean <- matrix(rep(0, length(m.sequence) * 2), nrow = 2)
colnames(iter.Mean) <- paste0("n.prov = ", m.sequence)
rownames(iter.Mean) <- c("pplasso with MM", "pplasso without MM")

iter.sd <- matrix(rep(0, length(m.sequence) * 2), nrow = 2)
colnames(iter.sd) <- paste0("n.prov = ", m.sequence)
rownames(iter.sd) <- c("pplasso with MM", "pplasso without MM")

data.loop <- 1:10
ind <- 0

for (j in m.sequence){
  start.time <- Sys.time()
  print(paste0("Iteration for n.prov = ", j, " starts.."))
  ind <- ind + 1
  cl <- makeCluster(5)
  registerDoParallel(cl) 
  Model.Comparison <- 
    foreach (i = data.loop, .packages = c("grpreg", "glmnet", "fastDummies", "RcppArmadillo", "MASS", "Matrix", "TmpLasso")) %dopar% {
      Y.char <- 'Y'
      prov.char <- 'Prov.ID'
      sim.parameters.Lasso <- list(m = j, n.beta = 50, n.groups = 50, prop.NonZero.group = 0.2,
                                   prop.outlier = 0, rho = 0.7)
      Sim_Lasso <- Simulation_data_GroupLasso(sim.parameters.Lasso, prov.size.mean = 200, unpenalized.beta = F)
      data_Lasso <- Sim_Lasso$sim.data
      Z.char <- colnames(data_Lasso)[3:(ncol(data_Lasso) - 1)]
      
      # results from Lasso with MM
      start <- Sys.time()
      cv.model_grp_lasso1 <- cv.pp.lasso(data_Lasso, Y.char, Z.char, prov.char,
                                         MM = T, nfolds = 10, trace.cv = F)
      end <- Sys.time()
      cv.process.time1 <- difftime(end, start, units = 'mins')  #runtime
      total.iter.pplasso1 <- sum((cv.model_grp_lasso1$fit)$iter)
      
      # results from Lasso without MM
      cv.model_grp_lasso2 <- cv.pp.lasso(data_Lasso, Y.char, Z.char, prov.char,
                                         MM = F, nfolds = 10, trace.cv = F)
      total.iter.pplasso2 <- sum((cv.model_grp_lasso2$fit)$iter)
      
      ## results from grpreg
      start <- Sys.time()
      # dummy data for grpreg
      dummy_data <- dummy_cols(data_Lasso, select_columns = prov.char, remove_selected_columns = TRUE, 
                               remove_first_dummy = TRUE)
      ID.char <- colnames(dummy_data)[(ncol(dummy_data) - sim.parameters.Lasso$m + 2):ncol(dummy_data)]
      cv.model_grpreg <- cv.grpreg(dummy_data[,c(Z.char, ID.char)], dummy_data[,Y.char], family = "binomial",
                                   penalty = "grLasso", group = c(1:length(Z.char), rep(0, length(ID.char))),
                                   alpha = 1, nfolds = 10, trace = F)
      end <- Sys.time()
      cv.process.time2 <-  difftime(end, start, units = 'mins') #runtime
      
      ## results from glmnet
      start <- Sys.time()
      # dummy data for glmnet
      dummy_data <- dummy_cols(data_Lasso, select_columns = prov.char, remove_selected_columns = TRUE, 
                               remove_first_dummy = TRUE)
      ID.char <- colnames(dummy_data)[(ncol(dummy_data) - sim.parameters.Lasso$m + 2):ncol(dummy_data)]
      cv.model_glmnet <- cv.glmnet(as.matrix(dummy_data[ ,c(Z.char, ID.char)]), dummy_data[,Y.char], family = "binomial", 
                                   penalty.factor = c(rep(1, length(Z.char)), rep(0, length(ID.char))), 
                                   alpha = 1, nfolds = 10)
      end <- Sys.time()
      cv.process.time3 <-  difftime(end, start, units = 'mins') #runtime
      
      
      result <- multiResultClass2()
      result$Runtime <- round(matrix(c(cv.process.time1, cv.process.time2, cv.process.time3), nrow = 3), digits = 3)
      result$total_iter <- round(matrix(c(total.iter.pplasso1, total.iter.pplasso2), nrow = 2), digits = 4)
      return(result)
    }
  stopCluster(cl)
  
  n.data.loop <- length(data.loop)
  
  Runtime <- matrix(rep(0, 3 * n.data.loop), nrow = 3)
  rownames(Runtime) <- c("pplasso", "grpreg", "glmnet")
  colnames(Runtime) <- paste0("Data_", data.loop)
  
  iter <- matrix(rep(0, 2 * n.data.loop), nrow = 2)
  rownames(iter) <- c("pplasso with MM", "pplasso without MM")
  colnames(iter) <- paste0("Data_", data.loop)
  
  for (i in data.loop){
    Runtime[, i] <- Model.Comparison[[i]]$Runtime
    iter[, i] <- Model.Comparison[[i]]$total_iter
  }
  
  Runtime.Mean[, ind] <- round(apply(Runtime, 1, mean), digits = 3)
  Runtime.sd[, ind] <- round(apply(Runtime, 1, sd), digits = 3)
  
  iter.Mean[, ind] <- round(apply(iter, 1, mean), digits = 3)
  iter.sd[, ind] <- round(apply(iter, 1, sd), digits = 3)
  
  end.time <- Sys.time()
  process.time <- difftime(end.time, start.time, units = 'mins')
  print(paste0("Processing time: ", round(process.time, digits = 4), " mins (current n.prov = ", j, ")"))
  
}

### Figures (Runtime, iter)
library(reshape2)

#### Figure1: Runtime
n.prov <- m.sequence
Runtime.Mean.lower <- Runtime.Mean - Runtime.sd
Runtime.Mean.upper <- Runtime.Mean + Runtime.sd


Runtime.figure.pplasso.df <- data.frame("model" = melt(Runtime.Mean)$Var1,
                                        "n.prov" = rep(n.prov, each = 3),
                                        "Mean" = melt(Runtime.Mean)$value,
                                        "lower" = melt(Runtime.Mean.lower)$value,
                                        "upper" = melt(Runtime.Mean.upper)$value)

Runtime.plot.pplasso <- ggplot(Runtime.figure.pplasso.df, aes(n.prov, group = factor(model))) +
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
  scale_linetype_manual(values = c("solid", "dotdash", "dotted"), name = "", labels = c("pplasso", "grpreg", "glmnet")) + 
  scale_color_manual(values = c("blue", "red", "black"), name = "", labels = c("pplasso", "grpreg", "glmnet")) + 
  scale_fill_manual(values = c("grey80", "grey85", "grey90"), name = "", labels = c("pplasso", "grpreg", "glmnet"))

Runtime.plot.pplasso <- annotate_figure(ggarrange(Runtime.plot.pplasso,
                                                  nrow = 1, ncol = 1,
                                                  labels = c("A")))

save(Runtime.figure.pplasso.df, Runtime.plot.pplasso,
     file = paste0("Runtime_pplasso_", Sys.Date(), ".RData"))

#### Figure2: iterations
n.prov <- m.sequence
iter.Mean.lower <- iter.Mean - iter.sd
iter.Mean.upper <- iter.Mean + iter.sd

iter.figure.pplasso.df <- data.frame("model" = melt(iter.Mean)$Var1,
                                     "n.prov" = rep(n.prov, each = 2),
                                     "Mean" = melt(iter.Mean)$value,
                                     "lower" = melt(iter.Mean.lower)$value,
                                     "upper" = melt(iter.Mean.upper)$value)

iter.plot.pplasso <- ggplot(iter.figure.pplasso.df, aes(n.prov, group = factor(model))) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(model)), alpha = 0.5) +
  geom_line(aes(y = Mean, color = factor(model), linetype = factor(model)), size = 1)  + 
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(size = 13, face = "bold", family = "serif"),
        axis.title = element_text(size = 13, family = "serif"),
        plot.caption = element_text(size = 10, face = "italic", family = "serif"),
        legend.text = element_text(size = 12, family = "serif")) + 
  theme(axis.text = element_text(face = "italic", size = 12, family = "serif")) + 
  theme(legend.position = c(0.22, 0.85),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank()) + 
  theme(legend.key.height = unit(0.6, 'cm'),
        legend.key.width = unit(1.2, 'cm'),
        legend.spacing.x = unit(0, 'cm')) + 
  labs(title = "", 
       x = "Number of providers", 
       y = "Total #iterations",
       caption = "( provider size varies from 50 to 400 )") +
  scale_x_continuous(breaks = seq(50, 400, 50)) + 
  scale_linetype_manual(values = c("dotdash", "solid"), name = "", labels = c("  pplasso with MM", "  pplasso without MM")) + 
  scale_color_manual(values = c("red", "blue"), name = "", labels = c("  pplasso with MM", "  pplasso without MM")) + 
  scale_fill_manual(values = c("grey85", "grey90"), name = "", labels = c("  pplasso with MM", "  pplasso without MM")) +
  ylim(0, 7000)


save(iter.figure.pplasso.df, iter.plot.pplasso,
     file = paste0("iter_num_pplasso_MM_", Sys.Date(), ".RData"))
