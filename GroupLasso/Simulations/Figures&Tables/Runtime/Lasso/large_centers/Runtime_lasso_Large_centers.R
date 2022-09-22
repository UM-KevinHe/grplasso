setwd("/home/ybshao/GroupLasso/Simulations/Figures&Tables/Runtime/Lasso/large_centers")
library(TmpLasso)
library(MASS)
library(RcppArmadillo)
library(Rcpp)
library(Matrix)
library(grpreg)
library(fastDummies)
library(foreach)
library(doParallel)
library(glmnet)

Y.char <- 'Y'
prov.char <- 'Prov.ID'
set.seed(1)
sim.parameters.Lasso <- list(m = 1000, n.beta = 50, n.groups = 50, prop.NonZero.group = 0.2,
                             prop.outlier = 0, rho = 0.7)
Sim_Lasso <- Simulation_data_GroupLasso(sim.parameters.Lasso, prov.size.mean = 200, unpenalized.beta = F)
data_Lasso <- Sim_Lasso$sim.data
Z.char <- colnames(data_Lasso)[3:(ncol(data_Lasso) - 1)]

# results from Lasso with MM
start <- Sys.time()
cv.model_grp_lasso <- cv.pp.lasso(data_Lasso, Y.char, Z.char, prov.char,
                                  MM = T, nfolds = 10, trace.cv = F)
end <- Sys.time()
cv.process.time.grplasso <- difftime(end, start, units = 'mins')  
save(cv.process.time.grplasso,
     file = paste0("large_centers_runtime_grplasso_m1000_", Sys.Date(), ".RData"))

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
cv.process.time.glmnet <-  difftime(end, start, units = 'mins') #runtime

save(cv.process.time.glmnet,
     file = paste0("large_centers_runtime_glmnet_m1000_", Sys.Date(), ".RData"))

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
cv.process.time.grpreg <-  difftime(end, start, units = 'mins') #runtime


save(cv.process.time.grpreg,
     file = paste0("large_centers_runtime_grpreg_m1000_", Sys.Date(), ".RData"))