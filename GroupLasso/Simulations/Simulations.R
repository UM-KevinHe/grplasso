library(MASS)
library(RcppArmadillo)
library(Rcpp)
library(Matrix)
library(grpreg)
library(fastDummies)
library(foreach)
library(doParallel)
library(glmnet)


######-------Import relevant functions-------------######
source("GroupLasso/Relevant Functions/Simulation_Functions.R")
source("GroupLasso/Relevant Functions/Relevant_Functions.R")
source("GroupLasso/grLasso_Code.R")

Rcpp::sourceCpp("GroupLasso/Relevant Functions/Relevant_Functions.cpp")
Rcpp::sourceCpp("GroupLasso/grLasso_Code.cpp")

######-----------------------------------------#####
######----------- 1. Estimation Comparison: multi-centers ------------######
set.seed(1)
Y.char <- 'Y'
prov.char <- 'Prov.ID'
sim.parameters.GrLasso <- list(m = 5, n.beta = 10, n.groups = 3, prop.NonZero.group = 0.2, 
                               prop.outlier = 0, rho = 0.7)
Sim_GrLasso <- Simulation_data_GroupLasso(sim.parameters.GrLasso, unpenalized.beta = F)
data_GrLasso <- Sim_GrLasso$sim.data
group <- Sim_GrLasso$group
true.beta <- Sim_GrLasso$beta
true.gamma <- Sim_GrLasso$gamma
Z.char <- colnames(data_GrLasso)[3:(ncol(data_GrLasso) - 1)]
#data_GrLasso <- fe.data.prep(data_GrLasso, Y.char, Z.char, prov.char, cutoff = 0, check = FALSE)
#data_GrLasso <- data_GrLasso[data_GrLasso$included == 1, ]
true.mu <- data_GrLasso$mu


# dummy data for grpreg
dummy_data <- dummy_cols(data_GrLasso, select_columns = prov.char, remove_selected_columns = TRUE, 
                         remove_first_dummy = TRUE)
ID.char <- colnames(dummy_data)[(ncol(dummy_data) - sim.parameters.GrLasso$m + 2):ncol(dummy_data)]


## 1. gr_ppLasso
start <- Sys.time()
model_grp_lasso <- grp.lasso(data_GrLasso, Y.char, Z.char, prov.char, group = group, 
                             trace.lambda = T, backtrack = TRUE)
end <- Sys.time()
process.time1 <- end - start
print(process.time1)

grLasso.beta <- model_grp_lasso$beta
grLasso.gamma <- model_grp_lasso$gamma
grLasso.lambda <- model_grp_lasso$lambda
grLasso.loss <- model_grp_lasso$loss

## 2. grpreg
start <- Sys.time()
model_grpreg <- grpreg(dummy_data[,c(Z.char, ID.char)], dummy_data[,Y.char], family = "binomial",
                       penalty = "grLasso", group = c(group, rep(0, length(ID.char))), alpha = 1)
end <- Sys.time()
process.time2 <- end - start
print(process.time2)

grpreg.beta <- model_grpreg$beta[2:(1 + sim.parameters.GrLasso$n.beta),]

if (ncol(model_grpreg$beta) == 1){
  grpreg.beta <- matrix(grpreg.beta, ncol = 1)
  grpreg.gamma <- c(model_grpreg$beta[1,],
                        model_grpreg$beta[1,] +
                          model_grpreg$beta[(sim.parameters.GrLasso$n.beta + 2):length(model_grpreg$beta[, 1]),])
  names(grpreg.gamma) <-  names(table(data_GrLasso$Prov.ID))
  grpreg.gamma <- matrix(grpreg.gamma, ncol = 1)
  n <- 1
} else {
  grpreg.gamma <- rbind(model_grpreg$beta[1,],
                        model_grpreg$beta[1,] +
                          model_grpreg$beta[(sim.parameters.GrLasso$n.beta + 2):length(model_grpreg$beta[, 1]),])
  rownames(grpreg.gamma) <- names(table(data_GrLasso$Prov.ID))
}

grpreg.lambda <- model_grpreg$lambda
grpreg.loss <- model_grpreg$loss

# estimation difference
n <- min(ncol(grLasso.beta), ncol(grpreg.beta)) - 1
beta.diff <- round(grLasso.beta[, 1:n] - grpreg.beta[, 1:n], digits = 4)
gamma.diff <- round(grLasso.gamma[, 1:n] - grpreg.gamma[, 1:n], digits = 4)


######----------- 2. Estimation Comparison: With only one provider (only intercept) ------------######
sim.parameters2 <- list(m = 1, n.beta = 20, n.groups = 5, prop.NonZero.group = 1, 
                        prop.outlier = 0, rho = 0.6)
set.seed(1)
## Note: If any beta are unpenalized, grpreg will add 1e-05 to the maximum lambda. 
#Sim_2 <- Simulation_data_GroupLasso(sim.parameters2, unpenalized.beta = T, prop.unpenalized.beta = 0.5)
Sim_2 <- Simulation_data_GroupLasso(sim.parameters2, unpenalized.beta = F)
data2 <- Sim_2$sim.data
Z.char <- colnames(data2)[3:(ncol(data2) - 1)]
Y.char <- 'Y'
prov.char <- 'Prov.ID'
#data2 <- fe.data.prep(data2, Y.char, Z.char, prov.char, cutoff = 0, check = FALSE)
original.group2 <- Sim_2$group

## 1. gr_ppLasso
start <- Sys.time()
model_grp_lasso_onlyIntercept <- grp.lasso(data2, Y.char, Z.char, trace.lambda = F, group = original.group2, 
                                           backtrack = F, tol = 1e-10, max.iter = 1e6, returnX = F)
end <- Sys.time()
process.time1 <- end - start
print(process.time1)

grLasso.beta.onlyIntercept <- model_grp_lasso_onlyIntercept$beta
grLasso.gamma.onlyIntercept <- model_grp_lasso_onlyIntercept$gamma
grLasso.lambda.onlyIntercept <- model_grp_lasso_onlyIntercept$lambda
grLasso.loss.onlyIntercept <- model_grp_lasso_onlyIntercept$loss


## 2. grpreg
start <- Sys.time()
model_grpreg_onlyIntercept <- grpreg(data2[,c(Z.char)], data2[,Y.char], family = "binomial", penalty = "grLasso",
                                     group = original.group2, alpha = 1, eps = 1e-10, max.iter = 1e6, returnX = F)
end <- Sys.time()
process.time2 <- end - start
print(process.time2)

grpreg.beta.onlyIntercept <- model_grpreg_onlyIntercept$beta[2:(1 + sim.parameters2$n.beta),]
grpreg.gamma.onlyIntercept <- matrix(model_grpreg_onlyIntercept$beta[1, ], nrow = 1)
rownames(grpreg.gamma.onlyIntercept) <- "Prov_1"
colnames(grpreg.gamma.onlyIntercept) <- round(model_grpreg_onlyIntercept$lambda, digits = 4)
grpreg.lambda <- model_grpreg_onlyIntercept$lambda
grpreg.loss <- model_grpreg_onlyIntercept$loss


#comparison

beta.diff <- round(grLasso.beta.onlyIntercept - grpreg.beta.onlyIntercept, digits = 5)
print(max(abs(beta.diff)))
gamma.diff <- round(grLasso.gamma.onlyIntercept - grpreg.gamma.onlyIntercept, digits = 5)
print(max(abs(gamma.diff)))



######----------- 3. Character provider names ------------######
set.seed(1)
Y.char <- 'Y'
prov.char <- 'Prov.ID'
sim.parameters.GrLasso <- list(m = 5, n.beta = 5, n.groups = 5, prop.NonZero.group = 1, 
                               prop.outlier = 0, rho = 0.7)
Sim_GrLasso <- Simulation_data_GroupLasso(sim.parameters.GrLasso, unpenalized.beta = F)
data_GrLasso <- Sim_GrLasso$sim.data
data_GrLasso[which(data_GrLasso$Prov.ID == 1), 2] <- "No.A"
data_GrLasso[which(data_GrLasso$Prov.ID == 2), 2] <- "B"
data_GrLasso[which(data_GrLasso$Prov.ID == 3), 2] <- "Facility.C"
data_GrLasso[which(data_GrLasso$Prov.ID == 4), 2] <- "Fourth"
data_GrLasso[which(data_GrLasso$Prov.ID == 5), 2] <- "No.5"
group <- Sim_GrLasso$group
true.beta <- Sim_GrLasso$beta
true.gamma <- Sim_GrLasso$gamma
Z.char <- colnames(data_GrLasso)[3:(ncol(data_GrLasso) - 1)]
#data_GrLasso <- fe.data.prep(data_GrLasso, Y.char, Z.char, prov.char, cutoff = 0, check = FALSE)
#data_GrLasso <- data_GrLasso[data_GrLasso$included == 1, ]

model_grp_lasso_char <- grp.lasso(data_GrLasso, Y.char, Z.char, prov.char, group = group, trace.lambda = F)

grLasso.beta <- model_grp_lasso_char$beta
grLasso.gamma <- model_grp_lasso_char$gamma

## compare with grpreg
n.show.lambda <- 1:10
dummy_data <- dummy_cols(data_GrLasso, select_columns = prov.char, remove_selected_columns = TRUE, 
                         remove_first_dummy = TRUE)
ID.char <- colnames(dummy_data)[(ncol(dummy_data) - sim.parameters.GrLasso$m + 2):ncol(dummy_data)]
model_grpreg_char <- grpreg(dummy_data[, c(Z.char, ID.char)], dummy_data[,Y.char], family = "binomial",
                            penalty = "grLasso", group = c(group, rep(0, length(ID.char))), alpha = 1)
grpreg.beta <- model_grpreg_char$beta[2:(1 + sim.parameters.GrLasso$n.beta),]
grpreg.gamma <- rbind(model_grpreg_char$beta[1,],
                      model_grpreg_char$beta[1,] +
                        model_grpreg_char$beta[(sim.parameters.GrLasso$n.beta + 2):length(model_grpreg_char$beta[, 1]),])


######----------- Cross Validation ------------######
set.seed(1)
Y.char <- 'Y'
prov.char <- 'Prov.ID'
sim.parameters.GrLasso <- list(m = 5, n.beta = 10, n.groups = 3, prop.NonZero.group = 0.2, 
                               prop.outlier = 0, rho = 0.7)
Sim_GrLasso <- Simulation_data_GroupLasso(sim.parameters.GrLasso, unpenalized.beta = F)
data_GrLasso <- Sim_GrLasso$sim.data
group <- Sim_GrLasso$group
true.beta <- Sim_GrLasso$beta
true.gamma <- Sim_GrLasso$gamma
Z.char <- colnames(data_GrLasso)[3:(ncol(data_GrLasso) - 1)]
true.mu <- data_GrLasso$mu


# dummy data for grpreg
dummy_data <- dummy_cols(data_GrLasso, select_columns = prov.char, remove_selected_columns = TRUE, 
                         remove_first_dummy = TRUE)
ID.char <- colnames(dummy_data)[(ncol(dummy_data) - sim.parameters.GrLasso$m + 2):ncol(dummy_data)]


## 1.gr_ppLasso
start <- Sys.time()
cv.model_grp_lasso <- cv.grp.lasso(data_GrLasso, Y.char, Z.char, prov.char, group = group, 
                                   trace.lambda = F, nfolds = 10, trace.cv = T)
end <- Sys.time()
cv.process.time1 <- end - start

cv_BestModel_grp_lasso <- cv.model_grp_lasso$fit
best.beta.grp_lasso <- cv_BestModel_grp_lasso$beta[, cv.model_grp_lasso$min]
best.lambda.grp_lasso <- cv.model_grp_lasso$lambda.min
best.eta.grp_lasso <- cv_BestModel_grp_lasso$linear.predictors[, cv.model_grp_lasso$min]

RME.grp_lasso <- sqrt(sum((best.beta.grp_lasso - true.beta)^2)/sim.parameters.GrLasso$n.beta)
RMSE.grp_lasso <- sqrt(sum((plogis(best.eta.grp_lasso) - true.mu)^2)/nrow(data_GrLasso))
cross_entropy.grp_lasso <- cv.model_grp_lasso$cve[cv.model_grp_lasso$min]
wrong.prediction.rate.grp_lasso <- cv.model_grp_lasso$pe[cv.model_grp_lasso$min]


## 2.grpreg
start <- Sys.time()
cv.model_grpreg <- cv.grpreg(dummy_data[,c(Z.char, ID.char)], dummy_data[,Y.char], family = "binomial",
                             penalty = "grLasso", group = c(group, rep(0, length(ID.char))), alpha = 1, 
                             nfolds = 10, trace.cv = T)
end <- Sys.time()
cv.process.time2 <- end - start

cv_BestModel_grpreg <- cv.model_grpreg$fit
best.beta.grpreg <- cv_BestModel_grpreg$beta[2:(1 + sim.parameters.GrLasso$n.beta), cv.model_grpreg$min]
best.lambda.grpreg <- cv.model_grpreg$lambda.min
best.eta.grpreg <- cv_BestModel_grpreg$linear.predictors[, cv.model_grpreg$min]


RME.grpreg <- sqrt(sum((best.beta.grpreg - true.beta)^2)/sim.parameters.GrLasso$n.beta)
RMSE.grpreg <- sqrt(sum((plogis(best.eta.grpreg) - true.mu)^2)/nrow(data_GrLasso))
cross_entropy.grpreg <- cv.model_grpreg$cve[cv.model_grpreg$min]
wrong.prediction.rate.grpreg <- cv.model_grpreg$pe[cv.model_grpreg$min]




######----------- Figure: Entire Regularization Path ------------######

# 1. Simulated at: 100 providers, totally 10 beta's with 3 groups (1 zero group, 2 penalized non-zero groups)
library(ggplot2)
Y.char <- 'Y'
prov.char <- 'Prov.ID'
set.seed(99)
sim.parameters.GrLasso <- list(m = 100, n.beta = 10, n.groups = 3, prop.NonZero.group = 0.6, 
                               prop.outlier = 0, rho = 0.5)
Sim_GrLasso <- Simulation_data_GroupLasso(sim.parameters.GrLasso, unpenalized.beta = F)
data_GrLasso <- Sim_GrLasso$sim.data
group <- Sim_GrLasso$group
true.beta <- Sim_GrLasso$beta
Z.char <- colnames(data_GrLasso)[3:(ncol(data_GrLasso) - 1)]

RegPath.grp_lasso <- grp.lasso(data_GrLasso, Y.char, Z.char, prov.char, group = group, trace.lambda = T)
grLasso.beta <- RegPath.grp_lasso$beta

iter.num <- rep(log(RegPath.grp_lasso$lambda), each = nrow(grLasso.beta))
y <- as.vector(grLasso.beta)
group <- rep(1:nrow(grLasso.beta), length(RegPath.grp_lasso$lambda))
path.figure.df <- as.data.frame(t(rbind(iter.num, y, group)))

Regularization.path1 <- ggplot(path.figure.df, aes(iter.num, y, group = factor(group))) + 
  geom_line(aes(color = factor(group)), size = 0.5) + 
  geom_abline(color = "red", linetype = 2, size = 0.5, intercept = 0, slope = 0) + 
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black")) + 
  theme(legend.text = element_text(colour = "black", size = 8, face = "italic"), 
        legend.text.align = 0, legend.title = 
          element_text(colour = "black", size = 10),
        legend.title.align = 0.5) +
  theme(plot.title = element_text(size = 13, face="bold", family = "serif"),
        axis.title = element_text(size = 12, family = "serif"),
        plot.caption = element_text(size = 10, face = "italic", family = "serif"),
        legend.title = element_text(size = 10, family = "serif"),
        legend.text = element_text(size = 8, family = "serif")) + 
  theme(axis.text = element_text(face = "italic", size = 11)) + 
  labs(title = "Regularization Coefficient Paths", 
       x = expression(log(lambda)), y = expression(paste(beta, " coefficients")),
       caption = "( Note: 2 of the 3 groups have true non-zero effect )") +
  scale_color_manual(values = 1:nrow(grLasso.beta), name = expression(paste(beta)), 
                     labels = c(1:nrow(grLasso.beta))) + 
  scale_x_continuous(trans = scales::reverse_trans(), breaks = round(seq(round(max(iter.num), 0), round(min(iter.num), 0), by = - 1), 1))
Regularization.path1


# 2. Simulated at: 100 providers, totally 10 beta's with 4 groups (1 unpenalized group, 1 zero group, 2 penalized non-zero groups)
Y.char <- 'Y'
prov.char <- 'Prov.ID'
set.seed(10)
sim.parameters.GrLasso <- list(m = 100, n.beta = 10, n.groups = 4, prop.NonZero.group = 0.5, 
                               prop.outlier = 0, rho = 0.5)
Sim_GrLasso <- Simulation_data_GroupLasso(sim.parameters.GrLasso, unpenalized.beta = T)
data_GrLasso <- Sim_GrLasso$sim.data
group <- Sim_GrLasso$group
true.beta <- Sim_GrLasso$beta
Z.char <- colnames(data_GrLasso)[3:(ncol(data_GrLasso) - 1)]


RegPath.grp_lasso <- grp.lasso(data_GrLasso, Y.char, Z.char, prov.char, group = group, trace.lambda = T)
grLasso.beta <- RegPath.grp_lasso$beta

iter.num <- rep(log(RegPath.grp_lasso$lambda), each = nrow(grLasso.beta))
y <- as.vector(grLasso.beta)
group <- rep(1:nrow(grLasso.beta), length(RegPath.grp_lasso$lambda))
path.figure.df <- as.data.frame(t(rbind(iter.num, y, group)))

Regularization.path2 <- ggplot(path.figure.df, aes(iter.num, y, group = factor(group))) + 
  geom_line(aes(color = factor(group)), size = 0.5) + 
  geom_abline(color = "red", linetype = 2, size = 0.5, intercept = 0, slope = 0) + 
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black")) + 
  theme(legend.text = element_text(colour = "black", size = 8, face = "italic"), 
        legend.text.align = 0, legend.title = 
          element_text(colour = "black", size = 10),
        legend.title.align = 0.5) +
  theme(plot.title = element_text(size = 13, face="bold", family = "serif"),
        axis.title = element_text(size = 12, family = "serif"),
        plot.caption = element_text(size = 10, face = "italic", family = "serif"),
        legend.title = element_text(size = 10, family = "serif"),
        legend.text = element_text(size = 8, family = "serif")) + 
  theme(axis.text = element_text(face = "italic", size = 11)) + 
  labs(title = "Regularization Coefficient Paths", x = expression(log(lambda)), 
       y = expression(paste(beta, " coefficients")),
       caption = "( Note: 3 of the 4 groups have true non-zero effect, and 1 of the 3 non-zero groups is unpenalized)") +
  scale_color_manual(values = 1:nrow(grLasso.beta), name = expression(paste(beta)), 
                     labels = c(1:nrow(grLasso.beta))) + 
  scale_x_continuous(trans = scales::reverse_trans(), breaks = round(seq(round(max(iter.num), 0), round(min(iter.num), 0), by = - 1), 1))
Regularization.path2

######----------- Figure: Cross-validation error at entire lambda sequence (Based on Real Data)------------######
data(Birthwt)
Z <- Birthwt$X
Y <- matrix(Birthwt$low, nrow = nrow(Z))
Y.char <- 'Y'
data <- as.data.frame(cbind(Y, Z))
Z.char <- colnames(Z)
colnames(data) <- c(Y.char, Z.char)
group <- Birthwt$group
cv.model_grp_lasso <- cv.grp.lasso(data, Y.char, Z.char, group = group, 
                                   trace.lambda = F, nfolds = 10, trace.cv = T)

CVE <- cv.model_grp_lasso$cve
CVE.upper <- cv.model_grp_lasso$cve + cv.model_grp_lasso$cvse
CVE.lower <- cv.model_grp_lasso$cve - cv.model_grp_lasso$cvse

log.lambda <- log(cv.model_grp_lasso$lambda)
CV.figure.df <- as.data.frame(cbind(log.lambda, CVE, CVE.upper, CVE.lower))

cv.plot.grLasso <- ggplot(CV.figure.df, aes(log.lambda, CVE)) +
  geom_line(aes(y = CVE), size = 0.05, color = "blue") + 
  geom_vline(xintercept = log(cv.model_grp_lasso$lambda.min), size = 0.5, linetype = "dashed", color = "blue") + 
  geom_point(size = 1, color = "red") +
  geom_errorbar(aes(ymin = CVE.lower, ymax = CVE.upper), width = 0.1, size = 0.1) +
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(size = 13, face="bold", family = "serif"),
        axis.title = element_text(size = 12, family = "serif"),
        plot.caption = element_text(size = 10, face = "italic", family = "serif"),
        legend.title = element_text(size = 12, family = "serif"),
        legend.text = element_text(size = 11, family = "serif"))  + 
  theme(axis.text = element_text(face = "italic", size = 11, family = "serif")) + 
  labs(title = "Cross Entropy Loss at Various Regularization Coefficients", 
       x = expression(log(lambda)), 
       y = "cross validation error",
       caption = "( Based on Birthwt dataset provided by grpreg package)") +
  scale_x_continuous(trans = scales::reverse_trans(), 
                     breaks = round(seq(round(max(log.lambda), 0), round(min(log.lambda), 0), by = - 1), 1))
cv.plot.grLasso


X.grpreg <- Birthwt$X
y.grpreg <- Birthwt$low
cv.model_grpreg <- cv.grpreg(X.grpreg, y.grpreg, family = "binomial", penalty = "grLasso", 
                             group = group, alpha = 1, nfolds = 10, trace.cv = T)
cv.model_grpreg$lambda.min
plot(cv.model_grpreg)




######----------- Real Data Estimation: (Based on "Birthwt") ------------######
data(Birthwt)
Z <- Birthwt$X
Y <- matrix(Birthwt$low, nrow = nrow(Z))
data <- as.data.frame(cbind(Y, Z))
Z.char <- colnames(Z)
Y.char <- 'Y'
colnames(data) <- c(Y.char, Z.char)
group <- Birthwt$group
RD.model_grp_lasso <- grp.lasso(data, Y.char, Z.char, group = group, trace.lambda = T, backtrack = F)
grLasso.beta <- RD.model_grp_lasso$beta
grLasso.gamma <- RD.model_grp_lasso$gamma
rownames(grLasso.gamma) <- "Intercept"
grLasso.lambda <- RD.model_grp_lasso$lambda
grLasso.loss <- RD.model_grp_lasso$loss


X.grpreg <- Birthwt$X
y.grpreg <- Birthwt$low
RD.model_grpreg <- grpreg(X.grpreg, y.grpreg, family = "binomial", penalty = "grLasso", 
                             group = group, alpha = 1)
grpreg.beta <- RD.model_grpreg$beta[2:(1 + length(Z.char)),]
grpreg.gamma <- matrix(RD.model_grpreg$beta[1, ], nrow = 1)
rownames(grpreg.gamma) <- "Intercept"
colnames(grpreg.gamma) <- round(RD.model_grpreg$lambda, digits = 4)
grpreg.lambda <- RD.model_grpreg$lambda
grpreg.loss <- RD.model_grpreg$loss

# estimation difference
n <- min(ncol(grLasso.beta), ncol(grpreg.beta)) 
beta.diff <- round(grLasso.beta[, 1:n] - grpreg.beta[, 1:n], digits = 4)
print(max(abs(beta.diff)))
gamma.diff <- round(grLasso.gamma[, 1:n] - grpreg.gamma[, 1:n], digits = 4)
print(max(abs(gamma.diff)))




######----------- Time Comparison ------------######
######----------- Time Comparison ------------######
######----------- Time Comparison ------------######

# Since the "foreach" package has problems on directly executing Rcpp functions,
# we need to create a simple package for our Rcpp functions! 

## After run "Rcpp.package.skeleton" line, we need to manually add ", RcppArmadillo" in "DESCRIPTION" file, 
## and replace "Rcpp.h" by "RcppArmadillo.h" in "src/RcppExports.cpp" file 
#Rcpp.package.skeleton("TmpGrlasso", path = "GroupLasso/Simulations/MyCpp",
#                      code_files = c("GroupLasso/Relevant Functions/Relevant_Functions.R",
#                                     "GroupLasso/grLasso_Code.R"),
#                      cpp_files = c("GroupLasso/grLasso_Code.cpp", 
#                                    "GroupLasso/Relevant Functions/Relevant_Functions.cpp"))

#RcppArmadillo.package.skeleton("compare", path = "GroupLasso/Simulations/MyCpp")
#install.packages("GroupLasso/Simulations/MyCpp/TmpGrlasso", repos = NULL, type="source") 
#library(TmpGrlasso)

# 1. n.beta = 100, n.groups = 20, non.zero.groups = 4, n.provider range from 200 to 3000

# Under each provider counts, 10 data are simulated; 
# "mean \pm sd runtime" based on 5-fold cross validation on entire lambda path will be reported; 
# RME and RMSE will be computed based on the "best lambda" selected by 10-fod cross validation;
# each data has 100 betas which have been randomly assigned into 20 groups. 4 groups contains non-zero beta, while other beta's are set to zero;
# no unpenalized beta's have been specified;

multiResultClass <- function(Runtime = NULL, RME = NULL, RMSE = NULL, cross_entropy = NULL, wrong_prediction_rate = NULL){
  result <- list(
    Runtime = Runtime,
    RME = RME,
    RMSE = RMSE,
    cross_entropy = cross_entropy,
    wrong_prediction_rate = wrong_prediction_rate
  )
  ## Set the name for the class
  class(result) <- append(class(result), "multiResultClass")
  return(result)
}


m.sequence <- seq(50, 2000, 50)

Runtime.Mean <- matrix(rep(0, length(m.sequence) * 2), nrow = 2)
colnames(Runtime.Mean) <- paste0("n.prov = ", m.sequence)
rownames(Runtime.Mean) <- c("grLasso", "grpreg")

Runtime.sd <- matrix(rep(0, length(m.sequence) * 2), nrow = 2)
colnames(Runtime.sd) <- paste0("n.prov = ", m.sequence)
rownames(Runtime.sd) <- c("grLasso", "grpreg")

RME.Mean <- matrix(rep(0, length(m.sequence) * 2), nrow = 2)
colnames(RME.Mean) <- paste0("n.prov = ", m.sequence)
rownames(RME.Mean) <- c("grLasso", "grpreg")

RME.sd <- matrix(rep(0, length(m.sequence) * 2), nrow = 2)
colnames(RME.sd) <- paste0("n.prov = ", m.sequence)
rownames(RME.sd) <- c("grLasso", "grpreg")

RMSE.Mean <- matrix(rep(0, length(m.sequence) * 2), nrow = 2)
colnames(RMSE.Mean) <- paste0("n.prov = ", m.sequence)
rownames(RMSE.Mean) <- c("grLasso", "grpreg")

RMSE.sd <- matrix(rep(0, length(m.sequence) * 2), nrow = 2)
colnames(RMSE.sd) <- paste0("n.prov = ", m.sequence)
rownames(RMSE.sd) <- c("grLasso", "grpreg")

cross_entropy.Mean <- matrix(rep(0, length(m.sequence) * 2), nrow = 2)
colnames(cross_entropy.Mean) <- paste0("n.prov = ", m.sequence)
rownames(cross_entropy.Mean) <- c("grLasso", "grpreg")

cross_entropy.sd <- matrix(rep(0, length(m.sequence) * 2), nrow = 2)
colnames(cross_entropy.sd) <- paste0("n.prov = ", m.sequence)
rownames(cross_entropy.sd) <- c("grLasso", "grpreg")

wrong_prediction_rate.Mean <- matrix(rep(0, length(m.sequence) * 2), nrow = 2)
colnames(wrong_prediction_rate.Mean) <- paste0("n.prov = ", m.sequence)
rownames(wrong_prediction_rate.Mean) <- c("grLasso", "grpreg")

wrong_prediction_rate.sd <- matrix(rep(0, length(m.sequence) * 2), nrow = 2)
colnames(wrong_prediction_rate.sd) <- paste0("n.prov = ", m.sequence)
rownames(wrong_prediction_rate.sd) <- c("grLasso", "grpreg")



data.loop <- 1:10
ind <- 0

for (j in m.sequence){ #outer loop for
  ind <- ind + 1
  
  cl.cores <- detectCores()
  cl <- makeCluster(cl.cores - 1)
  registerDoParallel(cl) 
  Model.Comparison <- 
    foreach (i = data.loop, .packages = c("grpreg", "fastDummies", "RcppArmadillo", "MASS", "Matrix", "myRcpp"),
             .noexport = c("Deviance", "grp_lasso", "logis_fe_prov", "SerBIN", "Z_max_grLasso")) %dopar% {
               source("GroupLasso/Relevant Functions/Relevant_Functions.R")
               source("GroupLasso/Relevant Functions/Simulated_Functions.R")
               source("GroupLasso/grLasso_Code.R")
               Y.char <- 'Y'
               prov.char <- 'Prov.ID'
               sim.parameters.GrLasso <- list(m = j, n.beta = 50, n.groups = 10, prop.NonZero.group = 0.2, 
                                              prop.outlier = 0.05, rho = 0.7)
               Sim_GrLasso <- Simulation_data_GroupLasso(sim.parameters.GrLasso, unpenalized.beta = F)
               data_GrLasso <- Sim_GrLasso$sim.data
               group <- Sim_GrLasso$group
               
               true.beta <- Sim_GrLasso$beta  #true beta's for computing RMSE
               Z.char <- paste0('Z_', 1:sim.parameters.GrLasso$n.beta)
               data_prep <- fe.data.prep(data_GrLasso, Y.char, Z.char, prov.char, cutoff = 0, check = FALSE)
               data_prep <- data_prep[data_prep$included == 1, ]
               true.mu <- data_prep$mu #true mu's for computing RME
               
               # dummy data for grpreg
               dummy_data <- dummy_cols(data_prep, select_columns = prov.char, remove_selected_columns = TRUE, 
                                        remove_first_dummy = TRUE)
               ID.char <- rep(NA, sim.parameters.GrLasso$m - 1)
               for (i in 1:(sim.parameters.GrLasso$m - 1)){
                 ID.char[i] <- paste0("Prov.ID_", i + 1)
               }
               
               # results from grLasso
               start <- Sys.time()
               cv.model_grp_lasso <- cv.grp.lasso(data_prep, Y.char, Z.char, prov.char, group = group, trace.lambda = T,
                                                  nfolds = 5, trace.cv = F)
               end <- Sys.time()
               cv.process.time1 <- difftime(end, start, units = 'mins')  #runtime
               cv_BestModel_grp_lasso <- cv.model_grp_lasso$fit
               best.beta.grp_lasso <- cv_BestModel_grp_lasso$beta[, cv.model_grp_lasso$min]
               best.eta.grp_lasso <- cv_BestModel_grp_lasso$linear.predictors[, cv.model_grp_lasso$min]
               
               RME.grp_lasso <- sqrt(sum((best.beta.grp_lasso - true.beta)^2)/sim.parameters.GrLasso$n.beta) #RME
               RMSE.grp_lasso <- sqrt(sum((plogis(best.eta.grp_lasso) - true.mu)^2)/nrow(data_prep)) #RMSE
               cross_entropy.grp_lasso <- cv.model_grp_lasso$cve[cv.model_grp_lasso$min] #CVE
               wrong.prediction.rate.grp_lasso <- cv.model_grp_lasso$pe[cv.model_grp_lasso$min] #PE
               
               
               ## results from grpreg
               start <- Sys.time()
               cv.model_grpreg <- cv.grpreg(dummy_data[,c(Z.char, ID.char)], dummy_data[,Y.char], family = "binomial",
                                            penalty = "grLasso", group = c(group, rep(0, length(ID.char))), alpha = 1, 
                                            nfolds = 5, trace.cv = T)
               end <- Sys.time()
               cv.process.time2 <-  difftime(end, start, units = 'mins') #runtime
               
               cv_BestModel_grpreg <- cv.model_grpreg$fit
               best.beta.grpreg <- cv_BestModel_grpreg$beta[2:(1 + sim.parameters.GrLasso$n.beta), cv.model_grpreg$min]
               best.eta.grpreg <- cv_BestModel_grpreg$linear.predictors[, cv.model_grpreg$min]
               
               RME.grpreg <- sqrt(sum((best.beta.grpreg - true.beta)^2)/sim.parameters.GrLasso$n.beta)  #RME
               RMSE.grpreg <- sqrt(sum((plogis(best.eta.grpreg) - true.mu)^2)/nrow(data_prep))  #RMSE
               cross_entropy.grpreg <- cv.model_grpreg$cve[cv.model_grpreg$min]  #CVE
               wrong.prediction.rate.grpreg <- cv.model_grpreg$pe[cv.model_grpreg$min]  #PE
               
               
               result <- multiResultClass()
               result$Runtime <- round(matrix(c(cv.process.time1, cv.process.time2), nrow = 2), digits = 3)
               result$RME <- round(matrix(c(RME.grp_lasso, RME.grpreg), nrow = 2), digits = 4)
               result$RMSE <- round(matrix(c(RMSE.grp_lasso, RMSE.grpreg), nrow = 2), digits = 4)
               result$cross_entropy <- round(matrix(c(cross_entropy.grp_lasso, cross_entropy.grpreg), nrow = 2), digits = 4)
               result$wrong_prediction_rate <- round(matrix(c(wrong.prediction.rate.grp_lasso, wrong.prediction.rate.grpreg), nrow = 2), digits = 4)
               
               return(result)
             }
  stopCluster(cl)
  
  n.data.loop <- length(data.loop)
  
  Runtime <- matrix(rep(0, 2 * n.data.loop), nrow = 2)
  rownames(Runtime) <- c("grLasso", "grpreg")
  colnames(Runtime) <- paste0("Data_", data.loop)
  
  RME <- matrix(rep(0, 2 * n.data.loop), nrow = 2)
  rownames(RME) <- c("grLasso", "grpreg")
  colnames(RME) <- paste0("Data_", data.loop)
  
  RMSE <- matrix(rep(0, 2 * n.data.loop), nrow = 2)
  rownames(RMSE) <- c("grLasso", "grpreg")
  colnames(RMSE) <- paste0("Data_", data.loop)
  
  cross_entropy <- matrix(rep(0, 2 * n.data.loop), nrow = 2)
  rownames(cross_entropy) <- c("grLasso", "grpreg")
  colnames(cross_entropy) <- paste0("Data_", data.loop)
  
  wrong_prediction_rate <- matrix(rep(0, 2 * n.data.loop), nrow = 2)
  rownames(wrong_prediction_rate) <- c("grLasso", "grpreg")
  colnames(wrong_prediction_rate) <- paste0("Data_", data.loop)
  
  for (i in data.loop){
    Runtime[, i] <- Model.Comparison[[i]]$Runtime
    RME[, i] <- Model.Comparison[[i]]$RME
    RMSE[, i] <- Model.Comparison[[i]]$RMSE
    cross_entropy[, i] <- Model.Comparison[[i]]$cross_entropy
    wrong_prediction_rate[, i] <- Model.Comparison[[i]]$wrong_prediction_rate
  }
  
  Runtime.Mean[, ind] <- round(apply(Runtime, 1, mean), digits = 3)
  Runtime.sd[, ind] <- round(apply(Runtime, 1, sd), digits = 3)
  
  RME.Mean[, ind] <- round(apply(RME, 1, mean), digits = 3)
  RME.sd[, ind] <- round(apply(RME, 1, sd), digits = 3)
  
  RMSE.Mean[, ind] <- round(apply(RMSE, 1, mean), digits = 3)
  RMSE.sd[, ind] <- round(apply(RMSE, 1, sd), digits = 3)
  
  cross_entropy.Mean[, ind] <- round(apply(cross_entropy, 1, mean), digits = 3)
  cross_entropy.sd[, ind] <- round(apply(cross_entropy, 1, sd), digits = 3)
  
  wrong_prediction_rate.Mean[, ind] <- round(apply(wrong_prediction_rate, 1, mean), digits = 3)
  wrong_prediction_rate.sd[, ind] <- round(apply(wrong_prediction_rate, 1, sd), digits = 3)
  
}


### Figures (Runtime, RME, RMSE, cross_entropy, wrong_prediction_rate)
library(reshape2)

#### Figure1: Runtime

n.prov <- m.sequence
Runtime.Mean.lower <- Runtime.Mean - Runtime.sd
Runtime.Mean.upper <- Runtime.Mean + Runtime.sd


Runtime.figure.df <- data.frame("model" = melt(Runtime.Mean)$Var1,
                                "n.prov" = rep(n.prov, each = 2),
                                "Mean" = melt(Runtime.Mean)$value,
                                "lower" = melt(Runtime.Mean.lower)$value,
                                "upper" = melt(Runtime.Mean.upper)$value)

Runtime.plot <- ggplot(Runtime.figure.df, aes(n.prov, group = factor(model))) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(model)), alpha = 0.5) +
  geom_line(aes(y = Mean, color = factor(model), linetype = factor(model)), size = 1)  + 
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(size = 13, face="bold", family = "serif"),
        axis.title = element_text(size = 12, family = "serif"),
        plot.caption = element_text(size = 10, face = "italic", family = "serif"),
        legend.text = element_text(size = 13, family = "serif")) + 
  theme(axis.text = element_text(face = "italic", size = 9, family = "serif")) + 
  theme(legend.position = c(0.12, 0.88)) + 
  theme(legend.key.height= unit(0.5, 'cm'),
        legend.key.width= unit(1.5, 'cm')) +
  labs(title = "Runtime of the Entire Regularization Path (with 10-fold cross validation)", 
       x = "Number of providers", 
       y = "Runtime (mins)",
       caption = "( provider size varies from 100 to 1100 )") +
  scale_x_continuous(breaks = n.prov) + 
  scale_linetype_manual(values = c("solid", "dotdash"), name = "", labels = c("grLasso", "grpreg")) + 
  scale_color_manual(values = c("blue", "red"), name = "", labels = c("grLasso", "grpreg")) + 
  scale_fill_manual(values = c("grey85", "grey90"), name = "", labels = c("grLasso", "grpreg"))


#### Figure2: RMSE
n.prov <- m.sequence
RMSE.Mean.lower <- RMSE.Mean - RMSE.sd
RMSE.Mean.upper <- RMSE.Mean + RMSE.sd


RMSE.figure.df <- data.frame("model" = melt(RMSE.Mean)$Var1,
                             "n.prov" = rep(n.prov, each = 2),
                             "Mean" = melt(RMSE.Mean)$value,
                             "lower" = melt(RMSE.Mean.lower)$value,
                             "upper" = melt(RMSE.Mean.upper)$value)

RMSE.plot <- ggplot(RMSE.figure.df, aes(n.prov, group = factor(model))) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(model)), alpha = 0.5) +
  geom_line(aes(y = Mean, color = factor(model), linetype = factor(model)), size = 1)  + 
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(size = 13, face="bold", family = "serif"),
        axis.title = element_text(size = 12, family = "serif"),
        plot.caption = element_text(size = 10, face = "italic", family = "serif"),
        legend.title = element_text(size = 12, family = "serif"),
        legend.text = element_text(size = 11, family = "serif")) + 
  theme(axis.text = element_text(face = "italic", size = 9, family = "serif")) + 
  theme(legend.key.height= unit(0.5, 'cm'),
        legend.key.width= unit(1.2, 'cm')) + 
  labs(title = "Root Mean Squared Error (RMSE) Comparison", 
       x = "Number of providers", 
       y = "RMSE") +
  scale_x_continuous(breaks = n.prov) + 
  scale_linetype_manual(values = c("solid", "dotdash"), name = "", labels = c("grLasso", "grpreg")) + 
  scale_color_manual(values = c("blue", "red"), name = "", labels = c("grLasso", "grpreg")) + 
  scale_fill_manual(values = c("grey85", "grey90"), name = "", labels = c("grLasso", "grpreg"))


#### Figure3: RME
n.prov <- m.sequence
RME.Mean.lower <- RME.Mean - RME.sd
RME.Mean.upper <- RME.Mean + RME.sd


RME.figure.df <- data.frame("model" = melt(RME.Mean)$Var1,
                            "n.prov" = rep(n.prov, each = 2),
                            "Mean" = melt(RME.Mean)$value,
                            "lower" = melt(RME.Mean.lower)$value,
                            "upper" = melt(RME.Mean.upper)$value)

RME.plot <- ggplot(RME.figure.df, aes(n.prov, group = factor(model))) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(model)), alpha = 0.5) +
  geom_line(aes(y = Mean, color = factor(model), linetype = factor(model)), size = 1)  + 
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(size = 13, face="bold", family = "serif"),
        axis.title = element_text(size = 12, family = "serif"),
        plot.caption = element_text(size = 10, face = "italic", family = "serif"),
        legend.title = element_text(size = 12, family = "serif"),
        legend.text = element_text(size = 11, family = "serif")) + 
  theme(axis.text = element_text(face = "italic", size = 9, family = "serif")) + 
  theme(legend.key.height= unit(0.5, 'cm'),
        legend.key.width= unit(1.2, 'cm')) + 
  labs(title = "Root Model Error (RME) Comparison", 
       x = "Number of providers", 
       y = "RME") +
  scale_x_continuous(breaks = n.prov) + 
  scale_linetype_manual(values = c("solid", "dotdash"), name = "", labels = c("grLasso", "grpreg")) + 
  scale_color_manual(values = c("blue", "red"), name = "", labels = c("grLasso", "grpreg")) + 
  scale_fill_manual(values = c("grey85", "grey90"), name = "", labels = c("grLasso", "grpreg"))

#### Figure4: cross entropy
n.prov <- m.sequence
cross_entropy.Mean.lower <- cross_entropy.Mean - cross_entropy.sd
cross_entropy.Mean.upper <- cross_entropy.Mean + cross_entropy.sd


cross_entropy.figure.df <- data.frame("model" = melt(cross_entropy.Mean)$Var1,
                                      "n.prov" = rep(n.prov, each = 2),
                                      "Mean" = melt(cross_entropy.Mean)$value,
                                      "lower" = melt(cross_entropy.Mean.lower)$value,
                                      "upper" = melt(cross_entropy.Mean.upper)$value)

cross_entropy.plot <- ggplot(cross_entropy.figure.df, aes(n.prov, group = factor(model))) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(model)), alpha = 0.5) +
  geom_line(aes(y = Mean, color = factor(model), linetype = factor(model)), size = 1)  + 
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(size = 13, face="bold", family = "serif"),
        axis.title = element_text(size = 12, family = "serif"),
        plot.caption = element_text(size = 10, face = "italic", family = "serif"),
        legend.title = element_text(size = 12, family = "serif"),
        legend.text = element_text(size = 11, family = "serif")) + 
  theme(legend.key.height= unit(0.5, 'cm'),
        legend.key.width= unit(1.2, 'cm')) + 
  theme(axis.text = element_text(face = "italic", size = 9, family = "serif")) + 
  labs(title = "Cross-entropy Loss Comparison", 
       x = "Number of providers", 
       y = "Cross entropy") +
  scale_x_continuous(breaks = n.prov) + 
  scale_linetype_manual(values = c("solid", "dotdash"), name = "", labels = c("grLasso", "grpreg")) + 
  scale_color_manual(values = c("blue", "red"), name = "", labels = c("grLasso", "grpreg")) + 
  scale_fill_manual(values = c("grey85", "grey90"), name = "", labels = c("grLasso", "grpreg"))



#### Figure5: wrong prediction rate
n.prov <- m.sequence
wrong_prediction_rate.Mean.lower <- wrong_prediction_rate.Mean - wrong_prediction_rate.sd
wrong_prediction_rate.Mean.upper <- wrong_prediction_rate.Mean + wrong_prediction_rate.sd


wrong_prediction_rate.figure.df <- data.frame("model" = melt(wrong_prediction_rate.Mean)$Var1,
                                              "n.prov" = rep(n.prov, each = 2),
                                              "Mean" = melt(wrong_prediction_rate.Mean)$value,
                                              "lower" = melt(wrong_prediction_rate.Mean.lower)$value,
                                              "upper" = melt(wrong_prediction_rate.Mean.upper)$value)

wrong_prediction_rate.plot <- ggplot(wrong_prediction_rate.figure.df, aes(n.prov, group = factor(model))) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(model)), alpha = 0.5) +
  geom_line(aes(y = Mean, color = factor(model), linetype = factor(model)), size = 1)  + 
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(size = 13, face="bold", family = "serif"),
        axis.title = element_text(size = 12, family = "serif"),
        plot.caption = element_text(size = 10, face = "italic", family = "serif"),
        legend.title = element_text(size = 12, family = "serif"),
        legend.text = element_text(size = 11, family = "serif")) + 
  theme(legend.key.height= unit(0.5, 'cm'),
        legend.key.width= unit(1.2, 'cm')) + 
  theme(axis.text = element_text(face = "italic", size = 9, family = "serif")) + 
  labs(title = "Prediction Error Rate Comparison", 
       x = "Number of providers", 
       y = "Prediction error rate (%)") +
  scale_x_continuous(breaks = n.prov) + 
  scale_linetype_manual(values = c("solid", "dotdash"), name = "", labels = c("grLasso", "grpreg")) + 
  scale_color_manual(values = c("blue", "red"), name = "", labels = c("grLasso", "grpreg")) + 
  scale_fill_manual(values = c("grey85", "grey90"), name = "", labels = c("grLasso", "grpreg"))




######----------- With only one Provider ------------######

multiResultClass2 <- function(Runtime = NULL, max.beta.diff = NULL, max.gamma.diff = NULL){
  result <- list(
    Runtime = Runtime,
    max.beta.diff = max.beta.diff,
    max.gamma.diff = max.gamma.diff
  )
  ## Set the name for the class
  class(result) <- append(class(result), "multiResultClass2")
  return(result)
}

num.beta <- seq(4, 50, 2)

Runtime.Mean <- matrix(rep(0, length(num.beta) * 2), nrow = 2)
colnames(Runtime.Mean) <- paste0("n.beta = ", num.beta)
rownames(Runtime.Mean) <- c("MyCode", "grpreg")

Runtime.sd <- matrix(rep(0, length(num.beta) * 2), nrow = 2)
colnames(Runtime.sd) <- paste0("n.beta = ", num.beta)
rownames(Runtime.sd) <- c("MyCode", "grpreg")

MBD.Mean <- matrix(rep(0, length(num.beta) * 1), nrow = 1)
colnames(MBD.Mean) <- paste0("n.beta = ", num.beta)
rownames(MBD.Mean) <- "MBD.Mean"

MBD.sd <- matrix(rep(0, length(num.beta) * 1), nrow = 1)
colnames(MBD.sd) <- paste0("n.beta = ", num.beta)
rownames(MBD.sd) <- "MBD.sd"

MGD.Mean <- matrix(rep(0, length(num.beta) * 1), nrow = 1)
colnames(MGD.Mean) <- paste0("n.beta = ", num.beta)
rownames(MGD.Mean) <- "MGD.Mean"

MGD.sd <- matrix(rep(0, length(num.beta) * 1), nrow = 1)
colnames(MGD.sd) <- paste0("n.beta = ", num.beta)
rownames(MGD.sd) <- "MGD.sd"


data.loop <- 1:10
ind <- 0

for (j in num.beta){ #outer loop for
  ind <- ind + 1
  
  cl.cores <- 6
  cl <- makeCluster(cl.cores - 1)
  registerDoParallel(cl) 
  Model.Comparison <- 
    foreach (i = data.loop, .packages = c("grpreg", "RcppArmadillo", "MASS", "Matrix", "TmpGrlasso")) %dopar% {
      source("GroupLasso/Relevant Functions/Simulation_Functions.R")
      sim.parameters <- list(m = 1, n.beta = j, n.groups = floor(j/3), prop.NonZero.group = 1, 
                             prop.outlier = 0, rho = 0.6)
      Sim <- Simulation_data_GroupLasso(sim.parameters, unpenalized.beta = F)
      data <- Sim$sim.data
      Z.char <- colnames(data)[3:(ncol(data) - 1)]
      Y.char <- 'Y'
      original.group <- Sim$group
      
      # Use exactly the same lambda sequence
      grplasso <- grp.lasso(data, Y.char, Z.char, group = original.group, nlambda = 50,
                            tol = 1e-2, max.iter = 1e5, trace.lambda = F)
      lambda.seq <- grplasso$lambda
      
      
      start <- Sys.time()
      grplasso2 <- grp.lasso(data, Y.char, Z.char, group = original.group, lambda = lambda.seq,
                             tol = 1e-12, max.iter = 1e5, trace.lambda = F)
      end <- Sys.time()
      process.time1 <- difftime(end, start, units = 'mins')
      
      start <- Sys.time()
      grpreg2 <- grpreg(data[,c(Z.char)], data[,Y.char], family = "binomial", lambda = lambda.seq,
                        penalty = "grLasso", group = original.group, alpha = 1, eps = 1e-12, max.iter = 1e5)
      end <- Sys.time()
      process.time2 <- difftime(end, start, units = 'mins')
      
      
      grLasso.beta <- grplasso2$beta
      grLasso.gamma <- grplasso2$gamma
      
      grpreg.beta <- grpreg2$beta[2:(1 + sim.parameters$n.beta), ,drop = F]
      grpreg.gamma <- matrix(grpreg2$beta[1, ], nrow = 1)
      
      n.show.lambda <- min(ncol(grLasso.beta), ncol(grpreg.beta)) - 1  # if reach maximum iterations
      beta.diff <- grLasso.beta[, 1:n.show.lambda] - grpreg.beta[, 1:n.show.lambda]
      gamma.diff <- grLasso.gamma[, 1:n.show.lambda] - grpreg.gamma[, 1:n.show.lambda]
      
      result <- multiResultClass2()
      result$Runtime <- round(matrix(c(process.time1, process.time2), nrow = 2), digits = 3)
      result$max.beta.diff <- max(abs(beta.diff))
      result$max.gamma.diff <- max(abs(gamma.diff))
      result
      return(result)
    }
  stopCluster(cl)
  
  n.data.loop <- length(data.loop)
  
  Runtime <- matrix(rep(0, 2 * n.data.loop), nrow = 2)
  rownames(Runtime) <- c("grLasso", "grpreg")
  colnames(Runtime) <- paste0("Data_", data.loop)
  
  MBD <- matrix(rep(0, n.data.loop), nrow = 1)
  rownames(MBD) <- "MBD"
  colnames(MBD) <- paste0("Data_", data.loop)
  
  MGD <- matrix(rep(0, n.data.loop), nrow = 1)
  rownames(MGD) <- "MGD"
  colnames(MGD) <- paste0("Data_", data.loop)
  
  
  for (i in data.loop){
    Runtime[, i] <- Model.Comparison[[i]]$Runtime
    MBD[, i] <- Model.Comparison[[i]]$max.beta.diff
    MGD[, i] <- Model.Comparison[[i]]$max.gamma.diff
  }
  
  Runtime.Mean[, ind] <- round(apply(Runtime, 1, mean), digits = 3)
  Runtime.sd[, ind] <- round(apply(Runtime, 1, sd), digits = 3)
  
  MBD.Mean[, ind] <- apply(MBD, 1, mean)
  MBD.sd[, ind] <- apply(MBD, 1, sd)
  
  MGD.Mean[, ind] <- apply(MGD, 1, mean)
  MGD.sd[, ind] <- apply(MGD, 1, sd)
}


### Figures (Runtime, RME, RMSE, cross_entropy, wrong_prediction_rate)
library(reshape2)

#### Figure1: Runtime
n.beta <- num.beta
Runtime.Mean.lower <- Runtime.Mean - Runtime.sd
Runtime.Mean.upper <- Runtime.Mean + Runtime.sd


Runtime.figure.df <- data.frame("model" = melt(Runtime.Mean)$Var1,
                                "n.beta" = rep(n.beta, each = 2),
                                "Mean" = melt(Runtime.Mean)$value,
                                "lower" = melt(Runtime.Mean.lower)$value,
                                "upper" = melt(Runtime.Mean.upper)$value)

Runtime.plot <- ggplot(Runtime.figure.df, aes(n.beta, group = factor(model))) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(model)), alpha = 0.5) +
  geom_line(aes(y = Mean, color = factor(model), linetype = factor(model)), size = 1)  + 
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(size = 13, face="bold", family = "serif"),
        axis.title = element_text(size = 12, family = "serif"),
        plot.caption = element_text(size = 10, face = "italic", family = "serif"),
        legend.text = element_text(size = 13, family = "serif")) + 
  theme(axis.text = element_text(face = "italic", size = 9, family = "serif")) + 
  theme(legend.position = c(0.12, 0.88)) + 
  theme(legend.key.height= unit(0.5, 'cm'),
        legend.key.width= unit(1.5, 'cm')) +
  labs(title = "Runtime of the Entire Regularization Path (with Single Provider)", 
       x = "Number of beta", 
       y = "Runtime (mins)",
       caption = "( Sample size is about 5,000; Number of beta varies from 4 to 50 )") +
  scale_x_continuous(breaks = n.beta) + 
  scale_linetype_manual(values = c("solid", "dotdash"), name = "", labels = c("MyCode", "grpreg")) + 
  scale_color_manual(values = c("blue", "red"), name = "", labels = c("MyCode", "grpreg")) + 
  scale_fill_manual(values = c("grey85", "grey90"), name = "", labels = c("MyCode", "grpreg"))


#### Figure2: maximum estimation difference
n.beta <- num.beta
MBD.Mean.lower <- MBD.Mean - MBD.sd
MBD.Mean.upper <- MBD.Mean + MBD.sd
MGD.Mean.lower <- MGD.Mean - MGD.sd
MGD.Mean.upper <- MGD.Mean + MGD.sd


MBD.figure.df <- data.frame("model" = melt(MBD.Mean)$Var1,
                            "n.beta" = rep(n.beta, each = 1),
                            "Mean" = melt(MBD.Mean)$value,
                            "lower" = melt(MBD.Mean.lower)$value,
                            "upper" = melt(MBD.Mean.upper)$value)
MGD.figure.df <- data.frame("model" = melt(MGD.Mean)$Var1,
                            "n.beta" = rep(n.beta, each = 1),
                            "Mean" = melt(MGD.Mean)$value,
                            "lower" = melt(MGD.Mean.lower)$value,
                            "upper" = melt(MGD.Mean.upper)$value)
esti.figure <- rbind(MBD.figure.df, MGD.figure.df)


MBD.plot <- ggplot(esti.figure, aes(n.beta, group = factor(model))) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(model)), alpha = 0.5) +
  geom_line(aes(y = Mean, color = factor(model), linetype = factor(model)), size = 1)  + 
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(size = 13, face="bold", family = "serif"),
        axis.title = element_text(size = 12, family = "serif"),
        plot.caption = element_text(size = 10, face = "italic", family = "serif"),
        legend.text = element_text(size = 13, family = "serif")) + 
  theme(axis.text = element_text(face = "italic", size = 9, family = "serif")) + 
  theme(legend.position = c(0.12, 0.88)) + 
  theme(legend.key.height= unit(0.5, 'cm'),
        legend.key.width= unit(1.2, 'cm')) + 
  labs(title = "Maximum Difference in Estimation (with Single Provider)", 
       x = "Number of beta", 
       y = "Max Difference",
       caption = "( Sample size is about 5,000; Number of beta varies from 4 to 50 )") +
  scale_x_continuous(breaks = n.beta) + 
  scale_linetype_manual(values = c("solid", "dotdash"), name = "", labels = c("beta", "intercept")) + 
  scale_color_manual(values = c("blue", "red"), name = "", labels = c("beta", "intercept")) + 
  scale_fill_manual(values = c("grey85", "grey90"), name = "", labels = c("beta", "intercept")) +
  ylim(-1e-9, 1e-9)