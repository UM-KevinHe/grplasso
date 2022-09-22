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
setwd("/home/ybshao/GroupLasso")
source("Relevant Functions/Simulation_Functions.R")
source("Relevant Functions/Relevant_Functions.R")
source("grLasso_Code.R")

Rcpp::sourceCpp("Relevant Functions/Relevant_Functions.cpp")
Rcpp::sourceCpp("grLasso_Code.cpp")


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
                             trace.lambda = T, backtrack = F)
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
sim.parameters2 <- list(m = 1, n.beta = 1, n.groups = 1, prop.NonZero.group = 1, 
                        prop.outlier = 0, rho = 0.6)
set.seed(1)
## Note: If any beta are unpenalized, grpreg will add 1e-05 to the maximum lambda. 
#Sim_2 <- Simulation_data_GroupLasso(sim.parameters2, , prov.size.mean = 2000， unpenalized.beta = T, prop.unpenalized.beta = 0.5)
Sim_2 <- Simulation_data_GroupLasso(sim.parameters2, prov.size.mean = 2000, unpenalized.beta = F)
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
model_grpreg_onlyIntercept <- grpreg(data2[,c(Z.char), drop = F], data2[,Y.char], family = "binomial", penalty = "grLasso",
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


######----------- 4. Cross Validation ------------######
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


######----------- 5. simple lasso ------------######
set.seed(123)
Y.char <- 'Y'
prov.char <- 'Prov.ID'
n.beta <- 10
sim.parameters.Lasso <- list(m = 50, n.beta = n.beta, n.groups = n.beta, prop.NonZero.group = 0.2,
                             prop.outlier = 0, rho = 0.7)
Sim_Lasso <- Simulation_data_GroupLasso(sim.parameters.Lasso, unpenalized.beta = F)
data_Lasso <- Sim_Lasso$sim.data
true.beta <- Sim_Lasso$beta
true.gamma <- Sim_Lasso$gamma
Z.char <- colnames(data_Lasso)[3:(ncol(data_Lasso) - 1)]
true.mu <- Sim_Lasso$mu


# dummy data for grpreg
dummy_data <- dummy_cols(data_Lasso, select_columns = prov.char, remove_selected_columns = TRUE, 
                         remove_first_dummy = TRUE)
ID.char <- colnames(dummy_data)[(ncol(dummy_data) - sim.parameters.Lasso$m + 2):ncol(dummy_data)]


## 1. ppLasso
### (1) MM = F
start <- Sys.time()
model_pp_lasso1 <- pp.lasso(data_Lasso, Y.char, Z.char, prov.char, MM = F)
end <- Sys.time()
process.time1 <- end - start
print(process.time1)

ppLasso.beta1 <- model_pp_lasso1$beta
ppLasso.gamma1 <- model_pp_lasso1$gamma
ppLasso.lambda1 <- model_pp_lasso1$lambda
ppLasso.loss1 <- model_pp_lasso1$loss
total.iter1 <- sum(model_pp_lasso1$iter)

### (2) MM = T
start <- Sys.time()
model_pp_lasso2 <- pp.lasso(data_Lasso, Y.char, Z.char, prov.char, MM = T)
end <- Sys.time()
process.time2 <- end - start
print(process.time2)

ppLasso.beta2 <- model_pp_lasso2$beta
ppLasso.gamma2<- model_pp_lasso2$gamma
ppLasso.lambda2 <- model_pp_lasso2$lambda
ppLasso.loss2 <- model_pp_lasso2$loss
total.iter2 <- sum(model_pp_lasso2$iter)

## 2. grpreg
start <- Sys.time()
model_grpreg <- grpreg(dummy_data[ ,c(Z.char, ID.char)], dummy_data[,Y.char], family = "binomial", 
                       penalty = "grLasso", group = c(1:length(Z.char), rep(0, length(ID.char))), 
                       alpha = 1)
end <- Sys.time()
process.time3 <- end - start
print(process.time3)

grpreg.beta <- model_grpreg$beta[2:(1 + sim.parameters.Lasso$n.beta),]

if (ncol(model_grpreg$beta) == 1){
  grpreg.beta <- matrix(grpreg.beta, ncol = 1)
  grpreg.gamma <- c(model_grpreg$beta[1,],
                    model_grpreg$beta[1,] +
                      model_grpreg$beta[(sim.parameters.Lasso$n.beta + 2):length(model_grpreg$beta[, 1]),])
  names(grpreg.gamma) <-  names(table(data_Lasso$Prov.ID))
  grpreg.gamma <- matrix(grpreg.gamma, ncol = 1)
  n <- 1
} else {
  grpreg.gamma <- rbind(model_grpreg$beta[1,],
                        model_grpreg$beta[1,] +
                          model_grpreg$beta[(sim.parameters.Lasso$n.beta + 2):length(model_grpreg$beta[, 1]),])
  rownames(grpreg.gamma) <- names(table(data_Lasso$Prov.ID))
}

grpreg.lambda <- model_grpreg$lambda
grpreg.loss <- model_grpreg$loss


## 3. glmnet
start <- Sys.time()
model_glmnet <- glmnet(dummy_data[ ,c(Z.char, ID.char)], dummy_data[,Y.char], family = "binomial", 
                       penalty.factor = c(rep(1, length(Z.char)), rep(0, length(ID.char))), 
                       alpha = 1)
end <- Sys.time()
process.time4 <- end - start
print(process.time4)

glmnet.beta <- as.matrix(model_glmnet$beta[1:(sim.parameters.Lasso$n.beta),])

## Estimation difference
# pplasso with MM
n <- min(ncol(ppLasso.beta1), ncol(grpreg.beta)) - 1
beta.diff1 <- round(ppLasso.beta1[, 1:n] - grpreg.beta[, 1:n], digits = 4)
# max(abs(beta.diff1))
gamma.diff1 <- round(ppLasso.gamma1[, 1:n] - grpreg.gamma[, 1:n], digits = 4)

# pplasso without MM
n <- min(ncol(ppLasso.beta2), ncol(grpreg.beta)) - 1
beta.diff2 <- round(ppLasso.beta2[, 1:n] - grpreg.beta[, 1:n], digits = 4)
# max(abs(beta.diff2))
gamma.diff2 <- round(ppLasso.gamma2[, 1:n] - grpreg.gamma[, 1:n], digits = 4)




######----------- 6.Figure: Entire Regularization Path ------------######

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

######----------- 7.Figure: Cross-validation error at entire lambda sequence (Based on Real Data)------------######
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


######----------- 8.Real Data Estimation: (Based on "Birthwt") ------------######
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






######----------- 9.Time Comparison: GroupLasso ------------######
# n.beta = 50, n.groups = 4, n.provider range from 50 to 500
multiResultClass <- function(Runtime = NULL){#, RME = NULL, RMSE = NULL, cross_entropy = NULL, wrong_prediction_rate = NULL){
  result <- list(
    Runtime = Runtime#,
    #RME = RME,
    #RMSE = RMSE#,
    #cross_entropy = cross_entropy,
    #wrong_prediction_rate = wrong_prediction_rate
  )
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

#RME.Mean <- matrix(rep(0, length(m.sequence) * 2), nrow = 2)
#colnames(RME.Mean) <- paste0("n.prov = ", m.sequence)
#rownames(RME.Mean) <- c("grLasso", "grpreg")

#RME.sd <- matrix(rep(0, length(m.sequence) * 2), nrow = 2)
#colnames(RME.sd) <- paste0("n.prov = ", m.sequence)
#rownames(RME.sd) <- c("grLasso", "grpreg")

#RMSE.Mean <- matrix(rep(0, length(m.sequence) * 2), nrow = 2)
#colnames(RMSE.Mean) <- paste0("n.prov = ", m.sequence)
#rownames(RMSE.Mean) <- c("grLasso", "grpreg")

#RMSE.sd <- matrix(rep(0, length(m.sequence) * 2), nrow = 2)
#colnames(RMSE.sd) <- paste0("n.prov = ", m.sequence)
#rownames(RMSE.sd) <- c("grLasso", "grpreg")

#cross_entropy.Mean <- matrix(rep(0, length(m.sequence) * 2), nrow = 2)
#colnames(cross_entropy.Mean) <- paste0("n.prov = ", m.sequence)
#rownames(cross_entropy.Mean) <- c("grLasso", "grpreg")

#cross_entropy.sd <- matrix(rep(0, length(m.sequence) * 2), nrow = 2)
#colnames(cross_entropy.sd) <- paste0("n.prov = ", m.sequence)
#rownames(cross_entropy.sd) <- c("grLasso", "grpreg")

#wrong_prediction_rate.Mean <- matrix(rep(0, length(m.sequence) * 2), nrow = 2)
#colnames(wrong_prediction_rate.Mean) <- paste0("n.prov = ", m.sequence)
#rownames(wrong_prediction_rate.Mean) <- c("grLasso", "grpreg")

#wrong_prediction_rate.sd <- matrix(rep(0, length(m.sequence) * 2), nrow = 2)
#colnames(wrong_prediction_rate.sd) <- paste0("n.prov = ", m.sequence)
#rownames(wrong_prediction_rate.sd) <- c("grLasso", "grpreg")


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
      cv.model_grp_lasso <- cv.grp.lasso(data_prep, Y.char, Z.char, prov.char, group = group, trace.lambda = T,
                                         nfolds = 10, trace.cv = F)
      end <- Sys.time()
      cv.process.time1 <- difftime(end, start, units = 'mins')  #runtime
      #cv_BestModel_grp_lasso <- cv.model_grp_lasso$fit
      #best.beta.grp_lasso <- cv_BestModel_grp_lasso$beta[, cv.model_grp_lasso$min]
      #best.eta.grp_lasso <- cv_BestModel_grp_lasso$linear.predictors[, cv.model_grp_lasso$min]
      
      #RME.grp_lasso <- sqrt(sum((best.beta.grp_lasso - true.beta)^2)/sim.parameters.GrLasso$n.beta) #RME
      #RMSE.grp_lasso <- sqrt(sum((plogis(best.eta.grp_lasso) - true.mu)^2)/nrow(data_prep)) #RMSE
      
      
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
      
      #cv_BestModel_grpreg <- cv.model_grpreg$fit
      #best.beta.grpreg <- cv_BestModel_grpreg$beta[2:(1 + sim.parameters.GrLasso$n.beta), cv.model_grpreg$min]
      #best.eta.grpreg <- cv_BestModel_grpreg$linear.predictors[, cv.model_grpreg$min]
      
      #RME.grpreg <- sqrt(sum((best.beta.grpreg - true.beta)^2)/sim.parameters.GrLasso$n.beta)  #RME
      #RMSE.grpreg <- sqrt(sum((plogis(best.eta.grpreg) - true.mu)^2)/nrow(data_prep))  #RMSE
      
      
      result <- multiResultClass()
      result$Runtime <- round(matrix(c(cv.process.time1, cv.process.time2), nrow = 2), digits = 3)
      #result$RME <- round(matrix(c(RME.grp_lasso, RME.grpreg), nrow = 2), digits = 4)
      #result$RMSE <- round(matrix(c(RMSE.grp_lasso, RMSE.grpreg), nrow = 2), digits = 4)
      #result$cross_entropy <- round(matrix(c(cross_entropy.grp_lasso, cross_entropy.grpreg), nrow = 2), digits = 4)
      #result$wrong_prediction_rate <- round(matrix(c(wrong.prediction.rate.grp_lasso, wrong.prediction.rate.grpreg), nrow = 2), digits = 4)
      
      return(result)
    }
  stopCluster(cl)
  
  n.data.loop <- length(data.loop)
  
  Runtime <- matrix(rep(0, 2 * n.data.loop), nrow = 2)
  rownames(Runtime) <- c("grLasso", "grpreg")
  colnames(Runtime) <- paste0("Data_", data.loop)
  
  #RME <- matrix(rep(0, 2 * n.data.loop), nrow = 2)
  #rownames(RME) <- c("grLasso", "grpreg")
  #colnames(RME) <- paste0("Data_", data.loop)
  
  #RMSE <- matrix(rep(0, 2 * n.data.loop), nrow = 2)
  #rownames(RMSE) <- c("grLasso", "grpreg")
  #colnames(RMSE) <- paste0("Data_", data.loop)
  
  #cross_entropy <- matrix(rep(0, 2 * n.data.loop), nrow = 2)
  #rownames(cross_entropy) <- c("grLasso", "grpreg")
  #colnames(cross_entropy) <- paste0("Data_", data.loop)
  
  #wrong_prediction_rate <- matrix(rep(0, 2 * n.data.loop), nrow = 2)
  #rownames(wrong_prediction_rate) <- c("grLasso", "grpreg")
  #colnames(wrong_prediction_rate) <- paste0("Data_", data.loop)
  
  for (i in data.loop){
    Runtime[, i] <- Model.Comparison[[i]]$Runtime
    #RME[, i] <- Model.Comparison[[i]]$RME
    #RMSE[, i] <- Model.Comparison[[i]]$RMSE
    #cross_entropy[, i] <- Model.Comparison[[i]]$cross_entropy
    #wrong_prediction_rate[, i] <- Model.Comparison[[i]]$wrong_prediction_rate
  }
  
  Runtime.Mean[, ind] <- round(apply(Runtime, 1, mean), digits = 3)
  Runtime.sd[, ind] <- round(apply(Runtime, 1, sd), digits = 3)
  
  #RME.Mean[, ind] <- round(apply(RME, 1, mean), digits = 3)
  #RME.sd[, ind] <- round(apply(RME, 1, sd), digits = 3)
  
  #RMSE.Mean[, ind] <- round(apply(RMSE, 1, mean), digits = 3)
  #RMSE.sd[, ind] <- round(apply(RMSE, 1, sd), digits = 3)
  
  #cross_entropy.Mean[, ind] <- round(apply(cross_entropy, 1, mean), digits = 3)
  #cross_entropy.sd[, ind] <- round(apply(cross_entropy, 1, sd), digits = 3)
  
  #wrong_prediction_rate.Mean[, ind] <- round(apply(wrong_prediction_rate, 1, mean), digits = 3)
  #wrong_prediction_rate.sd[, ind] <- round(apply(wrong_prediction_rate, 1, sd), digits = 3)
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
  theme(plot.title = element_text(size = 13, face = "bold", family = "serif"),
        axis.title = element_text(size = 13, family = "serif"),
        plot.caption = element_text(size = 10, face = "italic", family = "serif"),
        legend.text = element_text(size = 13, family = "serif", face = "bold")) + 
  theme(axis.text = element_text(face = "italic", size = 10, family = "serif")) + 
  theme(legend.position = c(0.22, 0.82),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank()) + 
  theme(legend.key.height = unit(0.6, 'cm'),
        legend.key.width = unit(1.2, 'cm'),
        legend.spacing.y = unit(0, 'cm')) +
  labs(title = "", 
       x = "", 
       y = "") +
  scale_x_continuous(breaks = seq(50, 400, 50)) + 
  scale_linetype_manual(values = c("solid", "dotdash"), name = "", labels = c("grLasso", "grpreg")) + 
  scale_color_manual(values = c("blue", "red"), name = "", labels = c("grLasso", "grpreg")) + 
  scale_fill_manual(values = c("grey85", "grey90"), name = "", labels = c("grLasso", "grpreg"))


#### Figure2: RMSE
#n.prov <- m.sequence
#RMSE.Mean.lower <- RMSE.Mean - RMSE.sd
#RMSE.Mean.upper <- RMSE.Mean + RMSE.sd


#RMSE.figure.df <- data.frame("model" = melt(RMSE.Mean)$Var1,
#                             "n.prov" = rep(n.prov, each = 2),
#                             "Mean" = melt(RMSE.Mean)$value,
#                             "lower" = melt(RMSE.Mean.lower)$value,
#                             "upper" = melt(RMSE.Mean.upper)$value)

#RMSE.plot <- ggplot(RMSE.figure.df, aes(n.prov, group = factor(model))) +
#  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(model)), alpha = 0.5) +
#  geom_line(aes(y = Mean, color = factor(model), linetype = factor(model)), size = 1)  + 
#  theme(panel.grid = element_blank(), panel.background = element_blank(),
#        axis.line = element_line(colour = "black")) + 
#  theme(plot.title = element_text(size = 13, face="bold", family = "serif"),
#        axis.title = element_text(size = 13, family = "serif"),
#        plot.caption = element_text(size = 10, face = "italic", family = "serif"),
#        legend.title = element_text(size = 12, family = "serif"),
#        legend.text = element_text(size = 13, family = "serif", face = "bold")) + 
#  theme(axis.text = element_text(face = "italic", size = 9, family = "serif")) + 
#  theme(legend.position = c(0.85, 0.92)) + 
#  theme(legend.key.height= unit(0.5, 'cm'),
#        legend.key.width= unit(1.2, 'cm')) + 
#  labs(title = "Root Mean Squared Error (RMSE) Comparison", 
#       x = "Number of providers", 
#       y = "RMSE") +
#  scale_x_continuous(breaks = seq(50, 400, 50)) + 
#  scale_linetype_manual(values = c("solid", "dotdash"), name = "", labels = c("grLasso", "grpreg")) + 
#  scale_color_manual(values = c("blue", "red"), name = "", labels = c("grLasso", "grpreg")) + 
#  scale_fill_manual(values = c("grey85", "grey90"), name = "", labels = c("grLasso", "grpreg"))


#### Figure3: RME
#n.prov <- m.sequence
#RME.Mean.lower <- RME.Mean - RME.sd
#RME.Mean.upper <- RME.Mean + RME.sd


#RME.figure.df <- data.frame("model" = melt(RME.Mean)$Var1,
#                            "n.prov" = rep(n.prov, each = 2),
#                            "Mean" = melt(RME.Mean)$value,
#                            "lower" = melt(RME.Mean.lower)$value,
#                            "upper" = melt(RME.Mean.upper)$value)

#RME.plot <- ggplot(RME.figure.df, aes(n.prov, group = factor(model))) +
#  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(model)), alpha = 0.5) +
#  geom_line(aes(y = Mean, color = factor(model), linetype = factor(model)), size = 1)  + 
#  theme(panel.grid = element_blank(), panel.background = element_blank(),
#        axis.line = element_line(colour = "black")) + 
#  theme(plot.title = element_text(size = 13, face = "bold", family = "serif"),
#        axis.title = element_text(size = 13, family = "serif"),
#        plot.caption = element_text(size = 10, face = "italic", family = "serif"),
#        legend.title = element_text(size = 12, family = "serif"),
#        legend.text = element_text(size = 13, family = "serif", face = "bold")) + 
#  theme(axis.text = element_text(face = "italic", size = 9, family = "serif")) + 
#  theme(legend.position = c(0.85, 0.92)) +
#  theme(legend.key.height= unit(0.5, 'cm'),
#        legend.key.width= unit(1.2, 'cm')) + 
#  labs(title = "Root Model Error (RME) Comparison", 
#       x = "Number of providers", 
#       y = "RME") +
#  scale_x_continuous(breaks = seq(50, 500, 50)) + 
#  scale_linetype_manual(values = c("solid", "dotdash"), name = "", labels = c("grLasso", "grpreg")) + 
#  scale_color_manual(values = c("blue", "red"), name = "", labels = c("grLasso", "grpreg")) + 
#  scale_fill_manual(values = c("grey85", "grey90"), name = "", labels = c("grLasso", "grpreg"))

#### Figure4: cross entropy
#n.prov <- m.sequence
#cross_entropy.Mean.lower <- cross_entropy.Mean - cross_entropy.sd
#cross_entropy.Mean.upper <- cross_entropy.Mean + cross_entropy.sd


#cross_entropy.figure.df <- data.frame("model" = melt(cross_entropy.Mean)$Var1,
#                                      "n.prov" = rep(n.prov, each = 2),
#                                      "Mean" = melt(cross_entropy.Mean)$value,
#                                      "lower" = melt(cross_entropy.Mean.lower)$value,
#                                      "upper" = melt(cross_entropy.Mean.upper)$value)

#cross_entropy.plot <- ggplot(cross_entropy.figure.df, aes(n.prov, group = factor(model))) +
#  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(model)), alpha = 0.5) +
#  geom_line(aes(y = Mean, color = factor(model), linetype = factor(model)), size = 1)  + 
#  theme(panel.grid = element_blank(), panel.background = element_blank(),
#        axis.line = element_line(colour = "black")) + 
#  theme(plot.title = element_text(size = 13, face="bold", family = "serif"),
#        axis.title = element_text(size = 12, family = "serif"),
#        plot.caption = element_text(size = 10, face = "italic", family = "serif"),
#        legend.title = element_text(size = 12, family = "serif"),
#        legend.text = element_text(size = 11, family = "serif")) + 
#  theme(legend.key.height= unit(0.5, 'cm'),
#        legend.key.width= unit(1.2, 'cm')) + 
#  theme(axis.text = element_text(face = "italic", size = 9, family = "serif")) + 
#  labs(title = "Cross-entropy Loss Comparison", 
#       x = "Number of providers", 
#       y = "Cross entropy") +
#  scale_x_continuous(breaks = n.prov) + 
#  scale_linetype_manual(values = c("solid", "dotdash"), name = "", labels = c("grLasso", "grpreg")) + 
#  scale_color_manual(values = c("blue", "red"), name = "", labels = c("grLasso", "grpreg")) + 
#  scale_fill_manual(values = c("grey85", "grey90"), name = "", labels = c("grLasso", "grpreg"))



#### Figure5: wrong prediction rate
#n.prov <- m.sequence
#wrong_prediction_rate.Mean.lower <- wrong_prediction_rate.Mean - wrong_prediction_rate.sd
#wrong_prediction_rate.Mean.upper <- wrong_prediction_rate.Mean + wrong_prediction_rate.sd


#wrong_prediction_rate.figure.df <- data.frame("model" = melt(wrong_prediction_rate.Mean)$Var1,
#                                              "n.prov" = rep(n.prov, each = 2),
#                                              "Mean" = melt(wrong_prediction_rate.Mean)$value,
#                                              "lower" = melt(wrong_prediction_rate.Mean.lower)$value,
#                                              "upper" = melt(wrong_prediction_rate.Mean.upper)$value)

#wrong_prediction_rate.plot <- ggplot(wrong_prediction_rate.figure.df, aes(n.prov, group = factor(model))) +
#  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(model)), alpha = 0.5) +
#  geom_line(aes(y = Mean, color = factor(model), linetype = factor(model)), size = 1)  + 
#  theme(panel.grid = element_blank(), panel.background = element_blank(),
#        axis.line = element_line(colour = "black")) + 
#  theme(plot.title = element_text(size = 13, face="bold", family = "serif"),
#        axis.title = element_text(size = 12, family = "serif"),
#        plot.caption = element_text(size = 10, face = "italic", family = "serif"),
#        legend.title = element_text(size = 12, family = "serif"),
#        legend.text = element_text(size = 11, family = "serif")) + 
#  theme(legend.key.height= unit(0.5, 'cm'),
#        legend.key.width= unit(1.2, 'cm')) + 
#  theme(axis.text = element_text(face = "italic", size = 9, family = "serif")) + 
#  labs(title = "Prediction Error Rate Comparison", 
#       x = "Number of providers", 
#       y = "Prediction error rate (%)") +
#  scale_x_continuous(breaks = n.prov) + 
#  scale_linetype_manual(values = c("solid", "dotdash"), name = "", labels = c("grLasso", "grpreg")) + 
#  scale_color_manual(values = c("blue", "red"), name = "", labels = c("grLasso", "grpreg")) + 
#  scale_fill_manual(values = c("grey85", "grey90"), name = "", labels = c("grLasso", "grpreg"))


save(Runtime.figure.df, Runtime.plot, 
     #RMSE.figure.df, RMSE.plot, 
     #RME.figure.df, RME.plot, 
     file = paste0("Runtime_", Sys.Date(), ".RData"))






######----------- 10.Time Comparison: Lasso ------------######
# n.beta = 50, n.provider range from 50 to 500
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

Runtime.Mean <- matrix(rep(0, length(m.sequence) * 4), nrow = 4)
colnames(Runtime.Mean) <- paste0("n.prov = ", m.sequence)
rownames(Runtime.Mean) <- c("pplasso with MM", "pplasso without MM", "grpreg", "glmnet")

Runtime.sd <- matrix(rep(0, length(m.sequence) * 4), nrow = 4)
colnames(Runtime.sd) <- paste0("n.prov = ", m.sequence)
rownames(Runtime.sd) <- c("pplasso with MM", "pplasso without MM", "grpreg", "glmnet")

iter.Mean <- matrix(rep(0, length(m.sequence) * 3), nrow = 3)
colnames(iter.Mean) <- paste0("n.prov = ", m.sequence)
rownames(iter.Mean) <- c("pplasso with MM", "pplasso without MM", "grpreg")

iter.sd <- matrix(rep(0, length(m.sequence) * 3), nrow = 3)
colnames(iter.sd) <- paste0("n.prov = ", m.sequence)
rownames(iter.sd) <- c("pplasso with MM", "pplasso without MM", "grpreg")

data.loop <- 1:10
ind <- 0

for (j in m.sequence){ 
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
      start <- Sys.time()
      cv.model_grp_lasso2 <- cv.pp.lasso(data_Lasso, Y.char, Z.char, prov.char,
                                        MM = F, nfolds = 10, trace.cv = F)
      end <- Sys.time()
      cv.process.time2 <- difftime(end, start, units = 'mins')  #runtime
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
      cv.process.time3 <-  difftime(end, start, units = 'mins') #runtime
      total.iter.grpreg <- sum((cv.model_grpreg$fit)$iter)
      
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
      cv.process.time4 <-  difftime(end, start, units = 'mins') #runtime

      
      result <- multiResultClass2()
      result$Runtime <- round(matrix(c(cv.process.time1, cv.process.time2, cv.process.time3, cv.process.time4), nrow = 4), digits = 3)
      result$total_iter <- round(matrix(c(total.iter.pplasso1, total.iter.pplasso2, total.iter.grpreg), nrow = 3), digits = 4)
      return(result)
    }
  stopCluster(cl)
  
  n.data.loop <- length(data.loop)
  
  Runtime <- matrix(rep(0, 4 * n.data.loop), nrow = 4)
  rownames(Runtime) <- c("pplasso with MM", "pplasso without MM", "grpreg", "glmnet")
  colnames(Runtime) <- paste0("Data_", data.loop)
  
  iter <- matrix(rep(0, 3 * n.data.loop), nrow = 3)
  rownames(iter) <- c("pplasso with MM", "pplasso without MM", "grpreg")
  colnames(iter) <- paste0("Data_", data.loop)
  
  for (i in data.loop){
    Runtime[, i] <- Model.Comparison[[i]]$Runtime
    iter[, i] <- Model.Comparison[[i]]$total_iter
  }
  
  Runtime.Mean[, ind] <- round(apply(Runtime, 1, mean), digits = 3)
  Runtime.sd[, ind] <- round(apply(Runtime, 1, sd), digits = 3)
  
  iter.Mean[, ind] <- round(apply(iter, 1, mean), digits = 3)
  iter.sd[, ind] <- round(apply(iter, 1, sd), digits = 3)
}

### Figures (Runtime, iter)
library(reshape2)

#### Figure1: Runtime
n.prov <- m.sequence
Runtime.Mean.lower <- Runtime.Mean - Runtime.sd
Runtime.Mean.upper <- Runtime.Mean + Runtime.sd


Runtime.figure.df <- data.frame("model" = melt(Runtime.Mean)$Var1,
                                "n.prov" = rep(n.prov, each = 4),
                                "Mean" = melt(Runtime.Mean)$value,
                                "lower" = melt(Runtime.Mean.lower)$value,
                                "upper" = melt(Runtime.Mean.upper)$value)

Runtime.plot <- ggplot(Runtime.figure.df, aes(n.prov, group = factor(model))) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(model)), alpha = 0.5) +
  geom_line(aes(y = Mean, color = factor(model), linetype = factor(model)), size = 1)  + 
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(size = 13, face = "bold", family = "serif"),
        axis.title = element_text(size = 13, family = "serif"),
        plot.caption = element_text(size = 11, face = "italic", family = "serif"),
        legend.text = element_text(size = 12, family = "serif", face = "italic")) + 
  theme(axis.text = element_text(face = "italic", size = 12, family = "serif")) + 
  theme(legend.position = c(0.20, 0.85),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank()) + 
  theme(legend.key.height = unit(0.6, 'cm'),
        legend.key.width = unit(1.2, 'cm')) +
  labs(title = "Comparison of Running Speed", 
       x = "Number of providers", 
       y = "Runtime (mins)",
       caption = "( provider size varies from 50 to 400 )") +
  scale_x_continuous(breaks = m.sequence) + 
  scale_linetype_manual(values = c("solid", "dotdash", "dotted", "twodash"), name = "", labels = c("pplasso with MM", "pplasso without MM", "grpreg", "glmnet")) + 
  scale_color_manual(values = c("black", "bisque4", "red", "blue"), name = "", labels = c("pplasso with MM", "pplasso without MM", "grpreg", "glmnet")) + 
  scale_fill_manual(values = c("grey75", "grey80", "grey85", "grey90"), name = "", labels = c("pplasso with MM", "pplasso without MM", "grpreg", "glmnet"))


#### Figure2: iter
n.prov <- m.sequence
iter.Mean.lower <- iter.Mean - iter.sd
iter.Mean.upper <- iter.Mean + iter.sd


iter.figure.df <- data.frame("model" = melt(iter.Mean)$Var1,
                             "n.prov" = rep(n.prov, each = 3),
                             "Mean" = melt(iter.Mean)$value,
                             "lower" = melt(iter.Mean.lower)$value,
                             "upper" = melt(iter.Mean.upper)$value)

iter.plot <- ggplot(iter.figure.df, aes(n.prov, group = factor(model))) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = factor(model)), alpha = 0.5) +
  geom_line(aes(y = Mean, color = factor(model), linetype = factor(model)), size = 1)  + 
  theme(panel.grid = element_blank(), panel.background = element_blank(),
        axis.line = element_line(colour = "black")) + 
  theme(plot.title = element_text(size = 13, face = "bold", family = "serif"),
        axis.title = element_text(size = 13, family = "serif"),
        plot.caption = element_text(size = 11, face = "italic", family = "serif"),
        legend.text = element_text(size = 12, family = "serif", face = "italic")) + 
  theme(axis.text = element_text(face = "italic", size = 12, family = "serif")) + 
  theme(legend.position = c(0.5, 0.95),
        legend.background = element_rect(fill = "white", color = "black"),
        legend.title = element_blank()) + 
  theme(legend.key.height = unit(0.4, 'cm'),
        legend.key.width = unit(1.2, 'cm'),
        legend.direction = "horizontal",
        legend.spacing.x = unit(0.3, 'cm')) + 
  labs(title = "Comparison of Total Iteration Numbers", 
       x = "Number of providers", 
       y = "Total #iterations",
       caption = "( provider size varies from 50 to 400 )") +
  scale_x_continuous(breaks = m.sequence) + 
  scale_linetype_manual(values = c("solid", "dotdash", "dotted"), name = "", labels = c("pplasso with MM", "pplasso without MM", "grpreg")) + 
  scale_color_manual(values = c("bisque4", "blue", "red"), name = "", labels = c("pplasso with MM", "pplasso without MM", "grpreg")) + 
  scale_fill_manual(values = c("grey80", "grey85", "grey90"), name = "", labels = c("pplasso with MM", "pplasso without MM", "grpreg")) +
  ylim(0, 8000)
  
  
setwd("/home/ybshao/GroupLasso/Simulations/Figures&Tables/Runtime/Lasso")
save(Runtime.figure.df, Runtime.plot, 
     iter.figure.df, iter.plot,
     file = paste0("Runtime_iter_RME", Sys.Date(), ".RData"))

######----------- 11. With only one Provider ------------######
multiResultClass3 <- function(Runtime = NULL, max.beta.diff = NULL, max.gamma.diff = NULL){
  result <- list(
    Runtime = Runtime,
    max.beta.diff = max.beta.diff,
    max.gamma.diff = max.gamma.diff
  )
  ## Set the name for the class
  class(result) <- append(class(result), "multiResultClass3")
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
      sim.parameters <- list(m = 1, n.beta = j, n.groups = floor(j/3), prop.NonZero.group = 1, 
                             prop.outlier = 0, rho = 0.6)
      Sim <- Simulation_data_GroupLasso(sim.parameters, prov.size.mean = 2000, unpenalized.beta = F)
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
      
      result <- multiResultClass3()
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



