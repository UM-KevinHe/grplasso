library(RcppArmadillo)
library(Rcpp)
library(mvtnorm)
library(discSurv)
library(matrixStats)
library(survival)
library(TmpDiscSurv)
library(fastDummies)
library(grpreg)

setwd("/home/ybshao/Discrete Survival")

#source("ppSurv.R")
#source('Discrete_survival.R')
#sourceCpp("/home/ybshao/GroupLasso/Relevant Functions/Relevant_Functions.cpp")
#source("/home/ybshao/GroupLasso/Relevant Functions/Relevant_Functions.R")

#sourceCpp("ppSurv.cpp")
#sourceCpp("NR_residuals.cpp")
#sourceCpp("Discrete_logit_NR_Di.cpp")

set.seed(100)
n <- 1000
#p <- 100 #total number of covariate
p <- 10
cens_upper <- 10
Z.char <- paste0('Z', 1:p)
#beta <- c(round(runif(10, -2, 2), digits = 3), rep(0, p - 10))
beta <- c(round(runif(p, -2, 2), digits = 3))
n.days <- 10
# assume baseline hazard follows weibull(lambda, gamma) distribution:
lam <- 0.1
ga <- 1.2
baseline.hazard <- lam * ga * (lam * 1:n.days)^(ga - 1)
day_effect <- log(baseline.hazard)/(1 - baseline.hazard)
data <- sim.disc(beta, day_effect, n, Z.char, cens_upper)
#length(table(data$time))
#f.c <- colSums(table(data$time, data$status))
table(data$time, data$status)
#f.c[1]/(sum(f.c))

#dim(data)
#table(data$time)

## use "pplasso" to solve logit-link discrete survival model  
##library(TmpLasso)
df_long <- dataLong(dataShort = data, timeColumn = "time", eventColumn = "status", timeAsFactor = TRUE)
df_long_re <- dplyr::select(df_long, -c("time", "status", "obj")) #y, Z is all the information we need for estimating beta and gamma
Y.char <- "y"
time.char <- "timeInt"
Event.char <- "status"
Time.char <- "time"

#dim(df_long_re)

#start <- Sys.time()
use.ppLasso <- pp.lasso(df_long_re, Y.char, Z.char, time.char, MM = T,
                        backtrack = T)
#use.ppLasso$beta

use.ppSurv.MM <- pp.Surv(data, Event.char, Z.char, Time.char, MM = F,
                        backtrack = T, returnX = F)
#use.ppSurv.MM$beta

Z.2 <- dplyr::select(data, -c("time", "status"))  #covariate matrix
NoPenalized.SurvModel <- discreSurv_logit(data$time, Z.2, data$status)

# use grpreg
dummy_data <- dummy_cols(df_long_re, select_columns = time.char, remove_selected_columns = TRUE, 
                         remove_first_dummy = TRUE)
ID.char <- colnames(dummy_data)[(length(Z.char) + 2):ncol(dummy_data)]

model_grpreg <- grpreg(dummy_data[,c(Z.char, ID.char)], dummy_data[,Y.char], family = "binomial",
                       penalty = "grLasso", group = c(1:length(Z.char), rep(0, length(ID.char))), 
                       alpha = 1)
grpreg.beta <- model_grpreg$beta[2:(1 + length(Z.char)),]

if (ncol(model_grpreg$beta) == 1){
  grpreg.beta <- matrix(grpreg.beta, ncol = 1)
  grpreg.gamma <- c(model_grpreg$beta[1,],
                    model_grpreg$beta[1,] +
                      model_grpreg$beta[(length(Z.char) + 2):length(model_grpreg$beta[, 1]),])
  names(grpreg.gamma) <-  paste0("gamma_", 1:n.days)
  grpreg.gamma <- matrix(grpreg.gamma, ncol = 1)
} else {
  grpreg.gamma <- rbind(model_grpreg$beta[1,],
                        model_grpreg$beta[1,] +
                          model_grpreg$beta[(length(Z.char) + 2):length(model_grpreg$beta[, 1]),])
  rownames(grpreg.gamma) <- paste0("gamma_", 1:n.days)
}


beta.df <- cbind(use.ppSurv.MM$beta, 
                 NoPenalized.SurvModel$beta_v,
                 grpreg.beta)
colnames(beta.df) <- c("ppSurv", "Di (glm)", "grpreg")
beta.df


## check entire lambda sequence
use.ppSurv.MM2 <- pp.Surv(data, Event.char, Z.char, Time.char, MM = F,
                         backtrack = T, returnX = T, standardize = F)
ppSurv.beta2 <- use.ppSurv.MM2$beta[, c(1:3, 98:100)]
ppSurv.returnX <- use.ppSurv.MM2$returnX

model_grpreg2 <- grpreg(dummy_data[,c(Z.char, ID.char)], dummy_data[,Y.char], family = "binomial",
                       penalty = "grLasso", group = c(1:length(Z.char), rep(0, length(ID.char))), 
                       alpha = 1, standardize = T, returnX = T)
grpreg.beta2 <- model_grpreg2$beta[2:(1 + length(Z.char)),][, c(1:3, 98:100)]
grpreg.returnX <- model_grpreg2$XG

## grpreg和ppSurv做标准化使用的数据不一样，
## 前者是展开后的数据，后者是原始数据，因此mean和std不一样，进而使“中心化”之后的结果不一样
## 即使不用standardize




