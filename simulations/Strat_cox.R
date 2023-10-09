load("data/Cox_Data.rda")
data(Cox_Data)
data <- as.data.frame(Cox_Data$data)
Event.char <- Cox_Data$Event.char
Z.char <- Cox_Data$Z.char
Time.char <- Cox_Data$Time.char
prov.char <- Cox_Data$prov.char


#(1) stratified cox model with no penalty terms; compare with survival package
library(survival)
Surv.model <- coxph(Surv(time, status) ~ Z1 + Z2 + Z3 + Z4 + Z5 + strata(Prov.ID), data = data)
summary(Surv.model)$coefficients[, 1]
##          Z1          Z2          Z3          Z4           Z5 
## -0.33364676 -0.01804808  0.03927048  0.09085354  -0.05872726 

library(ppLasso)
fit.pplasso.lambda0 <- Strat.cox(data, Event.char, Z.char, Time.char, prov.char, lambda = 0, tol = 1e-16)
fit.pplasso.lambda0$beta
##              0
## Z1 -0.33364676
## Z2 -0.01804808
## Z3  0.03927048
## Z4  0.09085354
## Z5 -0.05872726



#(2) cox model with penalty terms; compare with grpreg package
library(grpreg)
X <- data[, Z.char]
y <- Surv(data[,Time.char], data[,Event.char])

## (2.1) grpreg
start.time <- Sys.time()
fit.grpreg <- grpsurv(X, y, eps = 1e-7, max.iter = 1e7)
end.time <- Sys.time()
process.grpreg <- difftime(end.time, start.time, units = 'sec')
process.grpreg #Time difference of 2.376379 secs

fit.grpreg$lambda[1:10]
# [1] 0.17446120 0.16270310 0.15173746 0.14151086 0.13197350 0.12307893
# [7] 0.11478382 0.10704778 0.09983311 0.09310470
fit.grpreg$beta[, c(1:3, 98:100)]
#     0.1745     0.1627      0.1517       2e-04       2e-04       2e-04
#  Z1      0 0.00000000 0.000000000 -0.20147362 -0.20177717 -0.20206026
#  Z2      0 0.00000000 0.000000000  0.10224548  0.10233260  0.10241386
#  Z3      0 0.01554764 0.025152982  0.16221528  0.16230329  0.16238537
#  Z4      0 0.00000000 0.005436268  0.12490637  0.12498693  0.12506207
#  Z5      0 0.00000000 0.000000000  0.06119763  0.06126359  0.06132512

## (2.2) pplasso
start.time <- Sys.time()
fit.pplasso.no_prov <- Strat.cox(data, Event.char, Z.char, Time.char, tol = 1e-7)
end.time <- Sys.time()
process.pplasso <- difftime(end.time, start.time, units = 'sec')
process.pplasso #Time difference of 2.420852 secs

fit.pplasso.no_prov$lambda[1:10]
# [1] 0.17447120 0.16270310 0.15173746 0.14151086 0.13197350 0.12307893
# [7] 0.11478382 0.10704778 0.09983311 0.09310470

fit.pplasso.no_prov$beta[, c(1:3, 98:100)]
#     0.1745     0.1627      0.1517       2e-04       2e-04       2e-04
#  Z1      0 0.00000000 0.000000000 -0.20147362 -0.20177717 -0.20206026
#  Z2      0 0.00000000 0.000000000  0.10224548  0.10233260  0.10241386
#  Z3      0 0.01554764 0.025152982  0.16221528  0.16230329  0.16238537
#  Z4      0 0.00000000 0.005436268  0.12490637  0.12498693  0.12506207
#  Z5      0 0.00000000 0.000000000  0.06119763  0.06126359  0.06132512

## (2.3) compare with coxph()
Surv.model <- coxph(Surv(time, status) ~ Z1 + Z2 + Z3 + Z4 + Z5, data = data)
summary(Surv.model)$coefficients[, 1]
#           Z1          Z2          Z3          Z4          Z5 
#  -0.20597781  0.10353829  0.16352106  0.12610190  0.06217664 




#####-------- 3. Compare grpreg and pplasso based on large number of facilities --------#####
cox_large <- sim.cox(100, 40) ##300 facilities, 40 predictors
dim(cox_large$data) #[1] 8008   43
#head(cox_large$data)
Event.char <- cox_large$Event.char
Z.char <- cox_large$Z.char
Time.char <- cox_large$Time.char
prov.char <- cox_large$prov.char

data <- cox_large$data
X <- data[, Z.char]
y <- Surv(data[,Time.char], data[,Event.char])

## (3.1) grpreg
start.time <- Sys.time()
fit.grpreg <- grpsurv(X, y, eps = 1e-14, max.iter = 1e8)
end.time <- Sys.time()
process.grpreg <- difftime(end.time, start.time, units = 'sec')
process.grpreg #Time difference of 210.9417 secs

fit.grpreg$lambda[1:10]
# [1] 0.8306732 0.7746886 0.7224772 0.6737846 0.6283738 0.5860235 0.5465275
# [8] 0.5096934 0.4753418 0.4433053
fit.grpreg$beta[1:5, c(1:3, 98:100)]

#     0.8307     0.7747    0.7225        0.001        9e-04        8e-04
#  Z1      0 0.07562118 0.1398554  1.829699932  1.830461804  1.831172574
#  Z2      0 0.00000000 0.0000000  1.010561206  1.011017096  1.011442431
#  Z3      0 0.04966932 0.1087896  1.610528013  1.611249407  1.611922423
#  Z4      0 0.00000000 0.0000000 -0.057231231 -0.057890893 -0.058506160
#  Z5      0 0.00000000 0.0000000  0.001410265  0.001634993  0.001844608


## (3.2) pplasso
start.time <- Sys.time()
fit.pplasso.no_prov <- Strat.cox(data, Event.char, Z.char, Time.char, tol = 1e-14,
                                 max.each.iter = 50000)
end.time <- Sys.time()
process.pplasso <- difftime(end.time, start.time, units = 'sec')
process.pplasso #Time difference of 210.9398 secs

fit.pplasso.no_prov$lambda[1:10]
# [1] 0.8306832 0.7746886 0.7224772 0.6737846 0.6283738 0.5860235 0.5465275
# [8] 0.5096934 0.4753418 0.4433053

fit.pplasso.no_prov$beta[1:5, c(1:3, 98:100)]
#     0.8307     0.7747    0.7225        0.001        9e-04        8e-04
#  Z1      0 0.07562118 0.1398554  1.829699932  1.830461804  1.831172574
#  Z2      0 0.00000000 0.0000000  1.010561206  1.011017096  1.011442431
#  Z3      0 0.04966932 0.1087896  1.610528013  1.611249407  1.611922423
#  Z4      0 0.00000000 0.0000000 -0.057231231 -0.057890893 -0.058506160
#  Z5      0 0.00000000 0.0000000  0.001410265  0.001634993  0.001844608





