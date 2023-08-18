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

fit.pplasso.lambda0 <- Strat.cox(data, Event.char, Z.char, Time.char, prov.char, lambda = 0)
fit.pplasso.lambda0$beta
##              0
## Z1 -0.33310397
## Z2 -0.01879292
## Z3  0.03807474
## Z4  0.09036853
## Z5 -0.05824560



#(2) cox model with penalty terms; compare with grpreg package
library(grpreg)
library(survival)
X <- data[, Z.char]
y <- Surv(data[,Time.char], data[,Event.char])
fit.grpreg <- grpsurv(X, y)
fit.grpreg$lambda[1:10]
# [1] 0.17446120 0.16270310 0.15173746 0.14151086 0.13197350 0.12307893
# [7] 0.11478382 0.10704778 0.09983311 0.09310470
fit.grpreg$beta[, c(1:3, 98:100)]
#     0.1745     0.1627      0.1517   ...       2e-04       2e-04       2e-04
#  Z1      0 0.00000000  0.00000000   ... -0.19624072 -0.19680026 -0.19734155
#  Z2      0 0.00000000  0.00000000   ...  0.10515659  0.10523349  0.10529852
#  Z3      0 0.01552398  0.02854876   ...  0.16216819  0.16233325  0.16249261
#  Z4      0 0.00000000  0.00235295   ...  0.12227949  0.12245610  0.12263064
#  Z5      0 0.00000000  0.00000000   ...  0.05634573  0.05649829  0.05665031

fit.pplasso.no_prov <- Strat.cox(data, Event.char, Z.char, Time.char)
fit.pplasso.no_prov$lambda[1:10]
# [1] 0.17447120 0.16270310 0.15173746 0.14151086 0.13197350 0.12307893
# [7] 0.11478382 0.10704778 0.09983311 0.09310470

fit.pplasso.no_prov$beta[, c(1:3, 98:100)]
#     0.1745    0.1627       0.1517   ...       2e-04       2e-04       2e-04
#  Z1      0 0.0000000  0.000000000   ... -0.20085083 -0.20114023 -0.20142806
#  Z2      0 0.0000000  0.000000000   ...  0.10257059  0.10264942  0.10272778
#  Z3      0 0.0155369  0.025755673   ...  0.16218908  0.16227237  0.16235487
#  Z4      0 0.0000000  0.004892787   ...  0.12459976  0.12467795  0.12475556
#  Z5      0 0.0000000  0.000000000   ...  0.06065316  0.06071941  0.06078516


Surv.model <- coxph(Surv(time, status) ~ Z1 + Z2 + Z3 + Z4 + Z5, data = data)
summary(Surv.model)$coefficients[, 1]
#           Z1          Z2          Z3          Z4          Z5 
#  -0.20597781  0.10353829  0.16352106  0.12610190  0.06217664 



