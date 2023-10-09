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
fit.pplasso.lambda0 <- Strat.cox(data, Event.char, Z.char, Time.char, prov.char, lambda = 0, tol = 1e-10)
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
fit.grpreg <- grpsurv(X, y, max.iter = 1e7, eps = 1e-6)
end.time <- Sys.time()
process.grpreg <- difftime(end.time, start.time, units = 'sec')
process.grpreg #Time difference of 0.362515 secs

fit.grpreg$lambda[1:10]
# [1] 0.17446120 0.16270310 0.15173746 0.14151086 0.13197350 0.12307893
# [7] 0.11478382 0.10704778 0.09983311 0.09310470
fit.grpreg$beta[, c(1:3, 98:100)]
#     0.1745     0.1627      0.1517       2e-04       2e-04       2e-04
#  Z1      0 0.00000000 0.000000000 -0.20148653 -0.20178960 -0.20207254
#  Z2      0 0.00000000 0.000000000  0.10224314  0.10233088  0.10241326
#  Z3      0 0.01554718 0.025159162  0.16222508  0.16231301  0.16239556
#  Z4      0 0.00000000 0.005430695  0.12491628  0.12499623  0.12507075
#  Z5      0 0.00000000 0.000000000  0.06119388  0.06125949  0.06132004

## (2.2) pplasso
start.time <- Sys.time()
fit.pplasso.no_prov <- Strat.cox(data, Event.char, Z.char, Time.char, tol = 1e-6)
end.time <- Sys.time()
process.pplasso <- difftime(end.time, start.time, units = 'sec')
process.pplasso #Time difference of 0.42958 secs

fit.pplasso.no_prov$lambda[1:10]
# [1] 0.17447120 0.16270310 0.15173746 0.14151086 0.13197350 0.12307893 0.11478382
# [8] 0.10704778 0.09983311 0.09310470

fit.pplasso.no_prov$beta[, c(1:3, 98:100)]
#     0.1745     0.1627      0.1517       2e-04       2e-04       2e-04
#  Z1      0 0.00000000 0.000000000 -0.20148653 -0.20178960 -0.20207254
#  Z2      0 0.00000000 0.000000000  0.10224314  0.10233088  0.10241326
#  Z3      0 0.01554718 0.025159335  0.16222508  0.16231301  0.16239556
#  Z4      0 0.00000000 0.005430539  0.12491628  0.12499623  0.12507075
#  Z5      0 0.00000000 0.000000000  0.06119388  0.06125949  0.06132004

## (2.3) compare with coxph()
Surv.model <- coxph(Surv(time, status) ~ Z1 + Z2 + Z3 + Z4 + Z5, data = data)
summary(Surv.model)$coefficients[, 1]
#           Z1          Z2          Z3          Z4          Z5 
#  -0.20597781  0.10353829  0.16352106  0.12610190  0.06217664 


## (2.4) fit a penalized stratified cox model with "Strat.cox"
start.time <- Sys.time()
fit.pplasso.prov <- Strat.cox(data, Event.char, Z.char, Time.char, prov.char, tol = 1e-6)
end.time <- Sys.time()
process.pplasso <- difftime(end.time, start.time, units = 'sec')
process.pplasso #Time difference of 0.3507521 secs

fit.pplasso.prov$lambda[1:10]
# [1] 0.03944909 0.03678103 0.03430211 0.03199026 0.02983423 0.02782350 0.02594829
# [8] 0.02419946 0.02256850 0.02104746

fit.pplasso.prov$beta[, c(1:3, 98:100)]
#     0.0394      0.0368      0.0343           0           0           0
#  Z1      0 -0.02052363 -0.03967758 -0.33307298 -0.33311138 -0.33314714
#  Z2      0  0.00000000  0.00000000 -0.01747346 -0.01751091 -0.01754578
#  Z3      0  0.00000000  0.00000000  0.03850788  0.03855982  0.03860821
#  Z4      0  0.00000000  0.00000000  0.09009458  0.09014486  0.09019166
#  Z5      0  0.00000000  0.00000000 -0.05819218 -0.05822888 -0.05826310




#####-------- 3. Compare grpreg and pplasso based on large number of facilities --------#####
set.seed(1)
cox_large <- sim.cox(100, 40) ##300 facilities, 40 predictors
dim(cox_large$data) #[1] 8093   43
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
process.grpreg #Time difference of 246.3346 secs

fit.grpreg$lambda[1:10]
# [1] 0.8154473 0.7604889 0.7092345 0.6614344 0.6168560 0.5752820 0.5365099
# [8] 0.5003509 0.4666289 0.4351797
fit.grpreg$beta[1:5, c(1:3, 98:100)]
#     0.8154    0.7605    0.7092      9e-04      9e-04      8e-04
# Z1       0 0.0000000 0.0000000 0.82989803 0.83027592 0.83063472
# Z2       0 0.0000000 0.0000000 1.26807622 1.26860363 1.26910344
# Z3       0 0.0000000 0.0000000 0.93803943 0.93849742 0.93893008
# Z4       0 0.1297449 0.2516932 2.97848497 2.97961782 2.98067684
# Z5       0 0.0000000 0.0000000 0.04962957 0.04976157 0.04989158


## (3.2) pplasso
start.time <- Sys.time()
fit.pplasso.no_prov <- Strat.cox(data, Event.char, Z.char, Time.char, tol = 1e-14,
                                 max.each.iter = 50000)
end.time <- Sys.time()
process.pplasso <- difftime(end.time, start.time, units = 'sec')
process.pplasso #Time difference of 247.2532 secs

fit.pplasso.no_prov$lambda[1:10]
# [1] 0.8154573 0.7604889 0.7092345 0.6614344 0.6168560 0.5752820 0.5365099
# [8] 0.5003509 0.4666289 0.4351797

fit.pplasso.no_prov$beta[1:5, c(1:3, 98:100)]
#     0.8155    0.7605    0.7092      9e-04      9e-04      8e-04
# Z1       0 0.0000000 0.0000000 0.82989803 0.83027592 0.83063472
# Z2       0 0.0000000 0.0000000 1.26807622 1.26860363 1.26910344
# Z3       0 0.0000000 0.0000000 0.93803943 0.93849742 0.93893008
# Z4       0 0.1297449 0.2516932 2.97848497 2.97961782 2.98067684
# Z5       0 0.0000000 0.0000000 0.04962957 0.04976157 0.04989158





