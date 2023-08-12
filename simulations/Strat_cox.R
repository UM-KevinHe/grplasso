load("data/Cox_data.rda")
data(Cox_data)
data <- as.data.frame(Cox_data$data)
Event.char <- Cox_data$Event.char
Z.char <- Cox_data$Z.char
Time.char <- Cox_data$Time.char
prov.char <- Cox_data$prov.char


#use survival package
library(survival)
Surv.model <- coxph(Surv(time, status) ~ Z1 + Z2 + Z3 + Z4 + Z5 + strata(Prov.ID), data = data)
summary(Surv.model)
