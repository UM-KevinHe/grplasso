#####----- discrete survival without center information -----#####
load("data/Surv_Data.large_nCenter.rda")
data <- Surv_Data.large_nCenter$data
Event.char <- Surv_Data$Event.char
Z.char <- Surv_Data$Z.char
Time.char <- Surv_Data$Time.char

### use ppLasso::DiscSurv
start_time.DiscSurv <- Sys.time()
fit.DiscSurv <- DiscSurv(data, Event.char, Z.char, Time.char, standardize = T)
end_time.DiscSurv <- Sys.time()
process_time.DiscSurv <- difftime(end_time.DiscSurv, start_time.DiscSurv, units = 'secs')
fit.DiscSurv$beta[, 1:5]


### use ppLasso::pp.lasso
start_time.ppLasso <- Sys.time()
data_long <- discSurv::dataLong(dataShort = data, timeColumn = "time",
                                eventColumn = "status", timeAsFactor = TRUE)
data_long <- dplyr::select(data_long, -c("time", "status", "obj")) 
fit.pp.lasso <- pp.lasso(data_long, Y.char = "y", Z.char, standardize = T)
end_time.ppLasso <- Sys.time()
process_time.ppLasso <- difftime(end_time.ppLasso, start_time.ppLasso, units = 'secs')
fit.pp.lasso$beta[, 1:5]

### use grpreg
start_time.grpreg <- Sys.time()
data_long <- discSurv::dataLong(dataShort = data, timeColumn = "time",
                                eventColumn = "status", timeAsFactor = TRUE)
data_long <- dplyr::select(data_long, -c("time", "status", "obj")) 
fit.grpreg <- grpreg::grpreg(data_long[, Z.char], data_long[, "y"], family = "binomial",
                             penalty = "grLasso", group = c(1:length(Z.char)), alpha = 1)
end_time.grpreg <- Sys.time()
process_time.grpreg <- difftime(end_time.grpreg, start_time.grpreg, units = 'secs')
fit.grpreg$beta[, 1:5]
# note that the standardized Z matrix are different between grpreg and pplasso
# grpreg uses expanded data, while pplasso uses original data


### use glmnet
start_time.glmnet <- Sys.time()
data_long <- discSurv::dataLong(dataShort = data, timeColumn = "time",
                                eventColumn = "status", timeAsFactor = TRUE)
data_long <- dplyr::select(data_long, -c("time", "status", "obj")) 
fit.glmnet <- glmnet::glmnet(data_long[, Z.char], data_long[, "y"], family = "binomial",
                             alpha = 1)
end_time.glmnet <- Sys.time()
process_time.glmnet <- difftime(end_time.glmnet, start_time.glmnet, units = 'secs')
fit.glmnet$beta[, 1:5]

