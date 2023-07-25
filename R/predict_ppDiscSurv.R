#' Predictions of a ppDiscSurv object
#'
#' Return the model predictions of a \code{ppDiscSurv} object
#'
#' @param fit a \code{ppDiscSurv} object.
#'
#' @param data an `dataframe` or `list` object that contains the variables for prediction.
#'
#' @param prov.char name of provider IDs variable in `data` as a character string.
#'
#' @param Z.char names of covariates in `data` as vector of character strings.
#'
#' @param Time.char name of the observation time in `data` as a character string.
#'
#' @param lambda values of the regularization parameter lambda at which predictions are requested. For values of lambda not in the sequence of fitted models, linear interpolation is used.
#'
#' @param which indices of the penalty parameter lambda at which predictions are required. By default, all indices are returned. If lambda is specified, this will override which.
#'
#' @param which.lambda determine which lambda values are included in the output of the prediction. By default, its value is set to "all," resulting in a matrix of predicted values for each lambda presented as a list. However, if specific numeric values are provided, only the predicted matrix corresponding to those specified values will be included in the output.
#' 
#' @param type type of prediction: \code{response} provides the fitted value of each person at each time point ;
#' \code{vars} returns the indices for the non-zero coefficients; \code{nvars} returns the number of non-zero coefficients;
#'
#' @param ...
#'
#' @exportS3Method predict ppDiscSurv
#'
#' @examples
#' data(Surv_Data)
#' data <- Surv_Data$data
#' Event.char <- Surv_Data$Event.char
#' prov.char <- Surv_Data$prov.char
#' Z.char <- Surv_Data$Z.char
#' Time.char <- Surv_Data$Time.char
#' fit <- pp.DiscSurv(data, Event.char, prov.char, Z.char, Time.char)
#' predict(fit, data, Event.char, prov.char, Z.char, Time.char, lambda = fit$lambda, type = "response", which.lambda = fit$lambda[1])[[1]][1:5,]
#' predict(fit, data, Event.char, prov.char, Z.char, Time.char, lambda = 0.04, type = "vars")


predict.ppDiscSurv <- function(fit, data, Event.char, prov.char, Z.char, Time.char, lambda, which = 1:length(fit$lambda),
                               type = c("response", "vars", "nvars"), return.Array = TRUE, 
                               which.lambda = "all", ...){
  gamma <- coef.ppDiscSurv(fit, lambda = lambda, which = which, drop = FALSE)$gamma #center effect
  alpha <- coef.ppDiscSurv(fit, lambda = lambda, which = which, drop = FALSE)$alpha #time effect
  beta <- coef.ppDiscSurv(fit, lambda = lambda, which = which, drop = FALSE)$beta
  
  time.ref <- fit$time.ref
  
  if (sum(!(unique(data$time) %in% time.ref[, 2])) != 0) {
    stop("Specific time points in the new data cannot be located within the data used to fit the model, 
         rendering it impossible to determine the time effect.", call. = FALSE)
  }
  
  if (type == "vars"){
    return(drop(apply(beta != 0, 2, FUN = which)))
  }
  
  if (type == "nvars") {
    v <- as.list(apply(beta != 0, 2, FUN = which))
    nvars <- sapply(v, length)
    return(nvars)
  }
  
  # predict response
  if (missing(data) | is.null(data)) {
    stop("Must supply data for predictions", call. = FALSE)
  }
  
  # Recode time
  data <- dplyr::inner_join(data, time.ref, by = c("time"))
  data <- data[ , colnames(data) != Time.char]
  colnames(data)[ncol(data)] <- Time.char
  max.time.new_data <- max(data$time)
  
  if (type == "response") { # need expand data. The return value corresponds to the fitted value of each person at each time point 
    n.obs <- nrow(data)
    gamma.temp <- cbind(rownames(gamma), as.data.frame(gamma))
    colnames(gamma.temp)[1] <- prov.char
    gamma.individual <- as.matrix(merge(data[, prov.char, drop = F], gamma.temp, by = c(prov.char))[, 2:(length(lambda) + 1)])
    eta <- as.matrix(data[, Z.char, drop = F]) %*% beta + gamma.individual # eta = gamma + Z\beta
    sum.time <- sum(data[, Time.char, drop = F])
    time <- as.matrix(data[, Time.char, drop = F])
    nlambda <- ncol(eta)
    pred.prob <- predict_linear_predictor(nlambda, n.obs, sum.time, time, alpha, eta)  # return a matrix
    
    #convert into a 3-dim array
    expand.time_ID <- discSurv::dataLong(dataShort = data[, c(Time.char, Event.char)], timeColumn = Time.char,
                                         eventColumn = Event.char, timeAsFactor = TRUE)
    pred.prob.individual <- cbind(expand.time_ID[, 1:2], pred.prob)
    colnames(pred.prob.individual) <- c("Individual", "time", paste0("lambda = ", lambda))
    
    pred.prob.array <- list()
    for (i in 1:length(lambda)){
      pred.prob.array[[i]] <- as.data.frame(tidyr::pivot_wider(pred.prob.individual[, c(1, 2, i + 2)],
                                                               id_cols = Individual, names_from = Time.char,
                                                               values_from = paste0("lambda = ", lambda)[i]))
      colnames(pred.prob.array[[i]]) <- c("Individual", rownames(alpha)[1:max.time.new_data])
    }
    
    names(pred.prob.array) <- paste0("lambda = ", round(lambda, 5))
    
    
    if (return.Array == TRUE){
      if (which.lambda == "all") {
        return(pred.prob.array)
      } else if (!is.numeric(which.lambda)) {
        stop("which.lambda must specify a numeric variable!", call. = FALSE)
      } else {
        if (sum(!(which.lambda %in% lambda)) != 0){
          stop("which.lambda contains values that do not exist in the given lambda sequence!", call. = FALSE)
        } else {
          return(pred.prob.array[pmatch(which.lambda, lambda)])
        }
      }
    } else {
      return(pred.prob)
    }
  }
  
}