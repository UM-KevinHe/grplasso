#' Plot the cross entropy loss from a cv.DiscSurv object
#'
#' Return the plot of the cross entropy loss from a \code{cv.DiscSurv} object
#'
#' @param fit a \code{cv.DiscSurv} object.
#'
#' @param log.x whether the horizontal axis be on the log scale.
#'
#' @param vertical.line whether draws a vertical line at the value where cross-validaton error is minimized.
#'
#' @importFrom ggplot2 ggplot geom_line geom_vline geom_point geom_errorbar aes theme element_line element_text element_blank labs scale_x_continuous
#'
#' @exportS3Method plot cv.DiscSurv
#'
#' @examples
#' data(DiscTime)
#' data <- DiscTime$data
#' Event.char <- DiscTime$Event.char
#' Z.char <- DiscTime$Z.char
#' Time.char <- DiscTime$Time.char
#' cv.fit.DiscSurv <- cv.DiscSurv(data, Event.char, Z.char, Time.char, nfolds = 10)
#' plot(cv.fit.DiscSurv)


plot.cv.DiscSurv <- function(fit, log.x = T, vertical.line = T, col.vertical.line = "blue",
                               col.dot = "red"){
  CVE <- fit$cve
  CVE.upper <- fit$cve + fit$cvse
  CVE.lower <- fit$cve - fit$cvse
  
  if (log.x == T){
    lambda <- log(fit$lambda)
  } else {
    lambda <- fit$lambda
  }
  
  CV.figure.df <- as.data.frame(cbind(lambda, CVE, CVE.upper, CVE.lower))
  
  cv.plot <- ggplot(CV.figure.df, aes(lambda, CVE)) +
    geom_line(aes(y = CVE), linewidth = 0.05, color = "blue") +
    geom_point(size = 1, color = col.dot) +
    geom_errorbar(aes(ymin = CVE.lower, ymax = CVE.upper), width = 0.1, size = 0.1) +
    theme(panel.grid = element_blank(), panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    theme(plot.title = element_text(size = 13, face="bold", family = "serif"),
          axis.title = element_text(size = 12, family = "serif"))  +
    theme(axis.text = element_text(size = 11, family = "serif")) +
    scale_x_continuous(trans = scales::reverse_trans(),
                       breaks = round(seq(round(max(lambda), 0), round(min(lambda), 0), by = - 1), 1))
  
  if (vertical.line == T){
    if (log.x == T){
      xintercept <- log(fit$lambda.min)
    } else {
      xintercept <- fit$lambda.min
    }
    
    cv.plot <- cv.plot +
      geom_vline(xintercept = xintercept, size = 0.5, linetype = "dashed", color = col.vertical.line)
  }
  
  if (log.x == T){
    cv.plot <- cv.plot +
      labs(title = "",
           x = expression(log(lambda)),
           y = "Cross Validation Error")
  } else {
    cv.plot <- cv.plot +
      labs(title = "",
           x = expression(lambda),
           y = "Cross Validation Error")
  }
  return(cv.plot)
}