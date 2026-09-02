#' Plot the cross entropy loss from a cv.ppLasso or cv.gr_ppLasso object
#'
#' Return the plot of the cross entropy loss from a \code{cv.ppLasso} or \code{cv.gr_ppLasso} object
#'
#' @param x a \code{cv.ppLasso} object.
#'
#' @param log.x whether the horizontal axis be on the log scale.
#'
#' @param vertical.line whether draws a vertical line at the value where cross-validaton error is minimized.
#'
#' @importFrom ggplot2 ggplot geom_line geom_vline geom_point geom_errorbar aes theme element_line element_text element_blank labs scale_x_continuous
#'
#' @param ... further arguments passed to or from other methods.
#'
#' @exportS3Method plot cv.ppLasso
#'
#' @examples
#' data(BinaryData)
#' data <- BinaryData$data
#' Y.char <- BinaryData$Y.char
#' prov.char <- BinaryData$prov.char
#' Z.char <- BinaryData$Z.char
#' cv.fit.pplasso <- cv.pp.lasso(data, Y.char, Z.char, prov.char, nfolds = 10)
#' plot(cv.fit.pplasso)


#' @param col.vertical.line colour of the vertical line marking the selected lambda.
#' @param col.dot colour of the plotted cross-validation error points.
#'
#' @param ... further arguments passed to or from other methods.
#'
#' @export
plot.cv.ppLasso <- function(x, log.x = TRUE, vertical.line = TRUE, col.vertical.line = "blue",
                            col.dot = "red", ...){
  fit <- x
  CVE <- fit$cve
  CVE.upper <- fit$cve + fit$cvse
  CVE.lower <- fit$cve - fit$cvse

  if (log.x == TRUE){
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

  if (vertical.line == TRUE){
    if (log.x == TRUE){
      xintercept <- log(fit$lambda.min)
    } else {
      xintercept <- fit$lambda.min
    }
    
    cv.plot <- cv.plot +
      geom_vline(xintercept = xintercept, size = 0.5, linetype = "dashed", color = col.vertical.line)
  }

  if (log.x == TRUE){
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

#' @rdname plot.cv.ppLasso
#'
#' @param x a \code{cv.gr_ppLasso} object.
#'
#' @param log.x whether the horizontal axis be on the log scale.
#'
#' @param vertical.line whether draws a vertical line at the value where cross-validaton error is minimized.
#'
#' @importFrom ggplot2 ggplot geom_line geom_vline geom_point geom_errorbar aes theme element_line element_text element_blank labs scale_x_continuous
#'
#' @param ... further arguments passed to or from other methods.
#'
#' @exportS3Method plot cv.gr_ppLasso
#'
#' @examples
#' data(BinaryData)
#' data <- BinaryData$data
#' Y.char <- BinaryData$Y.char
#' prov.char <- BinaryData$prov.char
#' Z.char <- BinaryData$Z.char
#' group <- BinaryData$group
#' cv.fit.grplasso <- cv.grp.lasso(data, Y.char, Z.char, prov.char, group = group, nfolds = 10)
#' @param col.vertical.line colour of the vertical line marking the selected lambda.
#' @param col.dot colour of the plotted cross-validation error points.
#'
#' plot(cv.fit.grplasso)
#'

#' @param ... further arguments passed to or from other methods.
#'
#' @export
plot.cv.gr_ppLasso <- function(x, log.x = TRUE, vertical.line = TRUE,
                               col.vertical.line = "blue", col.dot = "red", ...){
  fit <- x
  CVE <- fit$cve
  CVE.upper <- fit$cve + fit$cvse
  CVE.lower <- fit$cve - fit$cvse

  if (log.x == TRUE){
    lambda <- log(fit$lambda)
  } else {
    lambda <- fit$lambda
  }

  CV.figure.df <- as.data.frame(cbind(lambda, CVE, CVE.upper, CVE.lower))

  cv.plot <- ggplot(CV.figure.df, aes(lambda, CVE)) +
    geom_line(aes(y = CVE), size = 0.05, color = "blue") +
    geom_point(size = 1, color = col.dot) +
    geom_errorbar(aes(ymin = CVE.lower, ymax = CVE.upper), width = 0.1, size = 0.1) +
    theme(panel.grid = element_blank(), panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    theme(plot.title = element_text(size = 13, face="bold", family = "serif"),
          axis.title = element_text(size = 12, family = "serif"))  +
    theme(axis.text = element_text(face = "italic", size = 11, family = "serif")) +
    scale_x_continuous(trans = scales::reverse_trans(),
                       breaks = round(seq(round(max(lambda), 0), round(min(lambda), 0), by = - 1), 1))

  if (vertical.line == TRUE){
    if (log.x == TRUE){
      xintercept <- log(fit$lambda.min)
    } else {
      xintercept <- fit$lambda.min
    }
    
    cv.plot <- cv.plot +
      geom_vline(xintercept = xintercept, size = 0.5, linetype = "dashed", color = col.vertical.line)
  }

  if (log.x == TRUE){
    cv.plot <- cv.plot +
      labs(title = "",
           x = expression(log(lambda)),
           y = "cross validation error")
  } else {
    cv.plot <- cv.plot +
      labs(title = "",
           x = expression(lambda),
           y = "cross validation error")
  }
  return(cv.plot)
}
