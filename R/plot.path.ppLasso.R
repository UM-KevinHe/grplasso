#' Plot regularization path of coefficients from a ppLasso or gr_ppLasso object
#'
#' Return the plot the regularization path from a \code{ppLasso} or \code{gr_ppLasso} object
#'
#' @param fit a \code{ppLasso} object.
#'
#' @param log.x whether the horizontal axis be on the log scale.
#'
#' @param label whether annotates the plot with labels.
#'
#' @importFrom ggplot2 ggplot geom_line geom_abline aes theme element_line element_text element_blank labs scale_x_continuous scale_color_manual
#'
#' @exportS3Method plot ppLasso
#'
#' @examples
#' data(GLM_Data)
#' data <- GLM_Data$data
#' Y.char <- GLM_Data$Y.char
#' prov.char <- GLM_Data$prov.char
#' Z.char <- GLM_Data$Z.char
#' fit <- pp.lasso(data, Y.char, Z.char, prov.char)
#' plot(fit, label = T)

plot.ppLasso <- function(fit, log.x = T, label = F){
  beta <- fit$beta
  if (log.x == T){
    iter.num <- rep(log(fit$lambda), each = nrow(beta))
  } else {
    iter.num <- rep(fit$lambda, each = nrow(beta))
  }
  y <- as.vector(beta)
  group <- rep(1:nrow(beta), length(fit$lambda))
  path.figure.df <- as.data.frame(t(rbind(iter.num, y, group)))

  Regularization.path <- ggplot(path.figure.df, aes(iter.num, y, group = factor(group))) +
    geom_line(aes(color = factor(group)), size = 0.5) +
    geom_abline(color = "red", linetype = 2, size = 0.5, intercept = 0, slope = 0) +
    theme(panel.grid = element_blank(), panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    theme(plot.title = element_text(size = 13, face="bold", family = "serif"),
          axis.title = element_text(size = 12, family = "serif")) +
    theme(axis.text = element_text(face = "italic", size = 11)) +
    scale_x_continuous(trans = scales::reverse_trans(), breaks = round(seq(round(max(iter.num), 0), round(min(iter.num), 0), by = - 1), 1))

  if (label == T) {
    Regularization.path <- Regularization.path +
      theme(legend.text = element_text(size = 8, face = "italic", family = "serif"),
            legend.text.align = 0, legend.title = element_blank(),
            legend.title.align = 0.5) +
      scale_color_manual(values = 1:nrow(beta), name = expression(paste(beta)),
                         labels = rownames(beta))
  } else {
    Regularization.path <- Regularization.path + theme(legend.position = "none")
  }

  if (log.x == T){
    Regularization.path <- Regularization.path +
      labs(title = "",
           x = expression(log(lambda)),
           y = expression(paste(beta, " coefficients")))
  } else {
    Regularization.path <- Regularization.path +
      labs(title = "",
           x = expression(lambda),
           y = expression(paste(beta, " coefficients")))
  }
  return(Regularization.path)
}


#' @rdname plot.ppLasso
#'
#' @param fit a \code{gr_ppLasso} object.
#'
#' @param log.x whether the horizontal axis be on the log scale.
#'
#' @param label whether annotates the plot with labels.
#'
#' @importFrom ggplot2 ggplot geom_line geom_abline aes theme element_line element_text element_blank labs scale_x_continuous scale_color_manual
#'
#' @exportS3Method plot gr_ppLasso
#'
#' @examples
#' data(GLM_Data)
#' data <- GLM_Data$data
#' Y.char <- GLM_Data$Y.char
#' prov.char <- GLM_Data$prov.char
#' Z.char <- GLM_Data$Z.char
#' group <- GLM_Data$group
#' fit <- grp.lasso(data, Y.char, Z.char, prov.char, group = group)
#' plot(fit, label = T)

plot.gr_ppLasso <- function(fit, log.x = T, label = F){
  beta <- fit$beta
  if (log.x == T){
    iter.num <- rep(log(fit$lambda), each = nrow(beta))
  } else {
    iter.num <- rep(fit$lambda, each = nrow(beta))
  }
  y <- as.vector(beta)
  group <- rep(1:nrow(beta), length(fit$lambda))
  path.figure.df <- as.data.frame(t(rbind(iter.num, y, group)))

  Regularization.path <- ggplot(path.figure.df, aes(iter.num, y, group = factor(group))) +
    geom_line(aes(color = factor(group)), linewidth = 0.5) +
    geom_abline(color = "red", linetype = 2, linewidth = 0.5, intercept = 0, slope = 0) +
    theme(panel.grid = element_blank(), panel.background = element_blank(),
          axis.line = element_line(colour = "black")) +
    theme(plot.title = element_text(size = 13, face="bold", family = "serif"),
          axis.title = element_text(size = 12, family = "serif")) +
    theme(axis.text = element_text(face = "italic", size = 11)) +
    scale_x_continuous(trans = scales::reverse_trans(), breaks = round(seq(round(max(iter.num), 0), round(min(iter.num), 0), by = - 1), 1))

  if (label == T) {
    Regularization.path <- Regularization.path +
      theme(legend.text = element_text(size = 8, face = "italic", family = "serif"),
            legend.text.align = 0, legend.title = element_blank(),
            legend.title.align = 0.5) +
      scale_color_manual(values = 1:nrow(beta), name = expression(paste(beta)),
                         labels = rownames(beta))
  } else {
    Regularization.path <- Regularization.path + theme(legend.position = "none")
  }

  if (log.x == T){
    Regularization.path <- Regularization.path +
      labs(title = "",
           x = expression(log(lambda)),
           y = expression(paste(beta, " coefficients")))
  } else {
    Regularization.path <- Regularization.path +
      labs(title = "",
           x = expression(lambda),
           y = expression(paste(beta, " coefficients")))
  }
  return(Regularization.path)
}
