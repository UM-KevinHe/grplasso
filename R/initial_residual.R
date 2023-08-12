# "SerBIN" for computing response residuals
SerBIN.residuals <- function(Y, Z, n.prov, gamma.prov, beta){
  fit <- SerBIN(Y, Z, n.prov, gamma.prov, beta)
  gamma.prov <- as.numeric(fit$gamma);
  beta <- as.numeric(fit$beta)
  gamma.obs <- rep(gamma.prov, n.prov)
  Pred <- as.numeric(plogis(gamma.obs+Z %*% beta))
  response.residual <- as.matrix(Y - Pred, ncol = 1)
  colnames(response.residual) <- "response residuals"
  ls <- list(beta = beta, gamma = gamma.prov, residual = response.residual)
  return(ls)
}

# ## discrete survival with provider information
pp.DiscSurv.residuals <- function(delta.obs, new.Z, time, alpha, beta.new, n.true_beta){
  fit <- NR_residuals(as.matrix(time), as.matrix(new.Z), as.matrix(delta.obs), as.matrix(alpha),
                      as.matrix(beta.new), tol = 1e-4, max_iter = 1e4)  #.cpp functions
  alpha <- as.numeric(fit$alpha)
  beta <- as.numeric(fit$beta)  #beta here includes original beta and dummy "center effect"
  eta <- as.matrix(new.Z) %*% as.matrix(beta)
  residuals <- DiscSurv_residuals(nrow(new.Z), delta.obs, time, alpha, eta)
  colnames(residuals) <- "Within person residuals"
  beta.out <- ifelse(n.true_beta != 0, 
                     beta[1:n.true_beta],
                     NA)
  ls <- list(beta = beta.out, 
             alpha = alpha, 
             gamma = beta[(1 + n.true_beta):length(beta)],
             residual = residuals)
  return(ls)
}

## discrete survival without provider information
DiscSurv.residuals <- function(delta.obs, Z, time, alpha, beta){
  fit <- NR_residuals(as.matrix(time), as.matrix(Z), as.matrix(delta.obs), as.matrix(alpha),
                      as.matrix(beta), tol = 1e-4, max_iter = 1e4)  
  alpha <- as.numeric(fit$alpha)
  beta <- as.numeric(fit$beta)  
  eta <- as.matrix(Z) %*% as.matrix(beta)
  residuals <- DiscSurv_residuals(nrow(Z), delta.obs, time, alpha, eta)
  colnames(residuals) <- "Within person residuals"
  ls <- list(beta = beta, 
             alpha = alpha, 
             residual = residuals)
  return(ls)
}