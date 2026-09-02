## Smoke tests for the four model families.
##
## These exist mainly to catch the failure mode that is otherwise invisible: the
## package installs cleanly even when its C++ entry points are unregistered
## (src/RcppExports.cpp missing), and every fit then dies at run time with
## "object '_grplasso_Z_max_grLasso' not found".  Anything that calls compiled
## code will catch that.

test_that("the C++ entry points are registered", {
  set.seed(1)
  Z <- matrix(rnorm(200), 50, 4)
  expect_type(grplasso:::Z_max_grLasso(Z, rnorm(50), c(0, 2, 4), c(1, 1)), "double")
})

test_that("pp.lasso fits a penalised logistic model with provider effects", {
  set.seed(2)
  sim <- Simulation_data_GroupLasso(
    list(m = 15, n.beta = 10, n.groups = 2, prop.NonZero.group = 0.5,
         prop.outlier = 0, rho = 0.5), prov.size.mean = 25)
  d <- sim$sim.data
  ## column POSITIONS are shuffled by the generator, so read the names off the
  ## data rather than constructing them
  zc <- colnames(d)[3:12]

  fit <- pp.lasso(data = d, Y.char = "Y", Z.char = zc,
                  prov.char = "Prov.ID", nlambda = 10)

  expect_s3_class(fit, "ppLasso")
  expect_equal(nrow(fit$beta), length(zc))
  expect_equal(ncol(fit$beta), length(fit$lambda))
  expect_equal(nrow(fit$gamma), length(unique(d$Prov.ID)))
  ## a lasso path must start at the null model and pick variables up as lambda falls
  expect_equal(sum(fit$beta[, 1] != 0), 0)
  expect_gt(sum(fit$beta[, ncol(fit$beta)] != 0), 0)
  expect_true(all(diff(fit$lambda) < 0))
})

test_that("grp.lasso respects the group structure", {
  set.seed(3)
  sim <- Simulation_data_GroupLasso(
    list(m = 15, n.beta = 10, n.groups = 5, prop.NonZero.group = 0.6,
         prop.outlier = 0, rho = 0.5), prov.size.mean = 25)
  d <- sim$sim.data
  zc <- colnames(d)[3:12]

  fit <- grp.lasso(data = d, Y.char = "Y", Z.char = zc,
                   prov.char = "Prov.ID", group = sim$group, nlambda = 10)
  expect_s3_class(fit, "gr_ppLasso")
  ## within a group, coefficients enter and leave together
  k <- ncol(fit$beta)
  nz_by_group <- tapply(fit$beta[, k] != 0, sim$group, function(v) length(unique(v)))
  expect_true(all(nz_by_group == 1))
})

test_that("Strat.cox fits a Cox model stratified by provider", {
  set.seed(4)
  ## Note: sim.cox() as shipped produces degenerate event rates (see its help),
  ## so generate here with gamma standardised by its own moments.
  m <- 15; p <- 8; per <- 25
  ps   <- rep(per, m)
  g    <- runif(m, 0, 2)
  gs   <- rep(g, ps); n <- sum(ps)
  prov <- rep(seq_len(m), ps)
  beta <- c(0.8, -0.6, rep(0, p - 2))
  rho  <- 0.5
  Zl <- lapply(seq_len(m), function(i)
    MASS::mvrnorm(per, mu = ((g[i] - 1) * rho / sqrt(1 / 3)) * rep(1, p),
                  Sigma = diag(1 - rho, p) + (rho - rho^2) * matrix(1, p, p)))
  Z  <- do.call(rbind, Zl)
  tt <- -log(runif(n)) / (gs * exp(Z %*% beta))
  cc <- pmin(runif(n, 1, 10), 3)
  del <- as.numeric(!(cc < tt))
  d <- data.frame(Prov.ID = prov, Z, status = del,
                  time = as.vector(tt * (del == 1) + cc * (del == 0)))
  colnames(d)[2:(p + 1)] <- paste0("Z", seq_len(p))

  expect_gt(mean(d$status), 0.3)   # the data must actually contain events

  fit <- Strat.cox(data = d, Event.char = "status", Z.char = paste0("Z", 1:p),
                   Time.char = "time", prov.char = "Prov.ID", nlambda = 10)
  expect_s3_class(fit, "strat_cox")
  expect_equal(nrow(fit$beta), p)
  expect_equal(sum(fit$beta[, 1] != 0), 0)
  expect_gt(sum(fit$beta[, ncol(fit$beta)] != 0), 0)
})

test_that("coef and predict dispatch and agree with the fitted path", {
  set.seed(5)
  sim <- Simulation_data_GroupLasso(
    list(m = 12, n.beta = 8, n.groups = 2, prop.NonZero.group = 0.5,
         prop.outlier = 0, rho = 0.5), prov.size.mean = 25)
  d <- sim$sim.data
  zc <- colnames(d)[3:10]
  fit <- pp.lasso(data = d, Y.char = "Y", Z.char = zc,
                  prov.char = "Prov.ID", nlambda = 8)

  ## coef() must dispatch positionally, which is what the S3 signature fix was for
  cf <- coef(fit)
  expect_true(is.list(cf) || is.matrix(cf) || is.numeric(cf))

  pr <- predict(fit, data = d, Z.char = zc, prov.char = "Prov.ID", type = "nvars")
  expect_length(pr, length(fit$lambda))
  expect_true(all(pr >= 0))
})

test_that("lasso coefficients are exactly zero, not merely small", {
  set.seed(6)
  sim <- Simulation_data_GroupLasso(
    list(m = 12, n.beta = 12, n.groups = 3, prop.NonZero.group = 0.34,
         prop.outlier = 0, rho = 0.5), prov.size.mean = 25)
  d <- sim$sim.data
  zc <- colnames(d)[3:14]
  fit <- pp.lasso(data = d, Y.char = "Y", Z.char = zc,
                  prov.char = "Prov.ID", nlambda = 15)
  b <- fit$beta[, 5]
  expect_gt(sum(b == 0), 0)                       # exact zeros exist
  expect_equal(sum(abs(b) > 0 & abs(b) < 1e-8), 0) # nothing hiding just above zero
})
