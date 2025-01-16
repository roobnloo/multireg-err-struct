library(ConditionalGGM)
library(capme)

antac <- function(X, Y) {
  n <- nrow(X)
  q <- ncol(Y)
  p <- ncol(X)
  # finite sample version
  smax1 <- sqrt(n) / log(q)
  b1 <- qt(1 - (smax1 / q)^(1 + log(p) / log(q)) / 2, df = n - 1)
  lambda1 <- b1 / sqrt(n - 1 + b1^2)

  # population version
  # lambda1 <- sqrt(2 * (1 + log(p) / log(q)) / nrow(X))

  # sample version
  smax2 <- sqrt(n) / log(p)
  b2 <- qt(1 - smax2 / (2 * p), df = n - 1)
  lambda2 <- b2 / sqrt(n - 1 + b2^2)

  # population version
  # lambda2 <- sqrt(2 * log(p) / n)

  antac_res <- ConditionalGGM_twostage(Y, X, lambda_star = lambda1, lambda = lambda2)
  est_pre <- antac_res$precision
  est_pre <- est_pre + t(est_pre)
  diag(est_pre) <- diag(est_pre) / 2
  est_pre <- ChangeEta(est_pre, n, 2)
  list(
    Bhat = t(antac_res$gamma),
    Omega = est_pre
  )
}

capme_wrap <- function(X, Y, nlambda, ntau) {
  capme_fit <- cv.capme(fold = 5, loss = "likelihood", x = X, y = Y, nlambda = nlambda, ntau = ntau, tau.max = 10)
  capme_final <- capme(
    x = X, y = Y, lambda = capme_fit$lambdaopt, tau = capme_fit$tauopt,
    linsolver.Gamma = "simplex", linsolver.Omega = "simplex"
  )
  list(
    Bhat = capme_final$Gammalist[[1]],
    Omega = capme_final$Omegalist[[1]]
  )
}

glasso_oracle <- function(X, Y, B_star) {
  centered <- Y - X %*% B_star
  cvg <- CVglasso::CVglasso(X = centered, nlam = 100, lam.min.ratio = 1e-4, trace = "none")
  omega <- cvg$Omega
  list(
    Bhat = B_star,
    Omega = omega
  )
}

glasso_naive <- function(X, Y) {
  Bhat <- MASS::ginv(crossprod(X)) %*% crossprod(X, Y)
  Y0 <- Y - X %*% Bhat
  cvg <- CVglasso::CVglasso(X = Y0, nlam = 100, lam.min.ratio = 1e-4, trace = "none")
  omega <- cvg$Omega
  list(
    Bhat = Bhat,
    Omega = omega
  )
}
