performance_B <- function(Bhat, B_star, sigma_X, X_new, Y_new) {
  metrics <- list()
  metrics$errB <- sqrt(sum((Bhat - B_star)^2))
  metrics$tprB <- sum(abs(Bhat) > .Machine$double.eps & abs(B_star) > .Machine$double.eps) /
    sum(abs(B_star) > .Machine$double.eps)
  metrics$fprB <- sum(abs(Bhat) > .Machine$double.eps & abs(B_star) <= .Machine$double.eps) /
    sum(abs(B_star) <= .Machine$double.eps)
  E <- Bhat - B_star
  metrics$err_mod <- sum(diag((crossprod(E, sigma_X %*% E))))
  metrics$err_pred <- sqrt(mean((Y_new - X_new %*% Bhat)^2))
  metrics
}

performance_Sigma <- function(fit, Sigma_star, Omega_star) {
  metrics <- list()
  if (!is.null(fit$Sigma)) {
    metrics$err_sigma <- sqrt(sum((fit$Sigma - Sigma_star)^2))
    if (!is.infinite(determinant(fit$Sigma)$modulus)) {
      Sigma_inv <- solve(fit$Sigma)
      metrics$err_omega <- sqrt(sum((Sigma_inv - Omega_star)^2))
    } else {
      metrics$err_omega <- NA
    }
  } else {
    metrics$err_sigma <- NA
    if (!is.null(fit$Omega)) {
      metrics$err_omega <- sqrt(sum((fit$Omega - Omega_star)^2))
      if (!is.infinite(determinant(fit$Omega)$modulus)) {
        Omega_inv <- solve(fit$Omega)
        metrics$err_sigma <- sqrt(sum((Omega_inv - Sigma_star)^2))
      } else {
        metrics$err_sigma <- NA
      }
      metrics$tpr_omega <- sum(abs(fit$Omega) > .Machine$double.eps & abs(Omega_star) > .Machine$double.eps) /
        sum(abs(Omega_star) > .Machine$double.eps)
      metrics$fpr_omega <- sum(abs(fit$Omega) > .Machine$double.eps & abs(Omega_star) <= .Machine$double.eps) /
        sum(abs(Omega_star) <= .Machine$double.eps)
    }
  }
  metrics
}
