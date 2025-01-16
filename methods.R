library(MRCE)
library(glmnet)
library(MSRL)
library(CVglasso)
library(dpglasso)
library(MultiVarSel)
# library(reticulate)
# use_condaenv("sklearn-env")
# sklearn <- import("sklearn")
# source_python("Noise-Covariance-Estimation/lib_using_einsum.py")

multi_elastic_net <- function(X, Y, nfold = 5) {
  n <- nrow(Y)
  folds <- c(rep(1:nfold, each = n / nfold), rep(nfold, length.out = n %% nfold))
  alpha <- seq(0.1, 1, 0.1)
  cvm_all <- matrix(NA, nrow = 10, ncol = 100)
  Bhat_list <- vector("list", length = length(alpha))
  for (ai in seq_along(alpha)) {
    fit <- cv.glmnet(X, Y,
      alpha = alpha[ai], family = "mgaussian", foldid = folds, lambda.min.ratio = 1e-4,
      standardize = FALSE, intercept = FALSE
    )
    cvm <- fit$cvm
    if (length(cvm) < 100) {
      cvm <- c(cvm, rep(Inf, 100 - length(cvm)))
    }
    cvm_all[ai, ] <- cvm
    Bhat_list[[ai]] <- Reduce(cbind, coef(fit, s = "lambda.min"))[-1, ]
  }

  min_alpha_id <- arrayInd(which.min(cvm_all), dim(cvm_all))[1]
  list(
    Bhat = Bhat_list[[min_alpha_id]], # transpose here ensures dimensions p x m
    Omega = NULL,
    cvm = cvm_all
  )
}

# L1 penalty only
multi_lasso <- function(X, Y) {
  vecy <- as.numeric(Y)
  kronx <- kronecker(Matrix::Diagonal(ncol(Y)), X)
  fit <- cv.glmnet(kronx, vecy, standardize = FALSE, intercept = FALSE, alpha = 1)
  b <- coef(fit, s = "lambda.min")[-1, ]
  list(
    Bhat = matrix(b, nrow = ncol(X), ncol = ncol(Y)),
    Omega = NULL
  )
}


mrce_rothman <- function(X, Y, pen_omega_seq, pen_b_seq, silent = TRUE, maxit = 1000) {
  mrce_fit <- mrce(X, Y,
    lam1.vec = pen_omega_seq,
    lam2.vec = pen_b_seq,
    method = "cv",
    silent = silent,
    maxit.out = maxit
  )
  cvm <- mrce_fit$cv.err
  list(
    Bhat = mrce_fit$Bhat,
    Omega = mrce_fit$omega,
    cvm = cvm,
    min_ind = arrayInd(which.min(cvm), dim(cvm))
  )
}

mrce_rothman_oracle <- function(X, Y, omega, pen_b_seq, silent = TRUE, maxit = 1000) {
  mrce_fit <- mrce(X, Y,
    omega = omega,
    lam2.vec = pen_b_seq,
    method = "fixed.omega.cv",
    silent = silent,
    maxit.out = maxit
  )
  cvm <- mrce_fit$cv.err
  list(
    Bhat = mrce_fit$Bhat,
    Omega = mrce_fit$omega,
    cvm = cvm,
    min_ind = arrayInd(which.min(cvm), dim(cvm))
  )
}

mrce_approx_cv <- function(X, Y, pen_omega_seq, pen_b_seq, nfold = 5, parallel = TRUE) {
  n <- nrow(Y)
  # Split the data into folds
  folds <- c(rep(1:nfold, each = n / nfold), rep(nfold, length.out = n %% nfold))

  # Initialize variables to store the best parameters and the minimum error
  best_rho <- NULL
  best_lam2 <- NULL

  # Initialize a matrix to store cross-validation errors
  cvm <- matrix(0, nrow = length(pen_omega_seq), ncol = length(pen_b_seq))

  # Fit the model on the training set
  b0 <- multi_lasso(X, Y)$Bhat

  for (i in seq_along(pen_omega_seq)) {
    for (j in seq_along(pen_b_seq)) {
      rho <- pen_omega_seq[i]
      lam2 <- pen_b_seq[j]

      process_fold <- function(fold) {
        # Split the data into training and validation sets
        train_idx <- which(folds != fold)
        val_idx <- which(folds == fold)

        X_train <- X[train_idx, ]
        Y_train <- Y[train_idx, ]
        X_val <- X[val_idx, ]
        Y_val <- Y[val_idx, ]
        center <- Y_train - X_train %*% b0
        sigma_hat <- t(center) %*% center / nrow(Y_train)
        g <- glasso::glasso(sigma_hat, rho = rho)
        mra <- mrce(X_train, Y_train, omega = g$wi, method = "fixed.omega", lam2 = lam2, tol.in = 1e-12)

        Y_pred <- rep(1, nrow(X_val)) %*% t(mra$muhat) + X_val %*% mra$Bhat
        fold_error <- sum((Y_val - Y_pred)^2)

        return(fold_error)
      }

      # Use mclapply to process folds in parallel
      fold_errors <- NULL
      if (parallel) {
        fold_errors <- parallel::mclapply(1:nfold, process_fold, mc.cores = nfold)
      } else {
        fold_errors <- lapply(1:nfold, process_fold)
      }

      # Average the validation error over the folds
      cv_error <- mean(unlist(fold_errors))

      # Store the cross-validation error in the matrix
      cvm[i, j] <- cv_error
    }
  }

  # Determine the minimum indices from the matrix
  min_ind <- arrayInd(which.min(cvm), dim(cvm))
  best_rho <- pen_omega_seq[min_ind[1]]
  best_lam2 <- pen_b_seq[min_ind[2]]

  # Fit the final model using the best parameters
  center <- Y - X %*% b0
  sigma_hat <- t(center) %*% center / nrow(Y)
  g <- glasso::glasso(sigma_hat, rho = best_rho)
  mra <- mrce(X, Y, omega = g$wi, method = "fixed.omega", lam2 = best_lam2, tol.in = 1e-12)

  list(
    Bhat = mra$Bhat,
    Omega = mra$omega,
    best_pen_omega = best_rho,
    best_pen_b = best_lam2,
    cvm = cvm,
    min_ind = arrayInd(which.min(cvm), dim(cvm))
  )
}

mrce_oracle_cv <- function(X, Y, omega_star, pen_b_seq, nfold = 5, parallel = TRUE) {
  n <- nrow(Y)
  # Split the data into folds
  folds <- c(rep(1:nfold, each = n / nfold), rep(nfold, length.out = n %% nfold))

  # Initialize a matrix to store cross-validation errors
  cvm <- numeric(length(pen_b_seq))

  for (j in seq_along(pen_b_seq)) {
    pen_b <- pen_b_seq[j]

    process_fold <- function(fold) {
      # Split the data into training and validation sets
      train_idx <- which(folds != fold)
      val_idx <- which(folds == fold)

      X_train <- X[train_idx, ]
      Y_train <- Y[train_idx, ]
      X_val <- X[val_idx, ]
      Y_val <- Y[val_idx, ]
      mra <- mrce(X_train, Y_train, omega = omega_star, method = "fixed.omega", lam2 = pen_b, tol.in = 1e-12)

      Y_pred <- rep(1, nrow(X_val)) %*% t(mra$muhat) + X_val %*% mra$Bhat
      fold_error <- norm(Y_val - Y_pred, type = "F")

      return(fold_error)
    }

    if (parallel) {
      fold_errors <- parallel::mclapply(1:nfold, process_fold, mc.cores = nfold)
    } else {
      fold_errors <- lapply(1:nfold, process_fold)
    }

    # Average the validation error over the folds
    cv_error <- mean(unlist(fold_errors))
    cvm[j] <- cv_error
  }

  # Determine the minimum indices from the matrix
  min_ind <- which.min(cvm)
  best_pen_b <- pen_b_seq[min_ind]

  # Fit the final model using the best parameters
  mra <- mrce(X, Y, omega = omega_star, method = "fixed.omega", lam2 = best_pen_b, tol.in = 1e-12)

  list(
    Bhat = mra$Bhat,
    Omega = mra$omega,
    best_pen_b = best_pen_b,
    cvm = cvm,
    min_ind = which.min(cvm)
  )
}

# Molstad et al. (2019)
# Goal is to estimate B
msrl <- function(X, Y, quiet = TRUE) {
  Xc <- scale(X, scale = FALSE)
  Yc <- scale(Y, scale = FALSE)
  fit <- MSRL::MSRL.cv(Xc, Yc, nlambda = 100, nfolds = 5, quiet = quiet, inner.quiet = quiet)
  b <- MSRL.coef(fit, lambda = fit$lam.min)$beta
  list(
    Bhat = b,
    Sigma = crossprod(Yc - Xc %*% b) / nrow(Y)
  )
}

tan_nce_R <- function(X, Y, Sigma_X = NULL) {
  if (is.null(Sigma_X)) {
    Sigma_X <- crossprod(X) / nrow(X)
  }
  result <- tan_nce(X, Y, Sigma_X)
  list(
    Bhat = result$Bhat,
    Sigma = result$Sigma
  )
}

mvs <- function(X, Y) {
  res <- lm(as.matrix(Y) ~ X - 1)$residuals
  S12_inv <- MultiVarSel::whitening(res, "nonparam")
  S12_inv <- tryCatch(
    MultiVarSel::whitening(res, "AR1"),
    error = function(e) MultiVarSel::whitening(res, "nonparam"),
    warning = function(w) MultiVarSel::whitening(res, "nonparam")
  )
  YY <- as.numeric(Y %*% S12_inv)
  XX <- kronecker(t(S12_inv), X)
  mvs_fit <- glmnet::cv.glmnet(XX, YY, alpha = 1, nfolds = 5, intercept = FALSE)

  b <- coef(mvs_fit, s = "lambda.min")[-1] # ignore intercept with -1
  dim(b) <- c(ncol(X), ncol(Y))
  list(
    Bhat = b,
    Omega = crossprod(S12_inv)
  )
}
