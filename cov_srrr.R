# This file implements Cov-SRRR from Chen, Huang (2016).
# It is separated from the rest of the implementation due to its complexity.

# From Chen, K., et al (2013) Reduced rank regression via adaptive nuclear norm penalization.
# Returns the estimated rank.
adaptive_rank_selection <- function(X, Y) {
  result <- rrpack::rrr(Y, X)
  return(result$rank)
}

srrr_wrap <- function(X, Y, pen_b, nrank = NULL, A0 = NULL, B0 = NULL) {
  if (is.null(nrank)) {
    nrank <- adaptive_rank_selection(X, Y)
  }
  result <- rrpack::srrr(Y, X, nrank = nrank, modstr = list(lamA = pen_b, nlam = 1))
  list(
    A = result$V %*% result$D,
    B = result$U,
    Bhat = coef(result)
  )
}

# srrr_custom <- function(X, Y, pen_b, nrank, B0 = NULL) {
#   B_diff_max <- Inf
#   B_diff_frob <- Inf
#   A <- matrix(0, nrow = ncol(Y), ncol = nrank)
#   B <- B0
#   if (is.null(B)) {
#     B <- matrix(rnorm(length(B0)), nrow = nrow(B0), ncol = ncol(B0))
#   }

#   while (B_diff_max > 1e-8 && B_diff_frob > 1e-8) {
#     # fixed B, solve A
#     procrustes <- svd(t(Y) %*% X %*% B)
#     A <- tcrossprod(procrustes$u, procrustes$v)
#     result <- rrpack::srrr(Y, X, nrank = nrank, modstr = list(lamA = pen_b, nlam = 1))
#     B_new <- result$U
#     B_diff_max <- max(abs(B_new - B0))
#     B_diff_frob <- norm(B_new - B0, type = "F")
#     B0 <- B_new
#   }
# }

srrr_cv <- function(X, Y, pen_b_seq, nrank = NULL) {
  if (is.null(nrank)) {
    nrank <- adaptive_rank_selection(X, Y)
  }
  result <- rrpack::cv.srrr(Y, X, nrank = nrank, nfold = 5, modstr = list(lamA = pen_b_seq))
  # result$A %*% t(result$V), result$V is A in Chen, Huang (2016) and result$A is B.
  # result$U %*% c(result$D) %*% t(result$V)
  list(
    A = result$V %*% result$D,
    B = result$U,
    Bhat = result$coef,
    Omega = NULL,
    min_ind = result$minid,
    cvm = result$cv.path
  )
}

cov_srrr_approx_cv <- function(
    X, Y, pen_omega_seq, pen_b_seq,
    nrank = NULL, nfold = 5, verbose = FALSE, ncores = nfold) {
  result <- cov_srrr_cv(X, Y, pen_omega_seq, pen_b_seq, nrank, nfold, max_iter = 2, verbose = verbose, ncores = ncores)
  list(
    A = result$A,
    B = result$B,
    Bhat = result$Bhat,
    Omega = result$Omega,
    cvm = result$cvm,
    min_ind = arrayInd(which.min(result$cvm), dim(result$cvm))
  )
}

cov_srrr_cv <- function(
    X, Y, pen_omega_seq, pen_b_seq, nrank = NULL,
    nfold = 5, max_iter = 1000, verbose = FALSE, ncores = 1) {
  n <- nrow(Y)
  fold <- c(rep(1:nfold, each = n / nfold), rep(nfold, length.out = n %% nfold))
  if (is.null(nrank)) {
    nrank <- adaptive_rank_selection(X, Y)
  }
  cvm <- matrix(0, nrow = length(pen_omega_seq), ncol = length(pen_b_seq))
  for (i in seq_along(pen_omega_seq)) {
    A0 <- NULL
    B0 <- NULL
    Omega0 <- NULL
    for (j in seq_along(pen_b_seq)) {
      cv_errors <- lapply(1:nfold, function(f) {
        if (verbose) {
          message(sprintf("pen_omega: %f, pen_b: %f, fold: %d", pen_omega_seq[i], pen_b_seq[j], f))
        }
        idx <- fold == f
        result <- cov_srrr(X[!idx, ], Y[!idx, ],
          pen_omega = pen_omega_seq[i], pen_b = pen_b_seq[j], nrank = nrank,
          A0 = A0, B0 = B0, Omega0 = Omega0, max_iter = max_iter, verbose = verbose
        )
        fold_err <- norm(Y[idx, ] - X[idx, ] %*% tcrossprod(result$B, result$A), type = "F")
        list(error = fold_err, result = result)
      })

      cvm[i, j] <- mean(sapply(cv_errors, function(x) x$error))
      A0 <- cv_errors[[1]]$result$A
      B0 <- cv_errors[[1]]$result$B
      Omega0 <- cv_errors[[1]]$result$Omega
    }
  }
  min_ind <- arrayInd(which.min(cvm), dim(cvm))
  result <- cov_srrr(X, Y, pen_omega_seq[min_ind[1]], pen_b_seq[min_ind[2]], nrank = nrank)
  list(
    A = result$A,
    B = result$B,
    Bhat = tcrossprod(result$B, result$A),
    Omega = result$O,
    cvm = cvm
  )
}

cov_srrr <- function(
    X, Y, pen_omega, pen_b, nrank = NULL,
    A0 = NULL, B0 = NULL, Omega0 = NULL,
    max_iter = 1000, max_best = 1000, tol = 1e-8, verbose = FALSE) {
  if (is.null(nrank)) {
    nrank <- adaptive_rank_selection(X, Y)
    message(sprintf("Estimated rank: %d", nrank))
  }
  #--------------------------------------------------------------
  # Initialize A, B, O, and variables used to monitor convergence:
  # - Initialize O as a diagonal matrix of inverse variances
  if (is.null(Omega0)) {
    Yvariances <- apply(Y, 2, var)
    omega <- diag(Yvariances)
  } else {
    omega <- Omega0
  }
  if (is.null(A0) || is.null(B0)) {
    ab <- srrr_wrap(X, Y, pen_b, nrank)
    A <- ab$A
    B <- ab$B # Estimated coef matrix is B %*% t(A) like in Chen, Huang (2016)
  } else {
    A <- A0
    B <- B0
  }
  A_best <- A
  B_best <- B
  O_best <- omega

  # - Initialize convergence criterion values to arbitrary large numbers
  obj_best <- obj_val <- obj_diff <- 1e10
  # - Initialize iteration counters to 1
  itera <- 1 # main iteration counter
  iterb <- 1 # number of iterations with current best value of obj. func.
  #--------------------------------------------------------------
  # While not converged, iterate between estimation of O by GLASSO and (A,B) by SRRR.
  while (obj_diff > tol && itera < max_iter && iterb < max_best) {
    # 1. GLASSO to estimate O for fixed A,B.
    Sigma_R <- crossprod(Y - X %*% tcrossprod(B, A)) / nrow(Y)
    # omega <- glasso::glasso(Sigma_R, rho = pen_omega, penalize.diagonal = FALSE)$wi
    omega <- dpglasso::dpglasso(Sigma_R, X = omega, rho = pen_omega)$X

    # 2. SRRR on Ytilde=Y omega^{1/2} and X to estimate Atilde and B
    cholU_Omega <- tryCatch(
      {
        chol(omega)
      },
      error = function(e) {
        if (verbose) {
          message("Omega is not positive semi-definite. Exiting loop.")
        }
      }
    )
    if (is.null(cholU_Omega)) {
      return(list(A = A_best, B = B_best, Omega = O_best, Bhat = tcrossprod(B_best, A_best), rank = nrank))
    }
    # cholU_Omega <- chol(omega)
    Ytilde <- tcrossprod(Y, cholU_Omega)
    ss <- srrr_wrap(X, Ytilde, pen_b, nrank, A0 = A, B0 = B)
    # if (sum(abs(ss$D)) <= .Machine$double.eps) {
    #   ss$U <- matrix(0, nrow = nrow(ss$U), ncol = ncol(ss$U))
    # }
    B <- ss$B
    A <- solve(cholU_Omega, ss$A) # Atilde=O^{1/2}A, so A=O^{-1/2}Atilde
    # check whether objective function has converged
    obj_newval <- nll_cov_srrr(Y, X, A, B, omega, pen_omega, pen_b)
    obj_diff <- abs(obj_val - obj_newval) / abs(obj_val)
    obj_val <- obj_newval

    # Check whether objective function has improved
    if (obj_val < obj_best) {
      obj_best <- obj_val
      O_best <- omega
      B_best <- B
      A_best <- A
      iterb <- 1
    } else {
      iterb <- iterb + 1
    }
    if (verbose) {
      message(sprintf("Iterations: %d, Objective: %f", itera, obj_val))
    }
    itera <- itera + 1
  }
  if (itera == max_iter && max_iter > 2) {
    warning(sprintf("Maximum number of iterations reached. (pen_omega = %.3f, pen_b = %.3f)", pen_omega, pen_b))
  }
  if (verbose) {
    message(sprintf("Iterations: %d, Objective: %f", itera, obj_val))
  }
  return(list(A = A_best, B = B_best, Omega = O_best, Bhat = tcrossprod(B_best, A_best), rank = nrank))
}

nll_cov_srrr <- function(Y, X, A, B, O, l1, l2) {
  E <- Y - X %*% B %*% t(A)
  Sigma_R <- crossprod(E) / nrow(Y)
  t1 <- sum(diag(Sigma_R %*% O))
  t2 <- -determinant(O, logarithm = TRUE)$modulus
  t3 <- l1 * sum(abs(O[row(O) != col(O)]))
  t4 <- l2 * sum(sqrt(rowSums(B^2)))
  return(t1 + t2 + t3 + t4)
}
