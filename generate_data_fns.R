library(igraph)
library(MASS)

generate_orthogonal_matrix <- function(n) {
  random_matrix <- matrix(rnorm(n^2), n, n)
  qr_decomp <- qr(random_matrix)
  Q <- qr.Q(qr_decomp)
  return(Q)
}

ar_cov <- function(p, rho) {
  abs_diff <- abs(row(matrix(1:p, p, p)) - col(matrix(1:p, p, p)))
  s <- rho^abs_diff
  return(s)
}

fgn_cov <- function(p, H) {
  abs_diff <- abs(row(matrix(1:p, p, p)) - col(matrix(1:p, p, p)))
  s <- 0.5 * ((abs_diff + 1)^(2 * H) - 2 * abs_diff^(2 * H) + (abs_diff - 1)^(2 * H))
  diag(s) <- 1
  return(s)
}

generate_sigma_E <- function(m, options, setting = c("identity", "molstad_compound", "molstad_orth", "ar", "fgn")) {
  setting <- match.arg(setting)
  switch(setting,
    "identity" = diag(m),
    "molstad_compound" = {
      s <- matrix(options$xi, m, m)
      diag(s) <- 1
      return(3 * s)
    },
    "molstad_orth" = {
      orth <- generate_orthogonal_matrix(m)
      s <- orth %*% diag(seq(1 / options$xi, 1, length.out = m)) %*% t(orth)
      return(2 * s)
    },
    "ar" = ar_cov(m, options$rho),
    "fgn" = fgn_cov(m, options$H) # H = 0.9 or 0.95
  )
}

generate_B_star <- function(
    p, m, options,
    setting = c("dense", "sparse_elt", "sparse_grp", "low_rank", "sparse_low_rank")) {
  switch(setting,
    "dense" = {
      b <- matrix(rnorm(p * m), p, m)
      return(b)
    },
    "sparse_elt" = {
      sparse <- matrix(rbinom(p * m, 1, options$s0), p, m)
      b <- matrix(rnorm(p * m), p, m)
      return(sparse * b)
    },
    "sparse_grp" = {
      W <- matrix(rnorm(p * m), p, m)
      row_indicator <- rbinom(p, 1, options$s0) # non-zero row with probability s0
      P <- t(matrix(row_indicator, m, p, byrow = TRUE))
      Q <- matrix(rbinom(p * m, 1, options$s1), p, m) # non-zero element with probability s1
      b <- W * P * Q
      return(b) # row-wise and element-wise sparsity
    },
    "low_rank" = {
      r <- options$r
      U <- matrix(rnorm(p * r), p, r)
      V <- matrix(rnorm(m * r), m, r)
      b <- U %*% t(V)
      return(b)
    },
    "sparse_low_rank" = {
      r <- options$r
      U <- matrix(rnorm(p * r), p, r)
      U[-sample(seq_len(nrow(U)), options$p0), ] <- 0
      V <- matrix(rnorm(m * r), m, r)
      b <- tcrossprod(U, V)
      return(b)
    }
  )
}

generate_power_Omega <- function(m, pwr0, lim = c(0.35, 0.5)) {
  tB <- matrix(0, m, m)
  l <- lim[1]
  u <- lim[2]

  g1 <- sample_pa(m, power = pwr0, directed = FALSE) # scale-free network
  A <- as_adjacency_matrix(g1, sparse = FALSE)
  rind <- sample(1:m, m, replace = FALSE)
  A <- A[rind, rind]
  tb <- matrix(0, m, m)
  tb[lower.tri(A) & A > 0] <- sample(
    c(
      runif(sum(A), -u, -l),
      runif(sum(A), l, u)
    ),
    sum(A) / 2,
    replace = FALSE
  )
  tB <- tb + t(tb)

  tB_temp <- matrix(0, m, m)
  for (j in 1:m) {
    tB_temp[, j] <- tB[, j] / (sum(abs(tB[, j])) * 1.5) # ensure diagonal dominance
  }

  tB <- (tB_temp + t(tB_temp)) / 2
  diag(tB) <- 1
  return(tB)
}

generate_gnp_Omega <- function(m, ve, lim = c(0.35, 0.5)) {
  tB <- matrix(0, m, m)
  l <- lim[1]
  u <- lim[2]

  g1 <- sample_gnp(m, ve, directed = FALSE) # scale-free network
  A <- as_adjacency_matrix(g1, sparse = FALSE)
  tb <- matrix(0, m, m)
  tb[lower.tri(A) & A > 0] <- sample(
    c(
      runif(sum(A), -u, -l),
      runif(sum(A), l, u)
    ),
    sum(A) / 2,
    replace = FALSE
  )
  tB <- tb + t(tb)

  tB_temp <- matrix(0, m, m)
  for (j in 1:m) {
    denom <- (sum(abs(tB[, j])) * 1.5)
    if (abs(denom) < 1e-6) {
      denom <- 1
    }
    tB_temp[, j] <- tB[, j] / denom # ensure diagonal dominance
  }

  tB <- (tB_temp + t(tB_temp)) / 2
  diag(tB) <- 1
  return(tB)
}


generate_X <- function(n, x_cov) {
  X <- mvrnorm(n, rep(0, ncol(x_cov)), x_cov)
  return(X)
}

generate_Y <- function(X, B, sigma_E) {
  Y <- X %*% B + mvrnorm(nrow(X), rep(0, ncol(sigma_E)), sigma_E)
  return(Y)
}
