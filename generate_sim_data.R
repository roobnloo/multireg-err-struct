source("generate_data_fns.R")

set.seed(48105)
m <- 50
p <- 25
n <- 200
nrep <- 100
ntest <- 100

dir.create("sim-data/vary-rho", recursive = TRUE, showWarnings = FALSE)
dir.create("sim-data/vary-s1", recursive = TRUE, showWarnings = FALSE)

sigma_X <- ar_cov(p, 0.7)
sigma_E_ar <- generate_sigma_E(m, list(rho = 0.7), "ar")
sigma_E_fgn <- generate_sigma_E(m, list(H = 0.9), "fgn")

## Generate and save the data
sparseB_arE <- vector("list", nrep)
for (rep in seq_len(nrep)) {
  B_star <- generate_B_star(p, m, list(s0 = 0.5, s1 = 0.7), "sparse_grp")
  X <- generate_X(n, sigma_X)
  Y <- generate_Y(X, B_star, sigma_E_ar)

  Xtest <- generate_X(ntest, sigma_X)
  Ytest <- generate_Y(Xtest, B_star, sigma_E_ar)

  sparseB_arE[[rep]] <- list(X = X, Y = Y, Xtest = Xtest, Ytest = Ytest, B_star = B_star)
}
sparseB_arE$setting <- list(
  sigma_X = sigma_X, sigma_E = sigma_E_ar,
  m = m, p = p, n = n, ntest = ntest, rho_e = 0.7, s0 = 0.5, s1 = 0.7
)
saveRDS(sparseB_arE, file = "sim-data/sparseB_arE.rds")

sparseB_fgnE <- vector("list", nrep)
for (rep in seq_len(nrep)) {
  B_star <- generate_B_star(p, m, list(s0 = 0.5, s1 = 0.7), "sparse_grp")
  X <- generate_X(n, sigma_X)
  Y <- generate_Y(X, B_star, sigma_E_fgn)

  Xtest <- generate_X(ntest, sigma_X)
  Ytest <- generate_Y(Xtest, B_star, sigma_E_fgn)

  sparseB_fgnE[[rep]] <- list(X = X, Y = Y, Xtest = Xtest, Ytest = Ytest, B_star = B_star)
}
sparseB_fgnE$setting <- list(
  sigma_X = sigma_X, sigma_E = sigma_E_fgn,
  m = m, p = p, n = n, ntest = ntest, H = 0.9, s0 = 0.5, s1 = 0.7
)
saveRDS(sparseB_fgnE, file = "sim-data/sparseB_fgnE.rds")

lowrankB_arE <- vector("list", nrep)
for (rep in seq_len(nrep)) {
  B_star <- generate_B_star(p, m, list(r = 5), "low_rank")
  X <- generate_X(n, sigma_X)
  Y <- generate_Y(X, B_star, sigma_E_ar)

  Xtest <- generate_X(ntest, sigma_X)
  Ytest <- generate_Y(Xtest, B_star, sigma_E_ar)

  lowrankB_arE[[rep]] <- list(X = X, Y = Y, Xtest = Xtest, Ytest = Ytest, B_star = B_star)
}
lowrankB_arE$setting <- list(
  sigma_X = sigma_X, sigma_E = sigma_E_ar,
  m = m, p = p, n = n, ntest = ntest, rho_e = 0.7, r = 5
)
saveRDS(lowrankB_arE, file = "sim-data/lowrankB_arE.rds")

lowrankB_fgnE <- vector("list", nrep)
for (rep in seq_len(nrep)) {
  B_star <- generate_B_star(p, m, list(r = 5), "low_rank")
  X <- generate_X(n, sigma_X)
  Y <- generate_Y(X, B_star, sigma_E_fgn)

  Xtest <- generate_X(ntest, sigma_X)
  Ytest <- generate_Y(Xtest, B_star, sigma_E_fgn)

  lowrankB_fgnE[[rep]] <- list(X = X, Y = Y, Xtest = Xtest, Ytest = Ytest, B_star = B_star)
}
lowrankB_fgnE$setting <- list(
  sigma_X = sigma_X, sigma_E = sigma_E_fgn,
  m = m, p = p, n = n, ntest = ntest, H = 0.9, r = 5
)
saveRDS(lowrankB_fgnE, file = "sim-data/lowrankB_fgnE.rds")

set.seed(134251)
nrep_rho <- 50
rho_e_seq <- seq(0, 0.9, by = 0.1)
for (rho in rho_e_seq) {
  sparseB_arE_rho <- vector("list", nrep_rho)
  for (rep in seq_len(nrep_rho)) {
    B_star <- generate_B_star(p, m, list(s0 = 0.5, s1 = 0.7), "sparse_grp")
    X <- generate_X(n, sigma_X)
    sigma_E <- generate_sigma_E(m, list(rho = rho), "ar")
    Y <- generate_Y(X, B_star, sigma_E)

    Xtest <- generate_X(ntest, sigma_X)
    Ytest <- generate_Y(Xtest, B_star, sigma_E)

    sparseB_arE_rho[[rep]] <- list(X = X, Y = Y, Xtest = Xtest, Ytest = Ytest, B_star = B_star)
  }
  sparseB_arE_rho$setting <- list(
    sigma_X = sigma_X, sigma_E = sigma_E,
    m = m, p = p, n = n, ntest = ntest, rho_e = rho, s0 = 0.5, s1 = 0.7
  )
  saveRDS(sparseB_arE_rho, file = paste0("sim-data/vary-rho/sparseB_arE_rho", rho, ".rds"))
}

set.seed(14152)
s0 <- 0.2
nrep_s1 <- 50
s1_seq <- seq(0.1, 0.9, by = 0.1)
for (s1 in s1_seq) {
  sparseB_arE_s1 <- vector("list", nrep_s1)
  for (rep in seq_len(nrep_s1)) {
    B_star <- generate_B_star(p, m, list(s0 = s0, s1 = s1), "sparse_grp")
    X <- generate_X(n, sigma_X)
    Y <- generate_Y(X, B_star, sigma_E_ar)

    Xtest <- generate_X(ntest, sigma_X)
    Ytest <- generate_Y(Xtest, B_star, sigma_E_ar)

    sparseB_arE_s1[[rep]] <- list(X = X, Y = Y, Xtest = Xtest, Ytest = Ytest, B_star = B_star)
  }
  sparseB_arE_s1$setting <- list(
    sigma_X = sigma_X, sigma_E = sigma_E_ar,
    m = m, p = p, n = n, ntest = ntest, rho_e = 0.7, s0 = s0, s1 = s1
  )
  saveRDS(sparseB_arE_s1, file = paste0("sim-data/vary-s1/sparseB_arE_s1-", s1, ".rds"))
}

set.seed(5124)
s0 <- 0.2
nrep_s1 <- 50
s1_seq <- seq(0.1, 0.9, by = 0.1)
for (s1 in s1_seq) {
  sparseB_fgnE_s1 <- vector("list", nrep_s1)
  for (rep in seq_len(nrep_s1)) {
    B_star <- generate_B_star(p, m, list(s0 = s0, s1 = s1), "sparse_grp")
    X <- generate_X(n, sigma_X)
    Y <- generate_Y(X, B_star, sigma_E_fgn)

    Xtest <- generate_X(ntest, sigma_X)
    Ytest <- generate_Y(Xtest, B_star, sigma_E_fgn)

    sparseB_fgnE_s1[[rep]] <- list(X = X, Y = Y, Xtest = Xtest, Ytest = Ytest, B_star = B_star)
  }
  sparseB_fgnE_s1$setting <- list(
    sigma_X = sigma_X, sigma_E = sigma_E_fgn,
    m = m, p = p, n = n, ntest = ntest, H = 0.9, s0 = s0, s1 = s1
  )
  saveRDS(sparseB_fgnE_s1, file = paste0("sim-data/vary-s1/sparseB_fgnE_s1-", s1, ".rds"))
}
