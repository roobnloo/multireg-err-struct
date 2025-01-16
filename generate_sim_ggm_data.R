source("generate_data_fns.R")

set.seed(232260)
m <- 25
p <- 25
n <- 200
nrep <- 200
ntest <- 100
sigma_X <- ar_cov(p, 0.7)

dir.create("sim-data/vary-p", recursive = TRUE, showWarnings = FALSE)

omega_power <- generate_power_Omega(m, 1)
sigma_power <- solve(omega_power)

sparseB_powerOmega <- vector("list", nrep)
for (rep in seq_len(nrep)) {
  B_star <- generate_B_star(p, m, list(s0 = 0.5, s1 = 0.7), "sparse_grp")
  X <- generate_X(n, sigma_X)
  Y <- generate_Y(X, B_star, sigma_power)

  Xtest <- generate_X(ntest, sigma_X)
  Ytest <- generate_Y(Xtest, B_star, sigma_power)

  sparseB_powerOmega[[rep]] <- list(X = X, Y = Y, Xtest = Xtest, Ytest = Ytest, B_star = B_star)
}
sparseB_powerOmega$setting <- list(
  sigma_X = sigma_X, omega = omega_power,
  m = m, p = p, n = n, ntest = ntest, graph = "power law 1.0", s0 = 0.5, s1 = 0.7
)
saveRDS(sparseB_powerOmega, file = "sim-data/sparseB_powerOmega.rds")

set.seed(1516161)
m <- 50
p <- 25
n <- 50

# omega_power_highdim <- generate_power_Omega(m, 1)
# sigma_power_highdim <- solve(omega_power_highdim)

# sparseB_powerOmega_highdim <- vector("list", nrep)
# for (rep in seq_len(nrep)) {
#   B_star <- generate_B_star(p, m, list(s0 = 0.5, s1 = 0.7), "sparse_grp")
#   X <- generate_X(n, sigma_X)
#   Y <- generate_Y(X, B_star, sigma_power_highdim)

#   Xtest <- generate_X(ntest, sigma_X)
#   Ytest <- generate_Y(Xtest, B_star, sigma_power_highdim)

#   sparseB_powerOmega_highdim[[rep]] <- list(X = X, Y = Y, Xtest = Xtest, Ytest = Ytest, B_star = B_star)
# }
# sparseB_powerOmega_highdim$setting <- list(
#   sigma_X = sigma_X, omega = omega_power_highdim,
#   m = m, p = p, n = n, ntest = ntest, graph = "power law 1.0", s0 = 0.5, s1 = 0.7
# )
# saveRDS(sparseB_powerOmega_highdim, file = "sim-data/sparseB_powerOmega_highdim.rds")

set.seed(6336)
p_seq <- seq(10, 100, by = 10)
m <- 25
n <- 200
nrep <- 50
omega_power_varyp <- generate_power_Omega(m, 1)
sigma_power_varyp <- solve(omega_power_varyp)

sparseB_powerOmega_varyp <- vector("list", nrep)
for (p in p_seq) {
  for (rep in seq_len(nrep)) {
    B_star <- generate_B_star(p, m, list(s0 = 0.5, s1 = 0.7), "sparse_grp")
    sigma_X <- ar_cov(p, 0.7)
    X <- generate_X(n, sigma_X)
    Y <- generate_Y(X, B_star, sigma_power_varyp)

    Xtest <- generate_X(ntest, sigma_X)
    Ytest <- generate_Y(Xtest, B_star, sigma_power_varyp)

    sparseB_powerOmega_varyp[[rep]] <- list(X = X, Y = Y, Xtest = Xtest, Ytest = Ytest, B_star = B_star)
  }
  sparseB_powerOmega_varyp$setting <- list(
    sigma_X = sigma_X, omega = omega_power_varyp,
    m = m, p = p, n = n, ntest = ntest, graph = "power law 1.0", s0 = 0.5, s1 = 0.7
  )
  saveRDS(sparseB_powerOmega_varyp, file = sprintf("sim-data/vary-p/sparseB_powerOmega_p%d.rds", p))
}
