source("methods.R")
source("performance.R")

in_dir <- normalizePath("sim-data/vary-rho")
dir.create("results/vary-rho", recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath("results/vary-rho")

in_dir <- normalizePath(in_dir)

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
out_dir <- normalizePath(out_dir, "")

method_funs <- list(
  "elastic" = multi_elastic_net,
  "mrce_approx" = \(X, Y) mrce_approx_cv(X, Y, 0.2, 10^seq(2, -5, length = 100)),
  "msrl" = msrl,
  "mvs" = mvs
)

rho_e_seq <- seq(0, 0.9, by = 0.1)
for (rho in rho_e_seq) {
  message("Running simulations for rho = ", rho)
  data_file <- file.path(in_dir, paste0("sparseB_arE_rho", rho, ".rds"))
  sim_data <- readRDS(data_file)
  nrep <- length(sim_data) - 1
  Sigma_E <- sim_data$setting$sigma_E
  Omega_E <- solve(Sigma_E)
  Sigma_X <- sim_data$setting$sigma_X

  # first run to determine structure of results
  repi <- 1
  d <- sim_data[[repi]]

  tictoc::tic(paste("rep", repi))
  fun_result <- lapply(method_funs, function(f) f(d$X, d$Y))
  tictoc::toc()

  perf_B <- lapply(fun_result, function(res) unlist(performance_B(res$Bhat, d$B_star, Sigma_X, d$Xtest, d$Ytest)))
  perf_Sigma <- lapply(fun_result, function(res) unlist(performance_Sigma(res, Sigma_E, Omega_E)))

  for (repi in seq(2, nrep)) {
    d <- sim_data[[repi]]

    tictoc::tic(paste("rep", repi))
    fun_result <- lapply(method_funs, function(f) f(d$X, d$Y))
    tictoc::toc()

    perf_B_rep <- lapply(fun_result, function(res) unlist(performance_B(res$Bhat, d$B_star, Sigma_X, d$Xtest, d$Ytest)))
    for (i in seq_along(perf_B_rep)) {
      perf_B[[i]] <- rbind(perf_B[[i]], perf_B_rep[[i]])
    }

    perf_Sigma_rep <- lapply(fun_result, function(res) unlist(performance_Sigma(res, Sigma_E, Omega_E)))
    for (i in seq_along(perf_Sigma)) {
      perf_Sigma[[i]] <- rbind(perf_Sigma[[i]], perf_Sigma_rep[[i]])
    }

    saveRDS(
      list(perf_B = perf_B, perf_Sigma = perf_Sigma),
      file = file.path(
        out_dir,
        sprintf("sparseB_arE_results_rho_%s.rds", rho)
      )
    )

    gc()
  }
}
