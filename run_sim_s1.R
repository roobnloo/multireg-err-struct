source("methods.R")
source("performance.R")

in_dir <- normalizePath("sim-data/vary-s1")
dir.create("results/vary-p", recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath("results/vary-s1")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}
out_dir <- normalizePath(out_dir)

pen_b_seq <- 10^seq(1, -3, length = 100)
method_funs <- list(
  "elastic" = \(X, Y, Omega_E) multi_elastic_net(X, Y),
  "mrce_approx" = \(X, Y, Omega_E) mrce_approx_cv(X, Y, 0.2, pen_b_seq),
  "mrce_oracle" = \(X, Y, Omega_E) mrce_oracle_cv(X, Y, Omega_E, pen_b_seq),
  "msrl" = \(X, Y, Omega_E) msrl(X, Y),
  "mvs" = \(X, Y, Omega_E) mvs(X, Y)
)

s1_seq <- seq(0.1, 0.9, by = 0.1)
for (s1 in 0.2) {
  message("Running simulations for s1 = ", s1)
  data_file <- file.path(in_dir, paste0("sparseB_fgnE_s1-", s1, ".rds"))
  sim_data <- readRDS(data_file)
  # nrep <- length(sim_data) - 1
  nrep <- 2
  Sigma_E <- sim_data$setting$sigma_E
  Omega_E <- solve(Sigma_E)
  Sigma_X <- sim_data$setting$sigma_X

  # first run to determine structure of results
  repi <- 12
  d <- sim_data[[repi]]

  tictoc::tic(paste("rep", repi))
  fun_result <- lapply(method_funs, function(f) f(d$X, d$Y, Omega_E))
  tictoc::toc()

  perf_B <- lapply(
    fun_result,
    function(res) {
      result <- performance_B(res$Bhat, d$B_star, Sigma_X, d$Xtest, d$Ytest)
      result$B_sparsity <- sum(abs(d$B_star) > 0)
      unlist(result)
    }
  )

  perf_Sigma <- lapply(fun_result, function(res) unlist(performance_Sigma(res, Sigma_E, Omega_E)))

  for (repi in seq(2, nrep)) {
    d <- sim_data[[repi]]

    tictoc::tic(paste("rep", repi))
    fun_result <- lapply(method_funs, function(f) f(d$X, d$Y, Omega_E))
    tictoc::toc()

    perf_B_rep <- lapply(
      fun_result,
      function(res) {
        result <- performance_B(res$Bhat, d$B_star, Sigma_X, d$Xtest, d$Ytest)
        result$B_sparsity <- sum(abs(d$B_star) > 0)
        unlist(result)
      }
    )
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
        sprintf("sparseB_arE_s1_%s.rds", s1)
      )
    )

    gc()
  }
}
