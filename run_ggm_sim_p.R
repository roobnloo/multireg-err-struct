source("methods.R")
source("methods_ggm.R")
source("performance.R")


in_dir <- normalizePath("sim-data/vary-p")
dir.create("results/vary-p", recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath("results/vary-p")

p_seq <- seq(10, 100, by = 10)

for (p in p_seq) {
  data_file <- file.path(in_dir, paste0("sparseB_powerOmega_p", p, ".rds"))
  sim_data <- readRDS(data_file)
  nrep <- length(sim_data) - 1
  Omega_E <- sim_data$setting$omega
  Sigma_E <- solve(Omega_E)
  Sigma_X <- sim_data$setting$sigma_X

  pen_omega_seq <- 0.2
  pen_b_seq <- 10^seq(1, -3, length = 100)
  method_funs <- list(
    "mrce_approx" = \(X, Y, B_star) mrce_approx_cv(X, Y, pen_omega_seq, pen_b_seq),
    "capme" = \(X, Y, B_star) capme_wrap(X, Y, 2, 50),
    "glasso_oracle" = glasso_oracle,
    "glasso_naive" = \(X, Y, B_star) glasso_naive(X, Y),
    "antac" = \(X, Y, B_star) antac(X, Y)
  )

  # first run to determine structure of results
  repi <- 1
  d <- sim_data[[repi]]

  tictoc::tic(paste("rep", repi))
  fun_result <- parallel::mclapply(method_funs, function(f) f(d$X, d$Y, d$B_star), mc.cores = length(method_funs))
  tictoc::toc()

  perf_B <- lapply(fun_result, function(res) unlist(performance_B(res$Bhat, d$B_star, Sigma_X, d$Xtest, d$Ytest)))
  perf_Sigma <- lapply(fun_result, function(res) unlist(performance_Sigma(res, Sigma_E, Omega_E)))

  for (repi in seq(2, nrep)) {
    d <- sim_data[[repi]]

    tictoc::tic(paste("rep", repi))
    fun_result <- lapply(method_funs, function(f) f(d$X, d$Y, d$B_star))
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
        sprintf("sparseB_powerOmega_results_p%s.rds", p)
      )
    )

    gc()
  }
}
