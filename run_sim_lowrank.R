source("methods.R")
source("performance.R")
source("cov_srrr.R")
RhpcBLASctl::blas_set_num_threads(10)

args <- commandArgs(trailingOnly = TRUE)
data_file <- args[1]
out_file <- args[2]

sim_data <- readRDS(data_file)
nrep <- length(sim_data) - 1
Sigma_E <- sim_data$setting$sigma_E
Omega_E <- solve(Sigma_E)
Sigma_X <- sim_data$setting$sigma_X

# first run to determine structure of results
pen_b_seq <- 10^seq(1, -1, length = 50)
pen_omega <- 0.2

method_funs <- list(
  "elastic" = multi_elastic_net,
  "srrr" = \(X, Y) srrr_cv(X, Y, pen_b_seq),
  "cov_srrr_approx" = \(X, Y) cov_srrr_approx_cv(X, Y, pen_omega, pen_b_seq)
)

time_method <- function(f, X, Y) {
  tictoc::tic()
  result <- f(X, Y)
  tt <- tictoc::toc(quiet = TRUE)
  result$time <- tt$toc - tt$tic
  return(result)
}

repi <- 1
d <- sim_data[[repi]]
cat("rep", repi, " ")
tictoc::tic()
fun_result <- lapply(method_funs, function(f) time_method(f, d$X, d$Y))
tictoc::toc(quiet = FALSE)

perf_B <- lapply(fun_result, function(res) unlist(performance_B(res$Bhat, d$B_star, Sigma_X, d$Xtest, d$Ytest)))
perf_Sigma <- lapply(fun_result, function(res) unlist(performance_Sigma(res, Sigma_E, Omega_E)))
times <- unlist(lapply(fun_result, function(res) res$time))
names(times) <- names(method_funs)

for (repi in seq(2, nrep)) {
  d <- sim_data[[repi]]

  cat("rep", repi, " ")
  tictoc::tic()
  fun_result <- lapply(method_funs, function(f) time_method(f, d$X, d$Y))
  tictoc::toc()

  perf_B_rep <- lapply(fun_result, function(res) unlist(performance_B(res$Bhat, d$B_star, Sigma_X, d$Xtest, d$Ytest)))
  for (i in seq_along(perf_B_rep)) {
    perf_B[[i]] <- rbind(perf_B[[i]], perf_B_rep[[i]])
  }

  perf_Sigma_rep <- lapply(fun_result, function(res) unlist(performance_Sigma(res, Sigma_E, Omega_E)))
  for (i in seq_along(perf_Sigma)) {
    perf_Sigma[[i]] <- rbind(perf_Sigma[[i]], perf_Sigma_rep[[i]])
  }

  times <- rbind(times, unlist(lapply(fun_result, function(res) res$time)))

  saveRDS(list(perf_B = perf_B, perf_Sigma = perf_Sigma, times = times), file = out_file)

  gc()
}
