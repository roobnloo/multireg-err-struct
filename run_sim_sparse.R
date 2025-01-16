source("methods.R")
source("performance.R")
# source("cov_srrr.R")
RhpcBLASctl::blas_set_num_threads(10)

args <- commandArgs(trailingOnly = TRUE)
data_file <- args[1]
out_file <- args[2]

sim_data <- readRDS(data_file)
# nrep <- length(sim_data) - 1
offset <- 0
nrep <- 100
Sigma_E <- sim_data$setting$sigma_E
Omega_E <- solve(Sigma_E)
Sigma_X <- sim_data$setting$sigma_X

pen_b_seq <- 10^seq(1, -3, length = 100)
pen_omega <- 0.2

method_funs <- list(
  "elastic" = multi_elastic_net,
  "mrce_approx" = \(X, Y) mrce_approx_cv(X, Y, pen_omega, pen_b_seq),
  "mrce_oracle" = \(X, Y) mrce_oracle_cv(X, Y, Omega_E, pen_b_seq),
  "mrce" = \(X, Y) mrce_rothman(X, Y, pen_omega, pen_b_seq),
  "msrl" = msrl,
  "mvs" = mvs
)

time_method <- function(f, X, Y) {
  tictoc::tic()
  result <- f(X, Y)
  tt <- tictoc::toc(quiet = TRUE)
  result$time <- tt$toc - tt$tic
  return(result)
}

# first run to determine structure of results
repi <- 1 + offset
d <- sim_data[[repi]]

tictoc::tic(paste("rep", repi))
fun_result <- lapply(method_funs, function(f) time_method(f, d$X, d$Y))
tictoc::toc(quiet = FALSE)

perf_B <- lapply(fun_result, function(res) unlist(performance_B(res$Bhat, d$B_star, Sigma_X, d$Xtest, d$Ytest)))
perf_Sigma <- lapply(fun_result, function(res) unlist(performance_Sigma(res, Sigma_E, Omega_E)))
times <- unlist(lapply(fun_result, function(res) res$time))
names(times) <- names(method_funs)

for (repi in seq(2, nrep)) {
  d <- sim_data[[repi + offset]]

  tictoc::tic(paste("rep", repi + offset))
  fun_result <- lapply(method_funs, function(f) time_method(f, d$X, d$Y))
  tictoc::toc(quiet = FALSE)

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
