source("methods.R")
source("methods_ggm.R")
source("performance.R")

args <- commandArgs(trailingOnly = TRUE)
data_file <- args[1]
out_file <- args[2]

sim_data <- readRDS(data_file)
nrep <- 100
Omega_E <- sim_data$setting$omega
Sigma_E <- solve(Omega_E)
Sigma_X <- sim_data$setting$sigma_X

pen_omega_seq <- 0.2
pen_b_seq <- 10^seq(1, -3, length = 100)
method_funs <- list(
  "mrce_approx" = \(X, Y, B_star) mrce_approx_cv(X, Y, pen_omega_seq, pen_b_seq),
  "capme" = \(X, Y, B_star) capme_wrap(X, Y, 10, 50),
  "glasso_oracle" = glasso_oracle,
  "glasso_naive" = \(X, Y, B_star) glasso_naive(X, Y),
  "antac" = \(X, Y, B_star) antac(X, Y)
)

time_method <- function(f, X, Y, B_star) {
  tictoc::tic()
  result <- f(X, Y, B_star)
  tt <- tictoc::toc(quiet = TRUE)
  result$time <- tt$toc - tt$tic
  return(result)
}

# first run to determine structure of results
repi <- 1
d <- sim_data[[repi]]

cat("rep", repi, " ")
tictoc::tic()
fun_result <- parallel::mclapply(method_funs, function(f) time_method(f, d$X, d$Y, d$B_star), mc.cores = 1)
tictoc::toc()

perf_B <- lapply(fun_result, function(res) unlist(performance_B(res$Bhat, d$B_star, Sigma_X, d$Xtest, d$Ytest)))
perf_Sigma <- lapply(fun_result, function(res) unlist(performance_Sigma(res, Sigma_E, Omega_E)))
times <- unlist(lapply(fun_result, function(res) res$time))
names(times) <- names(method_funs)

for (repi in seq(2, nrep)) {
  d <- sim_data[[repi]]

  cat("rep", repi, " ")
  tictoc::tic()
  fun_result <- lapply(method_funs, function(f) time_method(f, d$X, d$Y, d$B_star))
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
