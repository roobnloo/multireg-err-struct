library(tidyverse)

name_map <- c(
  "elastic" = "m-elnet",
  "mrce_approx" = "mrce-approx",
  "cov-sr3_approx" = "covSR3-approx",
  "msrl" = "msrl",
  "tan_nce" = "nce",
  "mvs" = "mvs",
  "mrce_oracle" = "mrce-oracle"
)
method_order <- names(name_map)

n_methods <- 6

mxB_all <- array(NA, dim = c(50, 6, 10))
mxmod_all <- array(NA, dim = c(50, 6, 10))
mxB <- array(NA, dim = c(6, 5, 10))
sd_mxB <- array(NA, dim = c(6, 5, 10))
rho_seq <- seq(0, 0.9, 0.1)
for (i in seq_along(rho_seq)) {
  rho <- rho_seq[i]
  data_file <- paste0("results/sparseB_arE_results_rho_", rho, ".rds")
  results <- readRDS(data_file)
  mxB[, , i] <- do.call(rbind, lapply(results$perf_B, colMeans)) |> round(3)
  sd_mxB[, , i] <- do.call(rbind, lapply(results$perf_B, \(mx) apply(mx, 2, sd))) |> round(3)
  mxB_all[, , i] <- do.call(cbind, lapply(results$perf_B, \(pb) pb[, "errB"])) |> round(3)
}
dimnames(mxB) <- list(names(results$perf_B), colnames(results$perf_B[[1]]), rho_seq)
df <- reshape2::melt(mxB, varnames = c("method", "metric", "rho")) |> as_tibble()

dimnames(mxB_all) <- list(NULL, names(results$perf_B), rho_seq)
dfB_all <- reshape2::melt(mxB_all, varnames = c("rep", "method", "rho")) |> as_tibble()

ggthemr::ggthemr("fresh")
dfB_all |>
  mutate(method = fct_recode(method, "m-elnet" = "elastic")) |>
  filter(method != "tan_nce" & method != "mrce_oracle") |>
  ggplot(mapping = aes(x = factor(rho), y = value, fill = method)) +
  geom_boxplot() +
  labs(x = expression(rho), y = expression(err[B])) +
  theme_minimal()

# dir.create("plots", showWarnings = FALSE)
# ggsave("plots/ar_results.png", device = "png")
