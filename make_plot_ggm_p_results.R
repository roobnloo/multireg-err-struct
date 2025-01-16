library(tidyverse)
library(patchwork)

name_map <- c(
  "mrce_approx" = "mrce-approx",
  "capme" = "capme",
  "antac" = "antac",
  "glasso_naive" = "glasso-naive",
  "glasso_oracle" = "glasso-oracle"
)
method_order <- names(name_map)
n_methods <- length(name_map)

p_seq <- seq(10, 100, by = 10)
mxErr_all <- array(NA, dim = c(50, n_methods, length(p_seq)))
mxTPR_all <- array(NA, dim = c(50, n_methods, length(p_seq)))
mxFPR_all <- array(NA, dim = c(50, n_methods, length(p_seq)))
mxF1_all <- array(NA, dim = c(50, n_methods, length(p_seq)))
for (i in seq_along(p_seq)) {
  p <- p_seq[i]
  data_file <- sprintf("results/vary-p/sparseB_powerOmega_results_p%d.rds", p)
  data_file_mrce <- sprintf("results/vary-p/sparseB_powerOmega_results_p%d_mrce.rds", p)
  results <- readRDS(data_file)
  results_mrce <- readRDS(data_file_mrce)
  mxErr_all[, , i] <- do.call(cbind, lapply(results$perf_Sigma, \(ps) ps[, "err_omega"]))
  mxTPR_all[, , i] <- do.call(cbind, lapply(results$perf_Sigma, \(ps) ps[, "tpr_omega"]))
  mxFPR_all[, , i] <- do.call(cbind, lapply(results$perf_Sigma, \(ps) ps[, "fpr_omega"]))

  mxErr_all[, 1, i] <- results_mrce$perf_Sigma[[1]][, "err_omega"]
  mxTPR_all[, 1, i] <- results_mrce$perf_Sigma[[1]][, "tpr_omega"]
  mxFPR_all[, 1, i] <- results_mrce$perf_Sigma[[1]][, "fpr_omega"]
}

dimnames(mxErr_all) <- list(NULL, names(results$perf_Sigma), p_seq)
dfErr_all <- reshape2::melt(mxErr_all, varnames = c("rep", "method", "p")) |> as_tibble()
dimnames(mxTPR_all) <- list(NULL, names(results$perf_Sigma), p_seq)
dfTPR_all <- reshape2::melt(mxTPR_all, varnames = c("rep", "method", "p")) |> as_tibble()
dimnames(mxFPR_all) <- list(NULL, names(results$perf_Sigma), p_seq)
dfFPR_all <- reshape2::melt(mxFPR_all, varnames = c("rep", "method", "p")) |> as_tibble()

ggthemr::ggthemr("pale")
dfFPR_all |>
  ggplot(mapping = aes(x = factor(p), y = value, fill = method)) +
  geom_boxplot() +
  labs(x = expression(p), y = expression(err[B])) +
  theme_minimal()

dfErr_all |>
  ggplot(mapping = aes(x = p, y = value, color = method, group = method, shape = method)) +
  geom_point(stat = "summary", fun = "mean", size = 2.5) +
  geom_line(stat = "summary", fun = "mean", mapping = aes(linetype = method)) +
  labs(x = expression(p), y = expression(err[Omega]))

tpr_plot <- dfTPR_all |>
  ggplot(mapping = aes(x = p, y = value, color = method, group = method, shape = method)) +
  geom_point(stat = "summary", fun = "mean", size = 2.5) +
  geom_line(stat = "summary", fun = "mean", mapping = aes(linetype = method)) +
  labs(x = expression(p), y = expression(tpr[Omega])) +
  ylim(0, 1) +
  theme(legend.position = "none")

fpr_plot <- dfFPR_all |>
  ggplot(mapping = aes(x = p, y = value, color = method, group = method, shape = method)) +
  geom_point(stat = "summary", fun = "mean", size = 2.5) +
  geom_line(stat = "summary", fun = "mean", mapping = aes(linetype = method)) +
  labs(x = expression(p), y = expression(fpr[Omega])) +
  ylim(0, 1)

tpr_plot + fpr_plot

# dir.create("plots", showWarnings = FALSE)
# ggsave("plots/ggm_p_plot.png", units = "px", device = "png")
