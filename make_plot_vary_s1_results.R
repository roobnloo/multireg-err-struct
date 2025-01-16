library(tidyverse)

name_map <- c(
  "elastic" = "m-elnet",
  "mrce_approx" = "mrce-approx",
  "msrl" = "msrl",
  "mvs" = "mvs",
  "mrce_oracle" = "mrce-oracle"
)

method_order <- names(name_map)

n_methods <- length(name_map)
nreps_each <- 50
n_s1 <- 9

mxB_ar_all <- array(NA, dim = c(nreps_each, n_methods, n_s1))
mxB_fgn_all <- array(NA, dim = c(nreps_each, n_methods, n_s1))
mxmod_all <- array(NA, dim = c(nreps_each, n_methods, n_s1))
mxB_ar <- array(NA, dim = c(n_methods, 5 + 1, n_s1))
mxB_fgn <- array(NA, dim = c(n_methods, 5 + 1, n_s1))
sd_mxB <- array(NA, dim = c(n_methods, 5 + 1, n_s1))
s1_seq <- seq(0.1, 0.9, 0.1)
b_sparsities <- NULL
for (i in seq_along(s1_seq)) {
  s1 <- s1_seq[i]
  results_ar <- readRDS(sprintf("results/vary-s1/sparseB_arE_s1_%s.rds", s1))
  results_fgn <- readRDS(sprintf("results/vary-s1/sparseB_fgnE_s1_%s.rds", s1))
  mxB_ar[, , i] <- do.call(rbind, lapply(results_ar$perf_B, colMeans)) |> round(3)
  mxB_fgn[, , i] <- do.call(rbind, lapply(results_fgn$perf_B, colMeans)) |> round(3)

  sd_mxB[, , i] <- do.call(rbind, lapply(results_ar$perf_B, \(mx) apply(mx, 2, sd))) |> round(3)

  mxB_ar_all[, , i] <- do.call(cbind, lapply(results_ar$perf_B, \(pb) pb[, "errB"])) |> round(3)
  mxB_fgn_all[, , i] <- do.call(cbind, lapply(results_fgn$perf_B, \(pb) pb[, "errB"])) |> round(3)

  b_sparsities <- c(b_sparsities, do.call(c, lapply(results_ar$perf_B, \(pb) pb[, "B_sparsity"])))
}

dimnames(mxB_ar) <- list(names(results_ar$perf_B), colnames(results_ar$perf_B[[1]]), s1_seq)
df_ar <- reshape2::melt(mxB_ar, varnames = c("method", "metric", "s2")) |> as_tibble()
dimnames(mxB_fgn) <- list(names(results_fgn$perf_B), colnames(results_fgn$perf_B[[1]]), s1_seq)
df_fgn <- reshape2::melt(mxB_fgn, varnames = c("method", "metric", "s2")) |> as_tibble()

dimnames(mxB_ar_all) <- list(NULL, names(results_ar$perf_B), s1_seq)
dfB_ar_all <- reshape2::melt(mxB_ar_all, varnames = c("rep", "method", "s2")) |> as_tibble()
dfB_ar_all$B_sparsity <- b_sparsities

dimnames(mxB_fgn_all) <- list(NULL, names(results_fgn$perf_B), s1_seq)
dfB_fgn_all <- reshape2::melt(mxB_fgn_all, varnames = c("rep", "method", "s2")) |> as_tibble()
dfB_fgn_all$B_sparsity <- b_sparsities

dfB_all <- bind_rows(dfB_ar_all, dfB_fgn_all)
dfB_all$e_type <- factor(rep(c("ar", "fgn"), each = nreps_each * n_methods * n_s1))

ggthemr::ggthemr("fresh")
e_type_labels <- c("ar" = "ar(0.7) error structure", "fgn" = "fgn(0.9) error structure")
dfB_all |>
  filter(value > 0) |>
  mutate(method = fct_recode(method, "m-elnet" = "elastic", "mrce-approx" = "mrce_approx")) |>
  filter(method != "mrce_oracle") |>
  ggplot(mapping = aes(x = factor(s2), y = value, fill = method)) +
  geom_boxplot() +
  facet_wrap(~e_type, nrow = 2, labeller = labeller(e_type = e_type_labels)) +
  labs(x = expression(s[2]), y = expression(err[B])) +
  theme_minimal()

# dfB_all |>
#   filter(B_sparsity > 0) |>
#   mutate(method = fct_recode(method, "m-elnet" = "elastic")) |>
#   # filter(method != "tan_nce" & method != "mrce_oracle") |>
#   ggplot(mapping = aes(x = B_sparsity, y = value, color = method)) +
#   geom_point() +
#   labs(x = expression(s2), y = expression(err[B])) +
#   theme_minimal()

# dir.create("plots", showWarnings = FALSE)
# ggsave("plots/vary-s2-plot.png", device = "png")
