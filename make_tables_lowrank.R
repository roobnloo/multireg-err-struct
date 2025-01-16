name_map <- c(
  "elastic" = "m-elnet",
  "srrr" = "srrr",
  "cov_srrr_approx" = "covSR3-approx"
)
method_order <- names(name_map)

i <- 0
ar <- readRDS(file.path("results", "lowrankB_arE_results.rds"))
fgn <- readRDS(file.path("results", "lowrankB_fgnE_results.rds"))
ar_mxB <- do.call(rbind, lapply(ar$perf_B, colMeans)) |> round(3)
colnames(ar_mxB) <- paste(colnames(ar_mxB), "ar", sep = "_")
fgn_mxB <- do.call(rbind, lapply(fgn$perf_B, colMeans)) |> round(3)
colnames(fgn_mxB) <- paste(colnames(fgn_mxB), "fgn", sep = "_")
mxB <- cbind(ar_mxB, fgn_mxB)

ar_sd_mxB <- do.call(rbind, lapply(ar$perf_B, \(mx) apply(mx, 2, sd))) |> round(3)
colnames(ar_sd_mxB) <- paste(colnames(ar_sd_mxB), "ar", sep = "_")
fgn_sd_mxB <- do.call(rbind, lapply(fgn$perf_B, \(mx) apply(mx, 2, sd))) |> round(3)
colnames(fgn_sd_mxB) <- paste(colnames(fgn_sd_mxB), "fgn", sep = "_")
sd_mxB <- cbind(ar_sd_mxB, fgn_sd_mxB)

generate_latex_table <- function(mx, sd_mx, model_num) {
  nmethod <- nrow(mx)
  latex_code <- r"(\begin{table}[ht])"
  latex_code <- paste(latex_code, r"(\centering)", sep = "\n")
  latex_code <- paste(latex_code, r"(\begin{tabular}{lllllll})", sep = "\n")
  latex_code <- paste(latex_code, r"(\toprule)", sep = "\n")
  latex_code <- paste(latex_code,
    r"(& \multicolumn{3}{c}{Model 3} & \multicolumn{3}{c}{Model 4}\\ )",
    sep = "\n"
  )
  latex_code <- paste(latex_code, r"(\cmidrule(lr){2-4}  \cmidrule(lr){5-7})", sep = "\n")
  latex_code <- paste(latex_code,
    paste0(
      r"(&\multicolumn{1}{c}{$\mathrm{err}_{\B}$})",
      r"(&\multicolumn{1}{c}{$\mathrm{err}_{\mathrm{mod}}$})",
      r"(&\multicolumn{1}{c}{$\mathrm{err}_{\mathrm{pred}}$})",
      r"(&\multicolumn{1}{c}{$\mathrm{err}_{\B}$})",
      r"(&\multicolumn{1}{c}{$\mathrm{err}_{\mathrm{mod}}$})",
      r"(&\multicolumn{1}{c}{$\mathrm{err}_{\mathrm{pred}}$})",
      r"(\\ \cmidrule(lr){2-4} \cmidrule(lr){5-7})"
    ),
    sep = "\n"
  )
  for (i in seq_len(nmethod)) {
    methodline <-
      sprintf(
        r"(\texttt{%s} & $%.3f$ ($%.3f$) & $%.3f$ ($%.3f$) & $%.3f$ ($%.3f$) & $%.3f$ ($%.3f$) & $%.3f$ ($%.3f$) & $%.3f$ ($%.3f$)  \\)",
        name_map[rownames(mx)[i]],
        mx[i, "errB_ar"], sd_mx[i, "errB_ar"], mx[i, "err_mod_ar"], sd_mx[i, "err_mod_ar"], mx[i, "err_pred_ar"], sd_mx[i, "err_pred_ar"],
        mx[i, "errB_fgn"], sd_mx[i, "errB_fgn"], mx[i, "err_mod_fgn"], sd_mx[i, "err_mod_fgn"], mx[i, "err_pred_fgn"], sd_mx[i, "err_pred_fgn"]
      )
    latex_code <- paste(latex_code, methodline, sep = "\n")
  }
  latex_code <- paste(latex_code, r"(\bottomrule)", sep = "\n")
  latex_code <- paste(latex_code, r"(\end{tabular})", sep = "\n")
  latex_code <- paste(latex_code, r"(\end{table})", sep = "\n")
  latex_code
}

cat(generate_latex_table(mxB, sd_mxB))

ar_time_mean <- colMeans(ar$times)[method_order]
ar_time_sd <- apply(ar$times, 2, sd)[method_order]
fgn_time_mean <- colMeans(fgn$times)[method_order]
fgn_time_sd <- apply(fgn$times, 2, sd)[method_order]

time_template <- r"(Model %d&$%.3f$ ($%.3f$) & $%.3f$ ($%.3f$) & $%.3f$ ($%.3f$)\\)"

cat(
  sprintf(
    time_template, 3,
    ar_time_mean["elastic"], ar_time_sd["elastic"],
    ar_time_mean["srrr"], ar_time_sd["srrr"],
    ar_time_mean["cov_srrr_approx"], ar_time_sd["cov_srrr_approx"]
  ),
  "\n"
)

cat(
  sprintf(
    time_template, 4,
    fgn_time_mean["elastic"], fgn_time_sd["elastic"],
    fgn_time_mean["srrr"], fgn_time_sd["srrr"],
    fgn_time_mean["cov_srrr_approx"], fgn_time_sd["cov_srrr_approx"]
  )
)
