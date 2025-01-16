name_map <- c(
  "mrce_approx" = "mrce-approx",
  "capme" = "capme",
  "antac" = "antac",
  "glasso_naive" = "glasso-naive",
  "glasso_oracle" = "glasso-oracle"
)
method_order <- names(name_map)

generate_latex_table <- function(mx, sd_mx, model_num) {
  nmethod <- nrow(mx)
  latex_code <- r"(\begin{table}[ht])"
  latex_code <- paste(latex_code, r"(\centering)", sep = "\n")
  latex_code <- paste(latex_code, r"(\begin{tabular}{@{}lllll@{}})", sep = "\n")
  latex_code <- paste(latex_code, r"(\toprule)", sep = "\n")
  latex_code <- paste(latex_code,
    sprintf(r"(\multicolumn{5}{c}{Model %d}\\ \midrule)", model_num),
    sep = "\n"
  )
  latex_code <- paste(latex_code,
    r"(&\multicolumn{1}{c}{$\mathrm{tpr}_{\bOmega}$}& \multicolumn{1}{c}{$\mathrm{fpr}_{\bOmega}$} & \multicolumn{1}{c}{$\mathrm{err}_{\bOmega}$} & \multicolumn{1}{c}{$\mathrm{err}_{\B}$}\\)",
    sep = "\n"
  )
  latex_code <- paste(latex_code,
    r"(\cmidrule(lr){2-4}  \cmidrule(lr){5-5})",
    sep = "\n"
  )
  for (i in seq_len(nmethod)) {
    methodline <-
      sprintf(
        r"(\texttt{%s} & $%.3f$ ($%.3f$) & $%.3f$ ($%.3f$) & $%.3f$ ($%.3f$) & $%.3f$ ($%.3f$) \\)",
        name_map[rownames(mx)[i]], mx[i, "tpr_omega"], sd_mx[i, "tpr_omega"], mx[i, "fpr_omega"], sd_mx[i, "fpr_omega"],
        mx[i, "err_omega"], sd_mx[i, "err_omega"], mx[i, "errB"], sd_mx[i, "errB"]
      )
    methodline <- sub("$NA$ ($NA$)", r"(\multicolumn{1}{c}{--})", methodline, fixed = TRUE)
    latex_code <- paste(latex_code, methodline, sep = "\n")
  }
  latex_code <- paste(latex_code, r"(\bottomrule)", sep = "\n")
  latex_code <- paste(latex_code, r"(\end{tabular})", sep = "\n")
  latex_code <- paste(latex_code, r"(\end{table})", sep = "\n")
  latex_code
}

# setting <- "sparseB_powerOmega_highdim_results.rds"
setting <- "sparseB_powerOmega_results.rds"
results <- readRDS(file.path("results", setting))
mxB <- do.call(rbind, lapply(results$perf_B, colMeans)) |> round(3)
sd_mxB <- do.call(rbind, lapply(results$perf_B, \(mx) apply(mx, 2, sd))) |> round(3)

mxS <- do.call(rbind, lapply(results$perf_Sigma, colMeans)) |> round(3)
sd_mxS <- do.call(rbind, lapply(results$perf_Sigma, \(mx) apply(mx, 2, sd))) |> round(3)

mx <- cbind(mxB, mxS)
mx <- mx[method_order, ]
sd_mx <- cbind(sd_mxB, sd_mxS)
sd_mx <- sd_mx[method_order, ]

cat(generate_latex_table(mx, sd_mx, 5))

times <- results$times[, method_order]
meantimes <- colMeans(times)
sdtimes <- apply(times, 2, sd)
sprintf(
  r"(Model 5&$%.3f$ ($%.3f$) & $%.3f$ ($%.3f$)  & $%.3f$ ($%.3f$) & $%.3f$ ($%.3f$)  & $%.3f$ ($%.3f$) \\)",
  meantimes[1], sdtimes[1], meantimes[2], sdtimes[2], meantimes[3], sdtimes[3], meantimes[4], sdtimes[4], meantimes[5], sdtimes[5]
) |> cat()
