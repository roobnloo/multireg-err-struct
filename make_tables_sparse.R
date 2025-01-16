name_map <- c(
  "elastic" = "m-elnet",
  "msrl" = "msrl",
  "mvs" = "mvs",
  "mrce_approx" = "mrce-approx",
  "mrce" = "mrce",
  "mrce_oracle" = "mrce-oracle"
)
method_order <- names(name_map)

i <- 0
ar <- readRDS(file.path("results", "sparseB_arE_results.rds"))
fgn <- readRDS(file.path("results", "sparseB_fgnE_results.rds"))

ar_mxB <- do.call(rbind, lapply(ar$perf_B, colMeans)) |> round(3)
colnames(ar_mxB) <- paste(colnames(ar_mxB), "ar", sep = "_")
ar_mxSigma <- do.call(rbind, lapply(ar$perf_Sigma, colMeans)) |> round(3)
colnames(ar_mxSigma) <- paste(colnames(ar_mxSigma), "ar", sep = "_")

fgn_mxB <- do.call(rbind, lapply(fgn$perf_B, colMeans)) |> round(3)
colnames(fgn_mxB) <- paste(colnames(fgn_mxB), "fgn", sep = "_")
fgn_mxSigma <- do.call(rbind, lapply(fgn$perf_Sigma, colMeans)) |> round(3)
colnames(fgn_mxSigma) <- paste(colnames(fgn_mxSigma), "fgn", sep = "_")

mxB <- cbind(ar_mxB, fgn_mxB, ar_mxSigma, fgn_mxSigma)

ar_sd_mxB <- do.call(rbind, lapply(ar$perf_B, \(mx) apply(mx, 2, sd))) |> round(3)
colnames(ar_sd_mxB) <- paste(colnames(ar_sd_mxB), "ar", sep = "_")
fgn_sd_mxB <- do.call(rbind, lapply(fgn$perf_B, \(mx) apply(mx, 2, sd))) |> round(3)
colnames(fgn_sd_mxB) <- paste(colnames(fgn_sd_mxB), "fgn", sep = "_")

ar_sd_mxSigma <- do.call(rbind, lapply(ar$perf_Sigma, \(mx) apply(mx, 2, sd))) |> round(3)
colnames(ar_sd_mxSigma) <- paste(colnames(ar_sd_mxSigma), "ar", sep = "_")
fgn_sd_mxSigma <- do.call(rbind, lapply(fgn$perf_Sigma, \(mx) apply(mx, 2, sd))) |> round(3)
colnames(fgn_sd_mxSigma) <- paste(colnames(fgn_sd_mxSigma), "fgn", sep = "_")

sd_mxB <- cbind(ar_sd_mxB, fgn_sd_mxB, ar_sd_mxSigma, fgn_sd_mxSigma)


mxB <- mxB[method_order, ]
sd_mxB <- sd_mxB[method_order, ]

generate_latex_table <- function(mx, sd_mx, model_num) {
  nmethod <- nrow(mx)
  latex_code <- r"(\begin{table}[ht])"
  latex_code <- paste(latex_code, r"(\centering)", sep = "\n")
  latex_code <- paste(latex_code, r"(\begin{tabular}{lcccccccc})", sep = "\n")
  latex_code <- paste(latex_code, r"(\toprule)", sep = "\n")
  latex_code <- paste(latex_code,
    r"(& \multicolumn{4}{c}{Model 1} & \multicolumn{4}{c}{Model 2}\\ )",
    sep = "\n"
  )
  latex_code <- paste(latex_code, r"(\cmidrule(lr){2-5}  \cmidrule(lr){6-9})", sep = "\n")
  latex_code <- paste(latex_code,
    paste0(
      r"(&{$\mathrm{err}_{\B}$})",
      r"(&{$\mathrm{err}_{\mathrm{mod}}$})",
      r"(&{$\mathrm{err}_{\mathrm{pred}}$})",
      r"(&{$\mathrm{err}_{\bOmega}$})",
      r"(&{$\mathrm{err}_{\B}$})",
      r"(&{$\mathrm{err}_{\mathrm{mod}}$})",
      r"(&{$\mathrm{err}_{\mathrm{pred}}$})",
      r"(&{$\mathrm{err}_{\bOmega}$})",
      r"(\\ \cmidrule(lr){2-4} \cmidrule(lr){5-5}\cmidrule(lr){6-8}\cmidrule(lr){9-9})"
    ),
    sep = "\n"
  )
  for (i in seq_len(nmethod)) {
    methodline <-
      sprintf(
        r"(\texttt{%s} & $%.3f$ ($%.3f$) & $%.3f$ ($%.3f$) & $%.3f$ ($%.3f$) & $%.3f$ ($%.3f$) & $%.3f$ ($%.3f$) & $%.3f$ ($%.3f$) & $%.3f$ ($%.3f$) & $%.3f$ ($%.3f$) \\)",
        name_map[rownames(mx)[i]],
        mx[i, "errB_ar"], sd_mx[i, "errB_ar"], mx[i, "err_mod_ar"], sd_mx[i, "err_mod_ar"], mx[i, "err_pred_ar"], sd_mx[i, "err_pred_ar"], mx[i, "err_omega_ar"], sd_mx[i, "err_omega_ar"],
        mx[i, "errB_fgn"], sd_mx[i, "errB_fgn"], mx[i, "err_mod_fgn"], sd_mx[i, "err_mod_fgn"], mx[i, "err_pred_fgn"], sd_mx[i, "err_pred_fgn"], mx[i, "err_omega_fgn"], sd_mx[i, "err_omega_fgn"]
      )
    methodline <- gsub("$NA$ ($NA$)", r"(\multicolumn{1}{c}{--})", methodline, fixed = TRUE)
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

time_template <- r"(Model %d&$%.3f$ ($%.3f$) & $%.3f$ ($%.3f$) & $%.3f$ ($%.3f$) & $%.3f$ ($%.3f$) & $%.3f$ ($%.3f$) & $%.3f$ ($%.3f$)\\)"

cat(
  sprintf(
    time_template, 1,
    ar_time_mean["elastic"], ar_time_sd["elastic"],
    ar_time_mean["msrl"], ar_time_sd["msrl"],
    ar_time_mean["mvs"], ar_time_sd["mvs"],
    ar_time_mean["mrce_approx"], ar_time_sd["mrce_approx"],
    ar_time_mean["mrce"], ar_time_sd["mrce"],
    ar_time_mean["mrce_oracle"], ar_time_sd["mrce_oracle"]
  ),
  "\n"
)

cat(
  sprintf(
    time_template, 2,
    fgn_time_mean["elastic"], fgn_time_sd["elastic"],
    fgn_time_mean["msrl"], fgn_time_sd["msrl"],
    fgn_time_mean["mvs"], fgn_time_sd["mvs"],
    fgn_time_mean["mrce_approx"], fgn_time_sd["mrce_approx"],
    fgn_time_mean["mrce"], fgn_time_sd["mrce"],
    fgn_time_mean["mrce_oracle"], fgn_time_sd["mrce_oracle"]
  )
)
