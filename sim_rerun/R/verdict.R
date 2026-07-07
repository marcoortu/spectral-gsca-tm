# ===================================================================
#  verdict.R  —  SE-validation flags, lambda verdict, output inventory
# ===================================================================

run_verdict <- function() {
  cat("\n=====================================================\n")
  cat("            BLOCK 1 RERUN VERDICT\n")
  cat("=====================================================\n")
  flag <- function(name, ok, detail) {
    cat(sprintf("  [%s] %-22s %s\n", if (isTRUE(ok)) "PASS" else "FAIL", name, detail))
    isTRUE(ok)
  }
  core <- c()
  b1p <- file.path(DATA_DIR, "block1.rds")
  if (file.exists(b1p)) {
    B1 <- readRDS(b1p); agg <- B1$agg; big <- agg[agg$M >= 1000, ]

    slope_se <- coef(lm(log(agg$mean_analytic_se) ~ log(agg$M)))[2]
    core["SE_SHRINKS"] <- flag("SE_SHRINKS",
      slope_se >= -0.65 && slope_se <= -0.35,
      sprintf("log-log slope(analytic SE) = %.3f", slope_se))

    core["SE_MATCHES_EMPIRICAL"] <- flag("SE_MATCHES_EMPIRICAL",
      all(big$se_analytic_over_empirical >= 0.85 & big$se_analytic_over_empirical <= 1.20),
      paste(sprintf("M%d:%.2f", big$M, big$se_analytic_over_empirical), collapse = " "))

    jr <- big$se_analytic_over_jack
    core["SE_MATCHES_JACKKNIFE"] <- flag("SE_MATCHES_JACKKNIFE",
      all(is.finite(jr)) && all(jr >= 0.8 & jr <= 1.25),
      paste(sprintf("M%d:%.2f", big$M, jr), collapse = " "))

    core["COVERAGE"] <- flag("COVERAGE",
      all(big$coverage95 >= 0.92 & big$coverage95 <= 0.97),
      paste(sprintf("M%d:%.3f", big$M, big$coverage95), collapse = " "))

    slope_r <- coef(lm(log(agg$RMSE) ~ log(agg$M)))[2]
    flag("RMSE_DECREASES", all(diff(agg$RMSE) < 0) && slope_r >= -0.65 && slope_r <= -0.35,
      sprintf("slope=%.3f monotone=%s", slope_r, all(diff(agg$RMSE) < 0)))

    flag("NOT_COLLAPSED", agg$ratio_to_null[which.max(agg$M)] < 0.5,
      sprintf("ratio_to_null=%.3f", agg$ratio_to_null[which.max(agg$M)]))

    wr <- B1$word_rho; relgap <- min(abs(diff(wr))) / max(wr)
    flag("DEGENERACY_GONE", B1$min_gap_cov > 0.1 && relgap > 0.02,
      sprintf("Cov(z) min gap=%.3f, word-rho min rel gap=%.3f", B1$min_gap_cov, relgap))
  } else cat("  Block 1 results missing.\n")

  overall <- length(core) == 4L && all(core)
  cat("-----------------------------------------------------\n")
  cat(sprintf("  OVERALL: %s\n", if (overall) "*** SE FORMULA VALIDATED ***"
              else "*** SE FORMULA NOT FULLY VALIDATED ***"))

  # --- lambda verdict (Part 3) -----------------------------------
  ldp <- file.path(DATA_DIR, "lambda_diag.rds")
  if (file.exists(ldp)) {
    df <- readRDS(ldp)
    b1 <- df$bias[df$setting == "lambda1"]; bA <- df$bias[df$setting == "lambdaA"]
    Ms <- df$M[df$setting == "lambda1"]
    r_lo <- b1[Ms == min(Ms)] / bA[Ms == min(Ms)]
    r_hi <- b1[Ms == max(Ms)] / bA[Ms == max(Ms)]
    inert <- r_hi < r_lo && r_hi < 1.5
    cat(sprintf("  LAMBDA VERDICT: bias(1)/bias(A) M=%d vs M=%d = %.2f , %.2f -> %s\n",
        min(Ms), max(Ms), r_lo, r_hi,
        if (inert) "PENALTY ASYMPTOTICALLY INERT (Section 4 lambda-free variance confirmed)"
        else "PENALTY NOT INERT (revise Section 4)"))
  }
  cat("=====================================================\n")

  # --- inventory --------------------------------------------------
  cat("\n--- Figure inventory (imgs/) ---\n")
  want_img <- c("block1_rmse_vs_M.pdf","block1_coverage_vs_M.pdf","block1_se_vs_M.pdf",
    "block1_qqplot.pdf","block2_linearisation.pdf","block3_mse_boxplot.pdf",
    "block3_timing.pdf","block5_disagree_vs_snr.pdf","block5_bestgain_vs_M.pdf",
    "block5_ndistinct.pdf")
  miss <- character(0)
  for (f in want_img) {
    p <- file.path(IMG_DIR, f); sz <- if (file.exists(p)) file.info(p)$size else 0
    cat(sprintf("  %-32s %s (%d bytes)  %s\n", f, if (sz>0) "OK" else "MISSING", sz,
                if (sz>0) normalizePath(p) else ""))
    if (sz == 0) miss <- c(miss, f)
  }
  cat("\n--- Table inventory (tables/) ---\n")
  want_tab <- c(paste0("block", 1:5, ".csv"), paste0("block", 1:5, ".tex"),
                "lambda_diagnostic.csv")
  for (f in want_tab) {
    p <- file.path(TAB_DIR, f); sz <- if (file.exists(p)) file.info(p)$size else 0
    cat(sprintf("  %-24s %s (%d bytes)  %s\n", f, if (sz>0) "OK" else "MISSING", sz,
                if (sz>0) normalizePath(p) else ""))
    if (sz == 0) miss <- c(miss, f)
  }
  if (length(miss)) cat("\n  !! MISSING:", paste(miss, collapse=", "), "\n")
  else cat("\n  All figures and tables rendered.\n")
  invisible(list(core = core, overall = overall))
}
