# ===================================================================
#  block1.R  —  Part 2: clean-DGP inference re-run, 3 variance estimates
# ===================================================================
#  analytic SE (sgscatm_vcov, corrected) vs empirical SD (gold standard)
#  vs delete-block jackknife.  Outputs table + 4 figures, all from `B1`.
# ===================================================================

run_block1 <- function() {
  cat("\n====== BLOCK 1 : clean-DGP inference (Part 2) ======\n")
  B0 <- build_clean_B0()
  saveRDS(B0, file.path(DATA_DIR, "B0.rds"))
  saveRDS(attr(B0, "R0"), file.path(DATA_DIR, "R0.rds"))
  saveRDS(attr(B0, "d"),  file.path(DATA_DIR, "score_d.rds"))
  cat("eig(Cov(z)) =", round(attr(B0, "eig_cov_z"), 3),
      "  min gap =", round(attr(B0, "min_gap"), 3), "\n")

  # data-driven lambda_A from a probe fit (truth-free)
  probe <- sgscatm(gen_clean(2000L, B0, 1L)$W, gen_clean(2000L, B0, 1L)$C,
                   K = K_TOPICS, lambda = 1, rotate = TRUE)
  lambda_A <- lambda_A_rule(probe)
  word_rho <- sgscatm_vcov(probe)$rho            # top-(K-1) O(1) word eigenvalues
  cat(sprintf("lambda_A (data-driven, (K-1)th word eigenvalue) = %.4g\n", lambda_A))
  cat("top-(K-1) word rho:", format(word_rho, digits = 3),
      " min rel gap =", round(min(abs(diff(word_rho))) / max(word_rho), 4), "\n")
  saveRDS(lambda_A, file.path(DATA_DIR, "lambda_A.rds"))

  cat(sprintf("Fitting pilot estimand B_star (M_pilot=%d)...\n", M_PILOT))
  B_star <- pilot_Bstar(B0, lambda_A)
  saveRDS(B_star, file.path(DATA_DIR, "B_star.rds"))
  null_base <- sqrt(sum(B_star^2) / length(B_star))

  Km1 <- K_TOPICS - 1L; P <- P_COV
  cells <- list()
  for (mi in seq_along(M_B1)) {
    M <- M_B1[mi]
    cat(sprintf("\n  M = %d  (%d reps)\n", M, N_REP_B1))
    reps <- vector("list", N_REP_B1)
    for (r in seq_len(N_REP_B1)) {
      seed_r <- SEED_B1 + mi * 10000L + r
      dat <- gen_clean(M, B0, seed_r)
      one <- tryCatch({
        cf <- corrected_fit2(dat, B0, lambda = lambda_A)
        al <- aligned_se2(cf, B_star)
        halfw <- ZQ * al$se
        list(B_al = al$B_al, se = al$se,
             cover = (B_star >= al$B_al - halfw) & (B_star <= al$B_al + halfw),
             mse = mean((al$B_al - B_star)^2),
             std = (al$B_al - B_star) / al$se,
             time = cf$time_fit)
      }, error = function(e) { message("rep ", r, " failed: ", e$message); NULL })
      reps[[r]] <- one
      if (r %% 10L == 0L) cat(sprintf("    rep %d/%d\n", r, N_REP_B1))
    }

    # empirical SD (gold standard) across reps
    ok <- Filter(Negate(is.null), reps)
    Bal_arr <- simplify2array(lapply(ok, `[[`, "B_al"))       # P x Km1 x nrep
    emp_sd  <- apply(Bal_arr, c(1L, 2L), sd)

    # delete-block jackknife on JK_REPS reps
    jk_ratio <- NA_real_
    if (M %in% JK_M) {
      cat(sprintf("    delete-block jackknife (%d blocks) x %d reps...\n",
                  JK_BLOCKS, JK_REPS))
      ratios <- c()
      for (r in seq_len(min(JK_REPS, length(ok)))) {
        seed_r <- SEED_B1 + mi * 10000L + r
        dat <- gen_clean(M, B0, seed_r)
        jse <- tryCatch(loo_block_jack(dat, B_star, lambda_A),
                        error = function(e) NULL)
        if (!is.null(jse)) ratios <- c(ratios, as.vector(ok[[r]]$se / jse))
      }
      jk_ratio <- if (length(ratios)) median(ratios, na.rm = TRUE) else NA
      cat(sprintf("    median analytic/jackknife = %.3f\n", jk_ratio))
    }

    cells[[as.character(M)]] <- list(M = M, reps = reps, emp_sd = emp_sd,
                                     jk_ratio = jk_ratio)
  }

  B1 <- list(B0 = B0, B_star = B_star, lambda_A = lambda_A, word_rho = word_rho,
             min_gap_cov = attr(B0, "min_gap"),
             null_base = null_base, M_values = M_B1, cells = cells)
  saveRDS(B1, file.path(DATA_DIR, "block1.rds"))

  agg <- do.call(rbind, lapply(cells, function(cl) {
    ok <- Filter(Negate(is.null), cl$reps)
    mse <- vapply(ok, `[[`, numeric(1), "mse")
    se_arr <- simplify2array(lapply(ok, `[[`, "se"))
    rmse <- sqrt(mean(mse))
    data.frame(
      M = cl$M, n_ok = length(ok), RMSE = rmse,
      coverage95 = mean(unlist(lapply(ok, `[[`, "cover"))),
      mean_analytic_se = mean(se_arr),
      mean_empirical_sd = mean(cl$emp_sd),
      se_analytic_over_empirical = mean(se_arr) / mean(cl$emp_sd),
      se_analytic_over_jack = cl$jk_ratio,
      null_baseline_RMSE = sqrt(sum(B1$B_star^2) / length(B1$B_star)),
      ratio_to_null = rmse / sqrt(sum(B1$B_star^2) / length(B1$B_star)),
      time_s = mean(vapply(ok, `[[`, numeric(1), "time")))
  }))
  rownames(agg) <- NULL
  B1$agg <- agg
  saveRDS(B1, file.path(DATA_DIR, "block1.rds"))
  cat("\n  Block 1 aggregate:\n"); print(round(agg, 4))

  write.csv(agg, file.path(TAB_DIR, "block1.csv"), row.names = FALSE)
  body <- apply(agg, 1L, function(row)
    sprintf("%s & %.4f & %.4f & %.4f & %.3f & %.2f & %.2f",
            fmt_int(row["M"]), as.numeric(row["RMSE"]),
            as.numeric(row["mean_analytic_se"]),
            as.numeric(row["mean_empirical_sd"]),
            as.numeric(row["coverage95"]),
            as.numeric(row["se_analytic_over_empirical"]),
            as.numeric(row["time_s"])))
  cap <- sprintf(paste0(
    "Corrected inference for $\\hat{\\mathbf{B}}_z$ on the clean, ",
    "well-separated design ($P=K-1=%d$, $\\mathrm{eig}(\\mathrm{Cov}(z))=",
    "(1,.7,.5,.35)$). Analytic SE from \\texttt{sgscatm\\_vcov} (corrected ",
    "three-term influence with the missing $M^{-1/2}$) versus the ",
    "across-replicate empirical SD (gold standard). $N=%d$, %d replicates, ",
    "$\\lambda=\\lambda_A$ (data-driven)."),
    P_COV, N_VOCAB, N_REP_B1)
  if (REDUCED) cap <- paste(cap, "\\emph{Reduced replicate count (runtime).}")
  write_booktabs(file.path(TAB_DIR, "block1.tex"),
    header_cells = c("$M$", "RMSE", "Analytic SE", "Empirical SD",
                     "Cover.\\ (95\\%)", "SE/SD", "Time (s)"),
    body_rows = body, caption = cap, label = "tab:block1", colspec = "rcccccc")
  cat("  Wrote tables/block1.csv, tables/block1.tex\n")

  make_block1_figures(B1)
  invisible(B1)
}

make_block1_figures <- function(B1) {
  agg <- B1$agg; M <- agg$M
  cB <- unname(OI["blue"]); cO <- unname(OI["orange"]); cG <- unname(OI["green"])

  ## RMSE vs M (log-log, dashed M^{-1/2}) --------------------------
  open_pdf(file.path(IMG_DIR, "block1_rmse_vs_M.pdf"), 5.2, 4.0)
  par(mar = c(4.2, 4.6, 1, 1))
  ref <- agg$RMSE[1] * sqrt(M[1] / M); slope <- coef(lm(log(agg$RMSE) ~ log(M)))[2]
  plot(M, agg$RMSE, log = "xy", type = "b", pch = 19, col = cB, lwd = 2,
       xlab = expression(italic(M)), ylab = expression("RMSE of"~hat(bold(B))[z]),
       ylim = range(agg$RMSE, ref) * c(0.9, 1.1))
  lines(M, ref, lty = 2, col = "grey40", lwd = 1.5)
  legend("topright", bty = "n", legend = c("RMSE", expression(italic(M)^{-1/2})),
         col = c(cB, "grey40"), lty = c(1, 2), pch = c(19, NA), lwd = 2)
  text(M[2], ref[2], bquote("slope" == .(round(slope, 3))), pos = 4, col = cB, cex = .9)
  dev.off()

  ## coverage vs M (dashed 0.95) -----------------------------------
  open_pdf(file.path(IMG_DIR, "block1_coverage_vs_M.pdf"), 5.2, 4.0)
  par(mar = c(4.2, 4.6, 1, 1))
  plot(M, agg$coverage95, log = "x", type = "b", pch = 19, col = cG, lwd = 2,
       ylim = c(0.85, 1.0), xlab = expression(italic(M)),
       ylab = "empirical coverage (nominal 95%)")
  abline(h = 0.95, lty = 2, col = "grey40", lwd = 1.5)
  dev.off()

  ## NEW: analytic SE and empirical SD vs M (log-log, must overlap)-
  open_pdf(file.path(IMG_DIR, "block1_se_vs_M.pdf"), 5.2, 4.0)
  par(mar = c(4.2, 4.6, 1, 1))
  yr <- range(agg$mean_analytic_se, agg$mean_empirical_sd)
  refm <- agg$mean_empirical_sd[1] * sqrt(M[1] / M)
  plot(M, agg$mean_analytic_se, log = "xy", type = "b", pch = 19, col = cB,
       lwd = 2, ylim = yr * c(0.85, 1.15), xlab = expression(italic(M)),
       ylab = "standard error scale")
  lines(M, agg$mean_empirical_sd, type = "b", pch = 17, col = cO, lwd = 2)
  lines(M, refm, lty = 2, col = "grey40", lwd = 1.5)
  legend("topright", bty = "n",
         legend = c("analytic SE", "empirical SD", expression(italic(M)^{-1/2})),
         col = c(cB, cO, "grey40"), lty = c(1, 1, 2), pch = c(19, 17, NA), lwd = 2)
  dev.off()

  ## QQ of standardized bias at M=2000 -----------------------------
  open_pdf(file.path(IMG_DIR, "block1_qqplot.pdf"), 4.6, 4.6)
  par(mar = c(4.2, 4.6, 1, 1))
  qq_key <- if (!is.null(B1$cells[["2000"]])) "2000" else
    tail(names(B1$cells), 1L)
  std <- unlist(lapply(Filter(Negate(is.null), B1$cells[[qq_key]]$reps),
                       function(x) as.vector(x$std)))
  std <- std[is.finite(std)]
  qqnorm(std, pch = 16, col = adjustcolor(cB, .4), cex = .6, main = "",
         xlab = "theoretical N(0,1) quantiles",
         ylab = expression("standardized bias"~(hat(B)[z]-B["*"])/SE))
  abline(0, 1, lty = 2, col = "grey40", lwd = 1.5)
  dev.off()
  cat("  Wrote imgs/block1_{rmse_vs_M,coverage_vs_M,se_vs_M,qqplot}.pdf\n")
}
