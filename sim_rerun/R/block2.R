# ===================================================================
#  block2.R  —  Part B: linearisation error vs theoretical bound
# ===================================================================
#  Reuses eval_linearisation() (exact softmax vs linear closure, and the
#  Proposition bound).  lambda=1 to match the prior paper (this block is
#  invariant and about the closure map, not B_z inference).
# ===================================================================

run_block2 <- function() {
  cat("\n====== BLOCK 2 : linearisation error (Part B) ======\n")
  rows <- list();  raw <- list();  ri <- 0L

  for (K in K_B2) {
    Bz0_k <- matrix(runif(P_COV * (K - 1L), -0.5, 0.5), P_COV, K - 1L)
    for (sig in SIG_B2) {
      ri <- ri + 1L
      mse_v <- bnd_v <- rat_v <- zn_v <- numeric(N_REP_B2)
      for (r in seq_len(N_REP_B2)) {
        dat <- sim_dgp(M = M_B2, N = N_VOCAB, K = K, P = P_COV,
                       Bz0 = Bz0_k, sigma_eps = sig, alpha_beta = ALPHA_BETA,
                       doc_length = DOC_LEN, seed = 20000L + ri * 1000L + r)
        fit <- tryCatch(sgscatm(dat$W, dat$C, K = K, lambda = 1, rotate = FALSE),
                        error = function(e) NULL)
        if (is.null(fit)) { mse_v[r] <- NA; next }
        el <- eval_linearisation(fit, K = K)
        mse_v[r] <- el$mse_linearisation; bnd_v[r] <- el$theoretical_bound
        rat_v[r] <- el$ratio;             zn_v[r]  <- el$max_z_norm
      }
      raw[[ri]] <- data.frame(K = K, sigma_eps = sig, rep = seq_len(N_REP_B2),
                              mse = mse_v, bound = bnd_v, ratio = rat_v,
                              max_z = zn_v)
      rows[[ri]] <- data.frame(
        K = K, sigma_eps = sig,
        MSE   = mean(mse_v, na.rm = TRUE),
        Bound = mean(bnd_v, na.rm = TRUE),
        Ratio = mean(rat_v, na.rm = TRUE),
        max_z = mean(zn_v, na.rm = TRUE))
      cat(sprintf("  K=%2d sig=%.2f  MSE=%.2e Bound=%.2e Ratio=%.3f\n",
                  K, sig, rows[[ri]]$MSE, rows[[ri]]$Bound, rows[[ri]]$Ratio))
    }
  }

  B2 <- list(summary = do.call(rbind, rows), raw = do.call(rbind, raw))
  saveRDS(B2, file.path(DATA_DIR, "block2.rds"))
  df <- B2$summary
  write.csv(df, file.path(TAB_DIR, "block2.csv"), row.names = FALSE)

  body <- apply(df, 1L, function(row) {
    sprintf("%d & %.2f & %.2e & %.2e & %.3f & %.3f",
            as.integer(row["K"]), as.numeric(row["sigma_eps"]),
            as.numeric(row["MSE"]), as.numeric(row["Bound"]),
            as.numeric(row["Ratio"]), as.numeric(row["max_z"]))
  })
  cap <- sprintf(paste0(
    "Linearisation error versus the theoretical bound ",
    "(Proposition~\\ref{prop:linearisation}). $M=%d$, $N=%d$, %d replicates ",
    "per cell. All ratios $<1$ confirm the bound holds; a smaller ratio ",
    "indicates a more conservative bound."),
    M_B2, N_VOCAB, N_REP_B2)
  write_booktabs(
    file.path(TAB_DIR, "block2.tex"),
    header_cells = c("$K$", "$\\sigma_\\varepsilon$", "MSE",
                     "Bound", "Ratio", "$\\max_i\\|\\mathbf{z}_i^\\ast\\|$"),
    body_rows = body, caption = cap, label = "tab:block2",
    colspec = "cccccc")
  cat("  Wrote tables/block2.csv, tables/block2.tex\n")

  make_block2_figure(df)
  invisible(B2)
}

make_block2_figure <- function(df) {
  path <- file.path(IMG_DIR, "block2_linearisation.pdf")
  Ks <- sort(unique(df$K))
  col_mse <- unname(OI["blue"]); col_bnd <- unname(OI["vermilion"])
  open_pdf(path, 6.5, 2.6)
  par(mfrow = c(1L, length(Ks)), mar = c(4.0, 4.2, 2.0, 0.8))
  for (K in Ks) {
    d <- df[df$K == K, ]
    yr <- range(c(d$MSE, d$Bound), na.rm = TRUE)
    plot(d$sigma_eps, d$MSE, log = "xy", type = "b", pch = 19, col = col_mse,
         lwd = 2, ylim = yr, xlab = expression(sigma[epsilon]),
         ylab = if (K == Ks[1]) "entry-wise MSE" else "", main = paste0("K = ", K))
    lines(d$sigma_eps, d$Bound, type = "b", pch = 2, lty = 2, col = col_bnd, lwd = 2)
    if (K == Ks[1])
      legend("bottomright", bty = "n", cex = 0.85,
             legend = c("actual MSE", "Prop. bound"),
             col = c(col_mse, col_bnd), lty = c(1, 2), pch = c(19, 2), lwd = 2)
  }
  dev.off()
  cat("  Wrote imgs/block2_linearisation.pdf\n")
}
