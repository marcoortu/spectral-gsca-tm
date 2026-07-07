# ===================================================================
#  block5.R  —  Part C: multimodality of the exact objective
# ===================================================================
#  Built from the earlier validation CSVs (results_corollary.csv,
#  results_multistart.csv).  If they are absent, prints where it looked
#  and STOPS Part C without fabricating anything.
# ===================================================================

run_block5 <- function() {
  cat("\n====== BLOCK 5 : multimodality of exact objective (Part C) ======\n")
  search_dirs <- c(file.path(ROOT, "validation"),
                   file.path(ROOT, "replication", "output"),
                   ROOT)
  find_csv <- function(name) {
    hits <- file.path(search_dirs, name)
    hits[file.exists(hits)][1L]
  }
  cor_path <- find_csv("results_corollary.csv")
  ms_path  <- find_csv("results_multistart.csv")

  if (is.na(cor_path) || is.null(cor_path)) {
    cat("  Part C ABORTED: results_corollary.csv not found. Looked in:\n")
    cat(paste0("    - ", file.path(search_dirs, "results_corollary.csv")), sep = "\n")
    return(invisible(NULL))
  }
  cat(sprintf("  Using corollary CSV : %s\n", cor_path))
  if (!is.na(ms_path)) cat(sprintf("  Using multistart CSV: %s\n", ms_path))

  cr <- read.csv(cor_path, stringsAsFactors = FALSE)
  cr$best_gain <- cr$gain_over_floor * cr$floor_num       # per-row objective gain

  # --- per (K,M) aggregation over gap_knob ------------------------
  key <- interaction(cr$K, cr$M, drop = TRUE)
  agg <- do.call(rbind, lapply(split(cr, key), function(d) {
    data.frame(
      K = d$K[1], M = d$M[1],
      mean_disagree      = mean(d$disagree, na.rm = TRUE),
      best_gain          = median(d$best_gain, na.rm = TRUE),
      n_distinct         = median(d$n_distinct_endpoints, na.rm = TRUE),
      lambda_min_H_ratio = mean(d$lambda_min_H_ratio, na.rm = TRUE),
      snr_bulk_min       = min(d$snr_bulk, na.rm = TRUE),
      snr_bulk_max       = max(d$snr_bulk, na.rm = TRUE))
  }))
  agg <- agg[order(agg$K, agg$M), ]; rownames(agg) <- NULL
  B5 <- list(agg = agg, raw = cr)
  saveRDS(B5, file.path(DATA_DIR, "block5.rds"))
  cat("\n  Block 5 aggregate:\n"); print(round(agg, 4))
  write.csv(agg, file.path(TAB_DIR, "block5.csv"), row.names = FALSE)

  body <- apply(agg, 1L, function(row) {
    sprintf("%d & %s & %.3f & %.2e & %.1f & %.3f & [%.2f, %.2f]",
            as.integer(row["K"]), fmt_int(row["M"]),
            as.numeric(row["mean_disagree"]), as.numeric(row["best_gain"]),
            as.numeric(row["n_distinct"]), as.numeric(row["lambda_min_H_ratio"]),
            as.numeric(row["snr_bulk_min"]), as.numeric(row["snr_bulk_max"]))
  })
  cap <- paste0(
    "Multimodality of the exact GSCA objective (built from the ",
    "multi-start / corollary validation runs). Per $(K,M)$ cell, averaged ",
    "over the gap knob: mean restart disagreement, median objective gain of ",
    "the best over the floor optimum, median number of distinct endpoints, ",
    "mean $\\lambda_{\\min}(H)$ ratio, and the SNR (bulk) range covered. ",
    "Disagreement is flat across SNR, indicating the competing optima are ",
    "an intrinsic feature of the objective, not a noise artefact.")
  write_booktabs(
    file.path(TAB_DIR, "block5.tex"),
    header_cells = c("$K$", "$M$", "Mean disagree", "Best gain",
                     "$n_{\\text{distinct}}$", "$\\lambda_{\\min}(H)$ ratio",
                     "SNR$_{\\text{bulk}}$ range"),
    body_rows = body, caption = cap, label = "tab:block5", colspec = "rrccccc")
  cat("  Wrote tables/block5.csv, tables/block5.tex\n")

  make_block5_figures(B5)

  # --- caption-ready console numbers ------------------------------
  imax <- which.max(cr$snr_bulk)
  cat("\n  --- Block 5 caption-ready numbers ---\n")
  cat(sprintf("  max snr_bulk reached      : %.3f\n", cr$snr_bulk[imax]))
  cat(sprintf("  disagree at that point    : %.3f\n", cr$disagree[imax]))
  cat(sprintf("  max n_distinct_endpoints  : %.1f\n", max(cr$n_distinct_endpoints, na.rm = TRUE)))
  for (K in sort(unique(agg$K))) {
    a1 <- agg$best_gain[agg$K == K & agg$M == 1000]
    a4 <- agg$best_gain[agg$K == K & agg$M == 4000]
    if (length(a1) && length(a4))
      cat(sprintf("  best_gain shrink K=%d (M=1000->4000): %.3g -> %.3g  (x%.2f)\n",
                  K, a1, a4, a4 / a1))
  }
  invisible(B5)
}

make_block5_figures <- function(B5) {
  cr <- B5$raw; agg <- B5$agg
  Kcol <- c("3" = unname(OI["blue"]), "5" = unname(OI["vermilion"]))

  ## disagree vs snr_bulk (x log), horizontal line at 0.5 ----------
  p1 <- file.path(IMG_DIR, "block5_disagree_vs_snr.pdf")
  open_pdf(p1, 5.5, 3.8)
  par(mar = c(4.2, 4.5, 1, 1))
  plot(cr$snr_bulk, cr$disagree, log = "x", pch = 19,
       col = adjustcolor(Kcol[as.character(cr$K)], 0.75), ylim = c(0, 1),
       xlab = expression(SNR[bulk]~"(log scale)"),
       ylab = "restart disagreement")
  abline(h = 0.5, lty = 2, col = "grey40", lwd = 1.5)
  legend("topright", legend = paste0("K = ", names(Kcol)), col = Kcol,
         pch = 19, bty = "n")
  dev.off()

  ## best_gain vs M (log y) for K=3,5, ~M^{-1} reference ----------
  p2 <- file.path(IMG_DIR, "block5_bestgain_vs_M.pdf")
  open_pdf(p2, 5.5, 3.8)
  par(mar = c(4.2, 4.8, 1, 1))
  yr <- range(agg$best_gain, na.rm = TRUE)
  plot(NA, xlim = range(agg$M), ylim = yr, log = "xy",
       xlab = expression(italic(M)), ylab = "median best objective gain")
  for (K in sort(unique(agg$K))) {
    d <- agg[agg$K == K, ]; d <- d[order(d$M), ]
    lines(d$M, d$best_gain, type = "b", pch = 19, col = Kcol[as.character(K)], lwd = 2)
  }
  # M^{-1} reference anchored at the smallest-M, K=3 point
  d3 <- agg[agg$K == 3, ]; d3 <- d3[order(d3$M), ]
  if (nrow(d3) >= 1) {
    Mr <- range(agg$M); ref <- d3$best_gain[1] * (Mr[1] / Mr)
    lines(Mr, ref, lty = 3, col = "grey40", lwd = 1.5)
    text(Mr[2], ref[2], expression(italic(M)^{-1}), pos = 3, col = "grey40", cex = 0.9)
  }
  legend("bottomleft", legend = paste0("K = ", names(Kcol)), col = Kcol,
         pch = 19, lty = 1, lwd = 2, bty = "n")
  dev.off()

  ## n_distinct vs K (points, averaged over M) --------------------
  p3 <- file.path(IMG_DIR, "block5_ndistinct.pdf")
  open_pdf(p3, 4.6, 3.8)
  par(mar = c(4.2, 4.5, 1, 1))
  nd <- aggregate(n_distinct ~ K, data = agg, FUN = mean)
  plot(nd$K, nd$n_distinct, type = "b", pch = 19, col = unname(OI["green"]),
       lwd = 2, xaxt = "n", xlab = "number of topics K",
       ylab = "median distinct endpoints", ylim = c(0, max(nd$n_distinct) * 1.15))
  axis(1, at = nd$K)
  dev.off()
  cat("  Wrote imgs/block5_disagree_vs_snr.pdf, _bestgain_vs_M.pdf, _ndistinct.pdf\n")
}
