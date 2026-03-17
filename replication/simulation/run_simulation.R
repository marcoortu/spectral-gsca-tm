#!/usr/bin/env Rscript
# ===================================================================
#  ILR-Spectral-GSCA Simulation Study — Biometrika Paper
# ===================================================================
#
#  Three blocks, each tied to a specific theoretical result:
#    Block 1: Consistency & asymptotic normality of Bz (Thm 12 & 14)
#    Block 2: Linearisation error bound (Prop 15)
#    Block 3: Structural comparison with STM (empirical positioning)
#
#  Output: /output/tables/  — LaTeX tables (.tex)
#          /output/figures/ — PDF plots (.pdf)
#          /output/data/    — raw results (.rds)
#
#  Usage:
#    Rscript run_simulation.R              # run all blocks
#    Rscript run_simulation.R --block 1    # run single block
#    Rscript run_simulation.R --quick      # reduced design for testing
# ===================================================================

# --- Parse arguments -----------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
RUN_BLOCK <- if ("--block" %in% args) {
  as.integer(args[which(args == "--block") + 1L])
} else {
  1:4
}
QUICK <- "--quick" %in% args

cat("=== ILR-Spectral-GSCA Simulation Study ===\n")
cat("Blocks:", paste(RUN_BLOCK, collapse = ", "), "\n")
cat("Mode:", if (QUICK) "QUICK (reduced design)\n" else "FULL\n")

# --- Setup ---------------------------------------------------------
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# Source sgscatm package (adjust path as needed)
# devtools::load_all(".")
source("R/sgscatm_fit.R")
source("R/ilr_contrast.R")
source("R/ilr_se.R")
source("R/methods.R")
source("R/utils.R")

# Source simulation infrastructure
source("scripts/simulation/sim_dgp.R")
source("scripts/simulation/sim_utils.R")

# Output directories
dir.create("output/tables",  recursive = TRUE, showWarnings = FALSE)
dir.create("output/figures", recursive = TRUE, showWarnings = FALSE)
dir.create("output/data",    recursive = TRUE, showWarnings = FALSE)

# --- Global design parameters -------------------------------------
if (QUICK) {
  N_REP       <- 50L
  M_VALUES_B1 <- c(500L, 1000L, 2000L)
  M_VALUES_B3 <- c(1000L, 5000L)
  M_VALUES_B4 <- c(500L, 2000L)
  N_REP_B4    <- 20L
  B_BOOT      <- 100L
} else {
  N_REP       <- 500L
  M_VALUES_B1 <- c(500L, 1000L, 2000L, 5000L, 10000L)
  M_VALUES_B3 <- c(1000L, 5000L, 20000L)
  M_VALUES_B4 <- c(500L, 1000L, 5000L, 20000L)
  N_REP_B4    <- 50L
  B_BOOT      <- 200L
}

# Block 4 scenario definitions (mirrored in the worker script)
# high_sep : sparse/distinct topics (alpha_beta=0.01), moderate signal
# low_sep  : diffuse/overlapping topics (alpha_beta=1.0), moderate signal
# weak_sig : medium topics, weak covariate signal (small b_max)
SCENARIOS_B4 <- c("high_sep", "low_sep", "weak_sig")

# Fixed across all blocks
N_VOCAB   <- 500L
K_TOPICS  <- 5L
P_COV     <- 3L
CONF      <- 0.95

# True Bz0 — fixed for reproducibility
set.seed(2026)
Bz0_TRUE <- matrix(c(
   0.40, -0.20,  0.10,  0.30,   # covariate 1 → 4 ILR components
  -0.15,  0.35, -0.25,  0.05,   # covariate 2
   0.20,  0.10,  0.40, -0.30    # covariate 3
), nrow = P_COV, ncol = K_TOPICS - 1L, byrow = TRUE)

# ggplot theme for paper
theme_paper <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )
theme_set(theme_paper)


# ===================================================================
#  BLOCK 1: Consistency & Asymptotic Normality of Bz
# ===================================================================

if (1L %in% RUN_BLOCK) {
  cat("\n====== BLOCK 1: Consistency & normality ======\n")

  results_b1 <- list()

  for (M in M_VALUES_B1) {
    cat(sprintf("  M = %d, %d replicates...\n", M, N_REP))

    mse_vec      <- numeric(N_REP)
    bias_mat     <- matrix(NA, N_REP, P_COV * (K_TOPICS - 1L))
    cov_anal_mat <- matrix(NA, N_REP, P_COV * (K_TOPICS - 1L))
    time_vec     <- numeric(N_REP)

    for (r in seq_len(N_REP)) {
      if (r %% 100L == 0L) cat(sprintf("    rep %d/%d\n", r, N_REP))

      dat <- sim_dgp(M = M, N = N_VOCAB, K = K_TOPICS, P = P_COV,
                     Bz0 = Bz0_TRUE, sigma_eps = 0.3,
                     alpha_beta = 0.1, doc_length = 200L,
                     seed = 10000L * which(M_VALUES_B1 == M) + r)

      res <- tryCatch({
        t0  <- proc.time()
        fit <- sgscatm(dat$W, dat$C, K = K_TOPICS, lambda = 1,
                       rotate = TRUE)
        tf  <- (proc.time() - t0)[3]

        # Procrustes alignment
        pa <- procrustes_align(fit$Bz, dat$Bz0)

        # Analytical SE
        se_res <- tryCatch(ilr_se_analytical(fit), error = function(e) NULL)

        list(mse = pa$mse, bias = as.vector(pa$Bz_aligned - dat$Bz0),
             time = tf, se_res = se_res, R = pa$R, Bz_al = pa$Bz_aligned)
      }, error = function(e) NULL)

      if (is.null(res)) next
      mse_vec[r]     <- res$mse
      bias_mat[r, ]  <- res$bias
      time_vec[r]    <- res$time

      # Coverage of analytical CI
      if (!is.null(res$se_res)) {
        se_rot <- .rotate_se(res$se_res, res$R, P_COV, K_TOPICS - 1L)
        z_q    <- qnorm(1 - (1 - CONF) / 2)
        ci_lo  <- res$Bz_al - z_q * se_rot
        ci_hi  <- res$Bz_al + z_q * se_rot
        covers <- (dat$Bz0 >= ci_lo) & (dat$Bz0 <= ci_hi)
        cov_anal_mat[r, ] <- as.vector(covers)
      }
    }

    results_b1[[as.character(M)]] <- list(
      M = M, mse = mse_vec, bias = bias_mat,
      coverage = cov_anal_mat, time = time_vec
    )
  }

  saveRDS(results_b1, "output/data/block1_results.rds")

  # --- Block 1 Summary Table (LaTeX) -------------------------------
  b1_summary <- do.call(rbind, lapply(results_b1, function(x) {
    ok_mse <- !is.na(x$mse) & x$mse > 0
    ok_cov <- apply(x$coverage, 1, function(r) all(!is.na(r)))
    data.frame(
      M               = x$M,
      mean_mse        = mean(x$mse[ok_mse]),
      median_mse      = median(x$mse[ok_mse]),
      rmse            = sqrt(mean(x$mse[ok_mse])),
      mean_bias       = mean(abs(x$bias[ok_mse, ])),
      coverage_95     = if (sum(ok_cov) > 0) mean(x$coverage[ok_cov, ]) else NA,
      mean_time_s     = mean(x$time[ok_mse]),
      n_ok            = sum(ok_mse)
    )
  }))

  # LaTeX table
  tex_b1 <- paste0(
    "\\begin{table}[t]\n",
    "\\centering\n",
    "\\caption{Consistency and coverage of $\\hat{\\mathbf{B}}_z$ ",
    "across corpus sizes. ",
    "Design: $K=", K_TOPICS, "$, $P=", P_COV, "$, $N=", N_VOCAB,
    "$, $\\sigma_\\varepsilon=0.3$, $", N_REP, "$ replicates. ",
    "Coverage is the empirical proportion of entries of $\\hat{\\mathbf{B}}_z$ ",
    "whose analytical 95\\% confidence interval contains the true value.}\n",
    "\\label{tab:block1}\n",
    "\\begin{tabular}{rcccccc}\n",
    "\\toprule\n",
    "$M$ & RMSE & Mean $|$bias$|$ & Coverage (95\\%) & ",
    "Time (s) & Replicates \\\\\n",
    "\\midrule\n",
    paste(apply(b1_summary, 1, function(row) {
      sprintf("%s & %.4f & %.4f & %.3f & %.2f & %d",
              format(as.integer(row["M"]), big.mark = ","),
              as.numeric(row["rmse"]),
              as.numeric(row["mean_bias"]),
              as.numeric(row["coverage_95"]),
              as.numeric(row["mean_time_s"]),
              as.integer(row["n_ok"]))
    }), collapse = " \\\\\n"),
    " \\\\\n",
    "\\bottomrule\n",
    "\\end{tabular}\n",
    "\\end{table}\n"
  )
  writeLines(tex_b1, "output/tables/block1_consistency.tex")
  cat("  Table written: output/tables/block1_consistency.tex\n")

  # --- Block 1 Figure 1: RMSE vs M (log-log) ----------------------
  b1_plot_df <- b1_summary
  b1_plot_df$M_num <- b1_plot_df$M

  p1a <- ggplot(b1_plot_df, aes(x = M_num, y = rmse)) +
    geom_point(size = 2.5, colour = "#534AB7") +
    geom_line(colour = "#534AB7", linewidth = 0.6) +
    # Reference line: M^{-1/2} rate
    geom_line(aes(y = rmse[1] * sqrt(M_num[1] / M_num)),
              linetype = "dashed", colour = "grey50", linewidth = 0.5) +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10() +
    labs(x = expression(italic(M)~"(corpus size)"),
         y = expression("RMSE of"~hat(bold(B))[z]),
         title = NULL) +
    annotate("text", x = max(b1_plot_df$M_num) * 0.7,
             y = min(b1_plot_df$rmse) * 1.5,
             label = "italic(M)^{-1/2}~rate",
             parse = TRUE,
             colour = "grey50", size = 3.5, hjust = 1)

  ggsave("output/figures/block1_rmse_vs_M.pdf", p1a,
         width = 5, height = 3.5)
  cat("  Figure written: output/figures/block1_rmse_vs_M.pdf\n")

  # --- Block 1 Figure 2: Coverage vs M -----------------------------
  p1b <- ggplot(b1_plot_df, aes(x = M_num, y = coverage_95)) +
    geom_point(size = 2.5, colour = "#1D9E75") +
    geom_line(colour = "#1D9E75", linewidth = 0.6) +
    geom_hline(yintercept = 0.95, linetype = "dashed",
               colour = "grey50", linewidth = 0.5) +
    scale_x_log10(labels = scales::comma) +
    scale_y_continuous(limits = c(0.80, 1.0),
                       breaks = seq(0.80, 1.0, 0.05)) +
    labs(x = expression(italic(M)~"(corpus size)"),
         y = "Empirical coverage (nominal 95%)",
         title = NULL) +
    annotate("text", x = min(b1_plot_df$M_num) * 1.2, y = 0.955,
             label = "Nominal 95%", colour = "grey50",
             size = 3.2, hjust = 0)

  ggsave("output/figures/block1_coverage_vs_M.pdf", p1b,
         width = 5, height = 3.5)
  cat("  Figure written: output/figures/block1_coverage_vs_M.pdf\n")

  # --- Block 1 Figure 3: Sampling distribution (QQ-plot) -----------
  # Use largest M to show normality
  M_large <- max(M_VALUES_B1)
  bias_large <- results_b1[[as.character(M_large)]]$bias
  ok_rows <- complete.cases(bias_large)
  bias_large <- bias_large[ok_rows, ]

  # Standardise each entry by its empirical SD
  bias_std <- scale(bias_large)
  # Take first column as representative
  qq_df <- data.frame(
    theoretical = qnorm(ppoints(nrow(bias_std))),
    empirical   = sort(bias_std[, 1])
  )

  p1c <- ggplot(qq_df, aes(x = theoretical, y = empirical)) +
    geom_point(size = 0.8, alpha = 0.5, colour = "#534AB7") +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "grey40") +
    coord_equal(xlim = c(-3.5, 3.5), ylim = c(-3.5, 3.5)) +
    labs(x = "Theoretical quantiles (standard normal)",
         y = "Standardised empirical quantiles",
         title = NULL,
         subtitle = bquote(italic(M) == .(format(M_large, big.mark = ","))))

  ggsave("output/figures/block1_qqplot.pdf", p1c,
         width = 4, height = 4)
  cat("  Figure written: output/figures/block1_qqplot.pdf\n")

  cat("  Block 1 complete.\n")
}


# ===================================================================
#  BLOCK 2: Linearisation Error Bound
# ===================================================================

if (2L %in% RUN_BLOCK) {
  cat("\n====== BLOCK 2: Linearisation error ======\n")

  M_FIX     <- 2000L
  K_VALUES  <- c(3L, 5L, 10L)
  SIG_VALUES <- c(0.05, 0.1, 0.2, 0.3, 0.5, 0.8)
  N_REP_B2  <- if (QUICK) 20L else 100L

  results_b2 <- list()
  row_idx <- 0L

  for (K in K_VALUES) {
    for (sig in SIG_VALUES) {
      row_idx <- row_idx + 1L
      cat(sprintf("  K=%d, sigma_eps=%.2f ...\n", K, sig))

      mse_lin <- numeric(N_REP_B2)
      bound_th <- numeric(N_REP_B2)
      max_znorm <- numeric(N_REP_B2)

      Bz0_k <- matrix(runif(P_COV * (K - 1L), -0.5, 0.5),
                       P_COV, K - 1L)

      for (r in seq_len(N_REP_B2)) {
        dat <- sim_dgp(M = M_FIX, N = N_VOCAB, K = K, P = P_COV,
                       Bz0 = Bz0_k, sigma_eps = sig,
                       alpha_beta = 0.1, doc_length = 200L,
                       seed = 20000L + row_idx * 1000L + r)

        fit <- tryCatch(
          sgscatm(dat$W, dat$C, K = K, lambda = 1, rotate = FALSE),
          error = function(e) NULL
        )
        if (is.null(fit)) next

        el <- eval_linearisation(fit, K = K)
        mse_lin[r]   <- el$mse_linearisation
        bound_th[r]  <- el$theoretical_bound
        max_znorm[r] <- el$max_z_norm
      }

      results_b2[[row_idx]] <- data.frame(
        K = K, sigma_eps = sig,
        mse_mean   = mean(mse_lin[mse_lin > 0], na.rm = TRUE),
        mse_median = median(mse_lin[mse_lin > 0], na.rm = TRUE),
        bound_mean = mean(bound_th[bound_th > 0], na.rm = TRUE),
        ratio_mean = mean(mse_lin[mse_lin > 0] / bound_th[bound_th > 0],
                          na.rm = TRUE),
        max_znorm_mean = mean(max_znorm[max_znorm > 0], na.rm = TRUE)
      )
    }
  }

  b2_df <- do.call(rbind, results_b2)
  saveRDS(b2_df, "output/data/block2_results.rds")

  # --- Block 2 LaTeX Table -----------------------------------------
  tex_b2 <- paste0(
    "\\begin{table}[t]\n",
    "\\centering\n",
    "\\caption{Linearisation error vs.\\ theoretical bound ",
    "(Proposition~\\ref{prop:linearisation_error}). ",
    "$M=", format(M_FIX, big.mark = ","), "$, $N=", N_VOCAB,
    "$, $", N_REP_B2, "$ replicates per cell. ",
    "Ratio $< 1$ confirms the bound holds; smaller ratio indicates ",
    "conservatism.}\n",
    "\\label{tab:block2}\n",
    "\\begin{tabular}{cccccc}\n",
    "\\toprule\n",
    "$K$ & $\\sigma_\\varepsilon$ & MSE (actual) & ",
    "Bound (Prop.~\\ref{prop:linearisation_error}) & ",
    "Ratio & $\\max_i\\|\\mathbf{z}_i^*\\|$ \\\\\n",
    "\\midrule\n",
    paste(apply(b2_df, 1, function(row) {
      sprintf("%d & %.2f & %.2e & %.2e & %.3f & %.3f",
              as.integer(row["K"]),
              as.numeric(row["sigma_eps"]),
              as.numeric(row["mse_mean"]),
              as.numeric(row["bound_mean"]),
              as.numeric(row["ratio_mean"]),
              as.numeric(row["max_znorm_mean"]))
    }), collapse = " \\\\\n"),
    " \\\\\n",
    "\\bottomrule\n",
    "\\end{tabular}\n",
    "\\end{table}\n"
  )
  writeLines(tex_b2, "output/tables/block2_linearisation.tex")
  cat("  Table written: output/tables/block2_linearisation.tex\n")

  # --- Block 2 Figure: MSE vs sigma_eps, faceted by K --------------
  b2_df$K_label <- paste0("K = ", b2_df$K)

  p2 <- ggplot(b2_df, aes(x = sigma_eps)) +
    geom_line(aes(y = mse_mean, colour = "Actual MSE"), linewidth = 0.7) +
    geom_point(aes(y = mse_mean, colour = "Actual MSE"), size = 2) +
    geom_line(aes(y = bound_mean, colour = "Theoretical bound"),
              linetype = "dashed", linewidth = 0.7) +
    geom_point(aes(y = bound_mean, colour = "Theoretical bound"),
               size = 2, shape = 2) +
    facet_wrap(~K_label, scales = "free_y") +
    scale_y_log10() +
    scale_colour_manual(values = c("Actual MSE" = "#534AB7",
                                   "Theoretical bound" = "#D85A30")) +
    labs(x = expression(sigma[epsilon]~"(residual noise)"),
         y = "Mean squared entry-wise error",
         colour = NULL)

  ggsave("output/figures/block2_linearisation.pdf", p2,
         width = 8, height = 3.5)
  cat("  Figure written: output/figures/block2_linearisation.pdf\n")
  cat("  Block 2 complete.\n")
}


# ===================================================================
#  BLOCK 3: Structural Comparison with STM
# ===================================================================

if (3L %in% RUN_BLOCK) {
  cat("\n====== BLOCK 3: Comparison with STM ======\n")

  # STM is run in a subprocess via the worker script block3_stm_worker.R,
  # which only loads stm (no tidyverse).
  R451 <- Sys.which("Rscript")
  if (nchar(R451) == 0L) R451 <- file.path(R.home("bin"), "Rscript")
  STM_WORKER <- "scripts/simulation/block3_stm_worker.R"
  STM_AVAILABLE <- nchar(R451) > 0L && file.exists(STM_WORKER) &&
    requireNamespace("stm", quietly = TRUE)
  if (!STM_AVAILABLE)
    cat("  NOTE: R 4.5.1 or worker script not found — STM skipped.\n")

  N_REP_B3    <- if (QUICK) 20L else 50L
  SIGNAL_LEVELS <- c(weak = 0.15, strong = 0.50)

  results_b3 <- list()
  row_idx <- 0L

  for (M in M_VALUES_B3) {
    for (sig_name in names(SIGNAL_LEVELS)) {
      b_max <- SIGNAL_LEVELS[sig_name]
      row_idx <- row_idx + 1L
      cat(sprintf("  M=%d, signal=%s (b_max=%.2f) ...\n",
                  M, sig_name, b_max))

      metrics <- data.frame(
        rep = integer(0), method = character(0),
        mse_Bz = numeric(0), time_s = numeric(0),
        stringsAsFactors = FALSE
      )

      for (r in seq_len(N_REP_B3)) {
        if (r %% 10L == 0L) cat(sprintf("    rep %d/%d\n", r, N_REP_B3))

        Bz0_r <- matrix(runif(P_COV * (K_TOPICS - 1L), -b_max, b_max),
                        P_COV, K_TOPICS - 1L)

        dat <- sim_dgp(M = M, N = N_VOCAB, K = K_TOPICS, P = P_COV,
                       Bz0 = Bz0_r, sigma_eps = 0.3,
                       alpha_beta = 0.1, doc_length = 200L,
                       seed = 30000L + row_idx * 1000L + r)

        # --- sgscatm ---
        t0 <- proc.time()
        fit_eg <- tryCatch(
          sgscatm(dat$W, dat$C, K = K_TOPICS, lambda = 1, rotate = TRUE),
          error = function(e) NULL
        )
        t_eg <- (proc.time() - t0)[3]

        if (!is.null(fit_eg)) {
          pa_eg <- procrustes_align(fit_eg$Bz, Bz0_r)
          metrics <- rbind(metrics, data.frame(
            rep = r, method = "sgscatm",
            mse_Bz = pa_eg$mse, time_s = t_eg
          ))
        }
      }

      # --- STM: run via subprocess (R 4.5.1 only) ---
      if (STM_AVAILABLE) {
        stm_rds <- tempfile(fileext = ".rds")
        cmd <- sprintf(
          '"%s" "%s" %d %s %.4f %d %d %d %d %d "%s"',
          R451, STM_WORKER,
          M, sig_name, b_max, N_REP_B3,
          K_TOPICS, P_COV, N_VOCAB, row_idx, stm_rds
        )
        ret <- system(cmd, ignore.stdout = FALSE, ignore.stderr = FALSE)
        if (ret == 0L && file.exists(stm_rds)) {
          stm_metrics <- readRDS(stm_rds)
          metrics <- rbind(metrics, stm_metrics)
          unlink(stm_rds)
        } else {
          cat(sprintf("    WARNING: STM worker failed (exit %d)\n", ret))
        }
      }

      results_b3[[row_idx]] <- list(
        M = M, signal = sig_name, b_max = b_max,
        metrics = metrics
      )
    }
  }

  saveRDS(results_b3, "output/data/block3_results.rds")

  # --- Block 3 Summary Table (LaTeX) -------------------------------
  b3_summary <- do.call(rbind, lapply(results_b3, function(x) {
    eg  <- x$metrics[x$metrics$method == "sgscatm", ]
    stm <- x$metrics[x$metrics$method == "STM", ]
    data.frame(
      M = x$M,
      signal = x$signal,
      mse_sgscatm  = if (nrow(eg) > 0) mean(eg$mse_Bz) else NA,
      mse_stm      = if (nrow(stm) > 0) mean(stm$mse_Bz) else NA,
      time_sgscatm = if (nrow(eg) > 0) mean(eg$time_s) else NA,
      time_stm     = if (nrow(stm) > 0) mean(stm$time_s) else NA,
      speedup      = if (nrow(stm) > 0 && nrow(eg) > 0)
                       mean(stm$time_s) / mean(eg$time_s) else NA,
      n_rep        = min(nrow(eg), nrow(stm))
    )
  }))

  tex_b3 <- paste0(
    "\\begin{table}[t]\n",
    "\\centering\n",
    "\\caption{Structural comparison: \\texttt{sgscatm} vs.\\ STM. ",
    "MSE of $\\hat{\\mathbf{B}}_z$ after Procrustes alignment, ",
    "and computation time. ",
    "$K=", K_TOPICS, "$, $P=", P_COV, "$, $N=", N_VOCAB, "$, ",
    "$\\sigma_\\varepsilon=0.3$.}\n",
    "\\label{tab:block3}\n",
    "\\begin{tabular}{clcccccc}\n",
    "\\toprule\n",
    "$M$ & Signal & \\multicolumn{2}{c}{MSE$(\\hat{\\mathbf{B}}_z)$} & ",
    "\\multicolumn{2}{c}{Time (s)} & Speedup & Reps \\\\\n",
    "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}\n",
    " & & sgscatm & STM & sgscatm & STM & & \\\\\n",
    "\\midrule\n",
    paste(apply(b3_summary, 1, function(row) {
      sprintf("%s & %s & %.4f & %.4f & %.1f & %.1f & %.1f$\\times$ & %d",
              format(as.integer(row["M"]), big.mark = ","),
              as.character(row["signal"]),
              as.numeric(row["mse_sgscatm"]),
              as.numeric(row["mse_stm"]),
              as.numeric(row["time_sgscatm"]),
              as.numeric(row["time_stm"]),
              as.numeric(row["speedup"]),
              as.integer(row["n_rep"]))
    }), collapse = " \\\\\n"),
    " \\\\\n",
    "\\bottomrule\n",
    "\\end{tabular}\n",
    "\\end{table}\n"
  )
  writeLines(tex_b3, "output/tables/block3_comparison.tex")
  cat("  Table written: output/tables/block3_comparison.tex\n")

  # --- Block 3 Figure: MSE comparison (boxplots) -------------------
  all_metrics <- do.call(rbind, lapply(results_b3, function(x) {
    x$metrics$M <- x$M
    x$metrics$signal <- x$signal
    x$metrics
  }))

  all_metrics$M_label <- paste0("M = ", format(all_metrics$M, big.mark = ","))
  all_metrics$M_label <- factor(all_metrics$M_label,
    levels = paste0("M = ", format(sort(unique(all_metrics$M)), big.mark = ",")))
  all_metrics$signal <- factor(all_metrics$signal,
    levels = c("weak", "strong"),
    labels = c("Signal: weak", "Signal: strong"))

  p3a <- ggplot(all_metrics, aes(x = method, y = mse_Bz, fill = method)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.8) +
    facet_grid(signal ~ M_label, scales = "free_y", switch = "y") +
    scale_fill_manual(values = c("sgscatm" = "#534AB7", "STM" = "#D85A30")) +
    scale_y_log10() +
    labs(x = NULL, y = expression("MSE of"~hat(bold(B))[z]),
         fill = NULL) +
    theme(
      legend.position = "none",
      strip.text.y.left = element_text(angle = 0, hjust = 1),
      strip.placement = "outside"
    )

  ggsave("output/figures/block3_mse_boxplot.pdf", p3a,
         width = 7, height = 5)
  cat("  Figure written: output/figures/block3_mse_boxplot.pdf\n")

  # --- Block 3 Figure: Timing comparison ---------------------------
  time_summary <- all_metrics %>%
    group_by(M, signal, method) %>%
    summarise(mean_time = mean(time_s), .groups = "drop")

  p3b <- ggplot(time_summary, aes(x = M, y = mean_time,
                                   colour = method, linetype = signal)) +
    geom_point(size = 2.5) +
    geom_line(linewidth = 0.7) +
    scale_x_log10(labels = scales::comma) +
    scale_y_log10() +
    scale_colour_manual(values = c("sgscatm" = "#534AB7",
                                    "STM" = "#D85A30")) +
    labs(x = expression(italic(M)~"(corpus size)"),
         y = "Mean computation time (s, log scale)",
         colour = "Method", linetype = "Signal")

  ggsave("output/figures/block3_timing.pdf", p3b,
         width = 6, height = 3.5)
  cat("  Figure written: output/figures/block3_timing.pdf\n")
  cat("  Block 3 complete.\n")
}


# ===================================================================
#  BLOCK 4: Four-Method Comparison — Quality and Speed
# ===================================================================
#
#  Methods compared:
#    sgscatm      – ILR-spectral baseline (closed-form, no EM)
#    sgscatm_ref  – spectral + one refine_phi step (k-means M-step)
#    stm          – variational EM with default Spectral init
#    stm_warm     – variational EM warm-started from sgscatm solution
#
#  Design axes:
#    M        – corpus size  (M_VALUES_B4)
#    scenario – "high_sep", "low_sep", "weak_sig"  (SCENARIOS_B4)
#
#  Metrics recorded per replicate:
#    mse_Bz  – MSE of path coefficients after Procrustes alignment
#    mse_phi – MSE of topic-word distributions under optimal row perm
#    time_s  – total wall time (stm_warm includes sgscatm time)
#    n_iter  – EM iterations to convergence (STM methods only)
# ===================================================================

if (4L %in% RUN_BLOCK) {
  cat("\n====== BLOCK 4: Four-Method Comparison ======\n")

  # Locate Rscript executable
  R_exe_b4 <- Sys.which("Rscript")
  if (nchar(R_exe_b4) == 0L) R_exe_b4 <- file.path(R.home("bin"), "Rscript")

  WORKER_B4 <- "scripts/simulation/block4_comparison_worker.R"
  STM_B4_OK <- nchar(R_exe_b4) > 0L &&
    file.exists(WORKER_B4) &&
    requireNamespace("stm", quietly = TRUE)

  if (!STM_B4_OK)
    cat("  WARNING: stm package or worker not found — Block 4 skipped.\n")

  results_b4 <- list()
  row_idx    <- 0L

  if (STM_B4_OK) {
    for (M in M_VALUES_B4) {
      for (scenario in SCENARIOS_B4) {
        row_idx <- row_idx + 1L
        cat(sprintf("  M=%d, scenario=%s ...\n", M, scenario))

        out_rds <- tempfile(fileext = ".rds")
        cmd <- sprintf(
          '"%s" "%s" %d %s %d %d %d %d %d "%s"',
          R_exe_b4, WORKER_B4,
          M, scenario, N_REP_B4,
          K_TOPICS, P_COV, N_VOCAB, row_idx, out_rds
        )
        ret <- system(cmd, ignore.stdout = FALSE, ignore.stderr = FALSE)

        cell_metrics <- NULL
        if (ret == 0L && file.exists(out_rds)) {
          cell_metrics <- readRDS(out_rds)
          unlink(out_rds)
        } else {
          cat(sprintf("    WARNING: worker failed (exit %d)\n", ret))
        }

        results_b4[[row_idx]] <- list(
          M        = M,
          scenario = scenario,
          metrics  = cell_metrics
        )
      }
    }
  }

  saveRDS(results_b4, "output/data/block4_results.rds")

  # --- Block 4 Summary Table (LaTeX) --------------------------------
  # One row per (M, scenario, method); grouped by scenario in the table.
  METHOD_ORDER <- c("sgscatm", "sgscatm_ref", "stm", "stm_warm")
  METHOD_LABEL <- c(
    sgscatm     = "\\texttt{sgscatm}",
    sgscatm_ref = "\\texttt{sgscatm}+refine",
    stm         = "\\texttt{stm}",
    stm_warm    = "\\texttt{stm}+warm"
  )
  SCENARIO_LABEL <- c(
    high_sep = "High sep.",
    low_sep  = "Low sep.",
    weak_sig = "Weak signal"
  )

  b4_rows <- do.call(rbind, lapply(results_b4, function(x) {
    if (is.null(x$metrics) || nrow(x$metrics) == 0L) return(NULL)
    do.call(rbind, lapply(METHOD_ORDER, function(mth) {
      m <- x$metrics[x$metrics$method == mth, ]
      if (nrow(m) == 0L) return(NULL)
      data.frame(
        M            = x$M,
        scenario     = x$scenario,
        method       = mth,
        mean_mse_Bz  = mean(m$mse_Bz,  na.rm = TRUE),
        sd_mse_Bz    = sd(m$mse_Bz,    na.rm = TRUE),
        mean_mse_phi = mean(m$mse_phi,  na.rm = TRUE),
        sd_mse_phi   = sd(m$mse_phi,    na.rm = TRUE),
        mean_time    = mean(m$time_s,   na.rm = TRUE),
        mean_n_iter  = mean(m$n_iter,   na.rm = TRUE),
        n_rep        = sum(!is.na(m$mse_Bz)),
        stringsAsFactors = FALSE
      )
    }))
  }))

  if (!is.null(b4_rows) && nrow(b4_rows) > 0L) {
    # Render one sub-block per scenario
    body_lines <- character(0)
    for (sc_name in SCENARIOS_B4) {
      sc_rows <- b4_rows[b4_rows$scenario == sc_name, ]
      if (nrow(sc_rows) == 0L) next

      body_lines <- c(body_lines,
        sprintf("\\multicolumn{7}{l}{\\textit{Scenario: %s}} \\\\",
                SCENARIO_LABEL[sc_name]),
        "\\midrule"
      )

      for (M_val in M_VALUES_B4) {
        m_rows <- sc_rows[sc_rows$M == M_val, ]
        first_m <- TRUE
        for (mth in METHOD_ORDER) {
          r_row <- m_rows[m_rows$method == mth, ]
          if (nrow(r_row) == 0L) next

          n_iter_str <- if (!is.na(r_row$mean_n_iter) &&
                            mth %in% c("stm", "stm_warm"))
            sprintf("%.0f", r_row$mean_n_iter) else "---"

          body_lines <- c(body_lines,
            sprintf(
              "%s & %s & %.4f$\\pm$%.4f & %.2e$\\pm$%.2e & %.1f & %s & %d \\\\",
              if (first_m) format(M_val, big.mark = ",") else "",
              METHOD_LABEL[mth],
              r_row$mean_mse_Bz, r_row$sd_mse_Bz,
              r_row$mean_mse_phi, r_row$sd_mse_phi,
              r_row$mean_time,
              n_iter_str,
              r_row$n_rep
            )
          )
          first_m <- FALSE
        }
        body_lines <- c(body_lines, "\\addlinespace")
      }
      body_lines <- c(body_lines, "\\midrule")
    }

    tex_b4 <- paste0(
      "\\begin{table}[t]\n",
      "\\centering\n",
      "\\small\n",
      "\\caption{Four-method comparison across corpus sizes and scenarios. ",
      "MSE$(\\hat{\\mathbf{B}}_z)$: mean squared error of path coefficients ",
      "after Procrustes alignment (mean$\\pm$sd). ",
      "MSE$(\\hat{\\boldsymbol{\\Phi}})$: topic-word MSE under optimal row ",
      "permutation. Time: total wall time in seconds. Iter: mean EM iterations ",
      "to convergence (STM methods only). ",
      "$K=", K_TOPICS, "$, $P=", P_COV, "$, $N=", N_VOCAB, "$, ",
      N_REP_B4, " replicates per cell.}\n",
      "\\label{tab:block4}\n",
      "\\begin{tabular}{rlccccc}\n",
      "\\toprule\n",
      "$M$ & Method & MSE$(\\hat{\\mathbf{B}}_z)$ & ",
      "MSE$(\\hat{\\boldsymbol{\\Phi}})$ & ",
      "Time (s) & Iter & Reps \\\\\n",
      "\\midrule\n",
      paste(body_lines, collapse = "\n"),
      "\n",
      "\\bottomrule\n",
      "\\end{tabular}\n",
      "\\end{table}\n"
    )
    writeLines(tex_b4, "output/tables/block4_comparison.tex")
    cat("  Table written: output/tables/block4_comparison.tex\n")
  }

  # --- Block 4 figures ----------------------------------------------
  if (!is.null(b4_rows) && nrow(b4_rows) > 0L) {

    all_b4_metrics <- do.call(rbind, lapply(results_b4, function(x) {
      if (is.null(x$metrics) || nrow(x$metrics) == 0L) return(NULL)
      x$metrics$M        <- x$M
      x$metrics$scenario <- x$scenario
      x$metrics
    }))

    # Factor ordering for display
    all_b4_metrics$method <- factor(all_b4_metrics$method,
      levels = METHOD_ORDER,
      labels = c("sgscatm", "sgscatm+refine", "stm", "stm+warm"))
    all_b4_metrics$M_label <- factor(
      paste0("M=", format(all_b4_metrics$M, big.mark = ",")),
      levels = paste0("M=", format(sort(unique(all_b4_metrics$M)),
                                   big.mark = ",")))
    all_b4_metrics$scenario <- factor(all_b4_metrics$scenario,
      levels = SCENARIOS_B4,
      labels = c("High sep.", "Low sep.", "Weak signal"))

    METHOD_COLOURS <- c(
      "sgscatm"       = "#534AB7",
      "sgscatm+refine"= "#1D9E75",
      "stm"           = "#D85A30",
      "stm+warm"      = "#C49A00"
    )

    # Figure 4a: MSE(Bz) boxplot — faceted by scenario x M
    p4a <- ggplot(all_b4_metrics,
                  aes(x = method, y = mse_Bz, fill = method)) +
      geom_boxplot(alpha = 0.7, outlier.size = 0.6) +
      facet_grid(scenario ~ M_label, scales = "free_y", switch = "y") +
      scale_fill_manual(values = METHOD_COLOURS) +
      scale_y_log10() +
      labs(x = NULL,
           y = expression("MSE of"~hat(bold(B))[z]~"(log scale)"),
           fill = NULL) +
      theme(
        axis.text.x  = element_text(angle = 35, hjust = 1, size = 8),
        legend.position = "none",
        strip.text.y.left = element_text(angle = 0, hjust = 1),
        strip.placement = "outside"
      )

    ggsave("output/figures/block4_mse_Bz.pdf", p4a,
           width = 8, height = 6)
    cat("  Figure written: output/figures/block4_mse_Bz.pdf\n")

    # Figure 4b: MSE(Phi) boxplot — same layout
    p4b <- ggplot(all_b4_metrics,
                  aes(x = method, y = mse_phi, fill = method)) +
      geom_boxplot(alpha = 0.7, outlier.size = 0.6) +
      facet_grid(scenario ~ M_label, scales = "free_y", switch = "y") +
      scale_fill_manual(values = METHOD_COLOURS) +
      scale_y_log10() +
      labs(x = NULL,
           y = expression("MSE of"~hat(bold(Phi))~"(log scale)"),
           fill = NULL) +
      theme(
        axis.text.x  = element_text(angle = 35, hjust = 1, size = 8),
        legend.position = "none",
        strip.text.y.left = element_text(angle = 0, hjust = 1),
        strip.placement = "outside"
      )

    ggsave("output/figures/block4_mse_phi.pdf", p4b,
           width = 8, height = 6)
    cat("  Figure written: output/figures/block4_mse_phi.pdf\n")

    # Figure 4c: Mean computation time vs M (log-log), by method and scenario
    time_b4 <- all_b4_metrics %>%
      group_by(M, scenario, method) %>%
      summarise(mean_time = mean(time_s, na.rm = TRUE), .groups = "drop")

    p4c <- ggplot(time_b4,
                  aes(x = M, y = mean_time, colour = method, shape = method)) +
      geom_point(size = 2.2) +
      geom_line(linewidth = 0.65) +
      facet_wrap(~ scenario, ncol = 3L) +
      scale_x_log10(labels = scales::comma) +
      scale_y_log10() +
      scale_colour_manual(values = METHOD_COLOURS) +
      scale_shape_manual(values = c(16L, 17L, 15L, 18L)) +
      labs(x = expression(italic(M)~"(corpus size)"),
           y = "Mean computation time (s, log scale)",
           colour = NULL, shape = NULL)

    ggsave("output/figures/block4_timing.pdf", p4c,
           width = 9, height = 3.5)
    cat("  Figure written: output/figures/block4_timing.pdf\n")

    # Figure 4d: EM iterations for stm vs stm+warm (convergence speedup)
    iter_b4 <- all_b4_metrics %>%
      filter(method %in% c("stm", "stm+warm"), !is.na(n_iter)) %>%
      group_by(M, scenario, method) %>%
      summarise(mean_iter = mean(n_iter, na.rm = TRUE),
                sd_iter   = sd(n_iter,   na.rm = TRUE),
                .groups = "drop")

    if (nrow(iter_b4) > 0L) {
      p4d <- ggplot(iter_b4,
                    aes(x = M, y = mean_iter,
                        colour = method, linetype = method)) +
        geom_ribbon(aes(ymin = pmax(mean_iter - sd_iter, 0),
                        ymax = mean_iter + sd_iter,
                        fill = method),
                    alpha = 0.15, colour = NA) +
        geom_line(linewidth = 0.7) +
        geom_point(size = 2.2) +
        facet_wrap(~ scenario, ncol = 3L) +
        scale_x_log10(labels = scales::comma) +
        scale_colour_manual(values = METHOD_COLOURS[c("stm", "stm+warm")]) +
        scale_fill_manual(  values = METHOD_COLOURS[c("stm", "stm+warm")]) +
        labs(x = expression(italic(M)~"(corpus size)"),
             y = "Mean EM iterations to convergence",
             colour = NULL, linetype = NULL, fill = NULL)

      ggsave("output/figures/block4_em_iters.pdf", p4d,
             width = 9, height = 3.5)
      cat("  Figure written: output/figures/block4_em_iters.pdf\n")
    }
  }

  cat("  Block 4 complete.\n")
}


# ===================================================================
#  Final summary
# ===================================================================

cat("\n====== Simulation complete ======\n")
cat("Output files:\n")
for (d in c("output/tables", "output/figures", "output/data")) {
  files <- list.files(d, full.names = TRUE)
  if (length(files) > 0) cat(paste("  ", files, collapse = "\n"), "\n")
}
cat("\nTo include tables in LaTeX:\n")
cat("  \\input{output/tables/block1_consistency}\n")
cat("  \\input{output/tables/block2_linearisation}\n")
cat("  \\input{output/tables/block3_comparison}\n")
cat("  \\input{output/tables/block4_comparison}\n")
