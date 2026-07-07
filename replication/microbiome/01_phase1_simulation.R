#!/usr/bin/env Rscript
# ===================================================================
# Phase 1 — Microbiome-calibrated simulation sweep
#
# Validates the corrected analytical SE (sgscatm_vcov) in the CoDa-favorable
# regime and traces the delocalization boundary (Proposition 16 crossover).
#
# Estimand handling: B_z is identified only up to rotation/sign and the
# estimator standardizes scores to unit variance, so B_z is NOT on the raw
# generative Bz0 scale. All SE/coverage metrics are therefore scale-free:
#   * calibration  = mean analytical SE vs empirical sampling SD across reps
#   * coverage     = of the generalized-Procrustes across-replicate mean
# Recovery of the true Bz0 is reported separately via (i) canonical
# correlations between estimated and true ILR scores and (ii) RMSE after an
# optimal linear (scale+rotation) score map onto the truth.
#
# Outputs: output/phase1_cells.rds, output/tables/phase1_*.csv
# ===================================================================
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
source("replication/simulation/sim_dgp.R")
source("replication/simulation/sim_utils.R")
dir.create("output/tables", recursive = TRUE, showWarnings = FALSE)

set.seed(2026)
lam <- 1

# NegBin library-size function (microbiome-like: overdispersed, mean ~2e4)
libsize_fun <- function(m) pmax(rnbinom(m, size = 3, mu = 2e4), 500L)

align_to <- function(B, Ref) { sv <- svd(crossprod(B, Ref)); B %*% (sv$u %*% t(sv$v)) }

# ---- evaluate one cell: R replicates, return per-cell summary + raw ----
run_cell <- function(M, N, K, P, b_max, R, cell_id,
                     do_boot = FALSE, B_boot = 80L, R_boot = 15L) {
  Km1 <- K - 1L
  # fixed true coefficients (scaled to b_max) so the estimand is stable
  set.seed(1000L + cell_id)
  Bz0 <- matrix(runif(P * Km1, -b_max, b_max), P, Km1)

  Bstd <- vector("list", R); fitobj <- vector("list", R)
  rvec <- numeric(R); cc_min <- numeric(R); rmse_rec <- numeric(R)
  for (r in seq_len(R)) {
    dat <- sim_dgp(M = M, N = N, K = K, P = P, Bz0 = Bz0,
                   sigma_eps = 0.3, alpha_beta = 0.05,
                   doc_length = libsize_fun, seed = cell_id * 10000L + r)
    fit <- sgscatm(dat$W, dat$C, K = K, lambda = lam, rotate = FALSE)
    fitobj[[r]] <- fit
    Bstd[[r]]   <- sgscatm_vcov(fit, identified = FALSE)$B
    rvec[r]     <- M * max(rowSums(fit$Z^2)) / Km1
    # score recovery vs truth
    cc <- tryCatch(cancor(fit$Z, dat$Z_true)$cor, error = function(e) NA_real_)
    cc_min[r] <- if (all(is.finite(cc))) min(cc) else NA_real_
    Amap <- tryCatch(solve(crossprod(fit$Z), crossprod(fit$Z, dat$Z_true)),
                     error = function(e) NULL)
    if (!is.null(Amap)) {
      Bz_rec <- solve(crossprod(dat$C), crossprod(dat$C, fit$Z %*% Amap))
      rmse_rec[r] <- sqrt(mean((Bz_rec - Bz0)^2))
    } else rmse_rec[r] <- NA_real_
  }

  # generalized Procrustes to common reference (estimand)
  Ref <- Bstd[[1]]
  for (pass in 1:3) {
    Bal <- lapply(Bstd, align_to, Ref = Ref); Ref <- Reduce(`+`, Bal) / R
  }
  Bal <- lapply(Bstd, align_to, Ref = Ref)
  arr    <- array(unlist(Bal), dim = c(P, Km1, R))
  emp_sd <- apply(arr, c(1, 2), sd)

  se_list <- vector("list", R); cov_arr <- array(NA, dim = c(P, Km1, R))
  for (r in seq_len(R)) {
    Rr <- { sv <- svd(crossprod(Bstd[[r]], Ref)); sv$u %*% t(sv$v) }
    vc <- sgscatm_vcov(fitobj[[r]], rotation = Rr, identified = TRUE)
    se_list[[r]]   <- vc$se
    cov_arr[, , r] <- (Ref >= Bal[[r]] - 1.96 * vc$se) &
                      (Ref <= Bal[[r]] + 1.96 * vc$se)
  }
  se_mean <- Reduce(`+`, se_list) / R
  ratio   <- se_mean / emp_sd

  # optional bootstrap-vs-analytical SE at this cell (reduced reps)
  boot_med_ratio <- NA_real_
  if (do_boot) {
    br <- numeric(R_boot)
    for (r in seq_len(R_boot)) {
      dat <- sim_dgp(M = M, N = N, K = K, P = P, Bz0 = Bz0,
                     sigma_eps = 0.3, alpha_beta = 0.05,
                     doc_length = libsize_fun, seed = cell_id * 77000L + r)
      fit <- sgscatm(dat$W, dat$C, K = K, lambda = lam, rotate = FALSE)
      # bootstrap: refit on resampled docs, sign/rotation-align each to fit
      Bref <- sgscatm_vcov(fit, identified = FALSE)$B
      Bb <- matrix(NA_real_, B_boot, P * Km1)
      for (b in seq_len(B_boot)) {
        idx <- sample.int(M, M, replace = TRUE)
        fb <- tryCatch(sgscatm(dat$W[idx, ], dat$C[idx, ], K = K,
                               lambda = lam, rotate = FALSE),
                       error = function(e) NULL)
        if (!is.null(fb)) {
          Bbs <- sgscatm_vcov(fb, identified = FALSE)$B
          Bb[b, ] <- as.vector(align_to(Bbs, Bref))
        }
      }
      boot_sd <- matrix(apply(Bb, 2, sd, na.rm = TRUE), P, Km1)
      vc <- sgscatm_vcov(fit, identified = TRUE)$se
      br[r] <- median(vc / boot_sd, na.rm = TRUE)
    }
    boot_med_ratio <- mean(br, na.rm = TRUE)
  }

  data.frame(
    cell_id = cell_id, M = M, N = N, K = K, P = P, b_max = b_max, R = R,
    r_deloc      = mean(rvec),
    r_deloc_sd   = sd(rvec),
    coverage     = mean(cov_arr),
    se_sd_ratio  = median(ratio),
    emp_sd_mean  = mean(emp_sd),
    ana_se_mean  = mean(se_mean),
    cc_min_mean  = mean(cc_min, na.rm = TRUE),
    rmse_recover = mean(rmse_rec, na.rm = TRUE),
    boot_se_ratio = boot_med_ratio
  )
}

# ================= design =================
P <- 4L
cells <- list()
cid <- 0L

# (A) primary b_max sweep at favorable dims -> G3 crossover figure
bmax_grid <- c(0.10, 0.25, 0.50, 0.75, 1.00, 1.50, 2.00)
for (bm in bmax_grid) {
  cid <- cid + 1L
  cells[[length(cells) + 1L]] <- list(M = 2000L, N = 200L, K = 5L, P = P,
                                       b_max = bm, R = 60L, cell_id = cid,
                                       do_boot = bm %in% c(0.25, 1.50))
}

# (B) dimension-robustness grid at favorable & stressed b_max
dim_grid <- expand.grid(N = c(200L, 500L), M = c(500L, 2000L), K = c(5L, 8L))
for (i in seq_len(nrow(dim_grid))) for (bm in c(0.25, 1.00)) {
  cid <- cid + 1L
  cells[[length(cells) + 1L]] <- list(M = dim_grid$M[i], N = dim_grid$N[i],
                                       K = dim_grid$K[i], P = P, b_max = bm,
                                       R = 50L, cell_id = cid, do_boot = FALSE)
}

cat(sprintf("Phase 1: %d cells\n", length(cells)))
res <- vector("list", length(cells))
for (i in seq_along(cells)) {
  cc <- cells[[i]]
  t0 <- proc.time()[3]
  res[[i]] <- do.call(run_cell, cc)
  cat(sprintf("[%2d/%2d] M=%d N=%d K=%d bmax=%.2f  r=%.2f cov=%.3f SE/SD=%.2f ccmin=%.2f  [%.1fs]\n",
              i, length(cells), cc$M, cc$N, cc$K, cc$b_max,
              res[[i]]$r_deloc, res[[i]]$coverage, res[[i]]$se_sd_ratio,
              res[[i]]$cc_min_mean, proc.time()[3] - t0))
  saveRDS(do.call(rbind, res[seq_len(i)]), "output/phase1_cells.rds")
}

df <- do.call(rbind, res)
write.csv(df, "output/tables/phase1_all_cells.csv", row.names = FALSE)
write.csv(df[df$M == 2000 & df$N == 200 & df$K == 5, ],
          "output/tables/phase1_bmax_sweep.csv", row.names = FALSE)
cat("Phase 1 DONE. Wrote output/phase1_cells.rds and tables.\n")
