#!/usr/bin/env Rscript
# ===================================================================
#  Block 3 STM Worker — runs on R 4.5.1 (only stm, no tidyverse)
#
#  Called by run_simulation.R via system():
#    Rscript block3_stm_worker.R <M> <sig_name> <b_max> <N_REP>
#           <K> <P> <N_VOCAB> <row_idx> <out_rds>
#
#  Writes a data.frame with columns:
#    rep, method="STM", mse_Bz, time_s
#  to <out_rds>.
# ===================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 9L) stop("Usage: block3_stm_worker.R M sig_name b_max N_REP K P N_VOCAB row_idx out_rds")

M_val    <- as.integer(args[1L])
sig_name <- args[2L]
b_max    <- as.numeric(args[3L])
N_REP    <- as.integer(args[4L])
K_TOPICS <- as.integer(args[5L])
P_COV    <- as.integer(args[6L])
N_VOCAB  <- as.integer(args[7L])
row_idx  <- as.integer(args[8L])
out_rds  <- args[9L]

# Source egscatm infrastructure (no ggplot2/dplyr dependencies)
for (.f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(.f)
source("scripts/simulation/sim_dgp.R")
source("scripts/simulation/sim_utils.R")

suppressPackageStartupMessages(library(stm))

metrics <- data.frame(
  rep = integer(0), method = character(0),
  mse_Bz = numeric(0), time_s = numeric(0),
  stringsAsFactors = FALSE
)

for (r in seq_len(N_REP)) {
  Bz0_r <- matrix(runif(P_COV * (K_TOPICS - 1L), -b_max, b_max),
                  P_COV, K_TOPICS - 1L)

  dat <- sim_dgp(M = M_val, N = N_VOCAB, K = K_TOPICS, P = P_COV,
                 Bz0 = Bz0_r, sigma_eps = 0.3,
                 alpha_beta = 0.1, doc_length = 200L,
                 seed = 30000L + row_idx * 1000L + r)

  # --- Convert DTM to STM list format (filter zero-count terms) ---
  active_cols <- which(colSums(dat$W) > 0L)
  W_stm <- dat$W[, active_cols, drop = FALSE]
  stm_docs <- lapply(seq_len(nrow(W_stm)), function(i) {
    idx <- which(W_stm[i, ] > 0L)
    rbind(as.integer(idx), as.integer(W_stm[i, idx]))
  })
  stm_vocab <- as.character(active_cols)

  prevalence_formula <- as.formula(
    paste0("~ ", paste0("V", seq_len(P_COV), collapse = " + "))
  )
  cov_df <- as.data.frame(dat$C)
  colnames(cov_df) <- paste0("V", seq_len(P_COV))

  t0 <- proc.time()
  fit_stm <- tryCatch(
    stm(documents = stm_docs, vocab = stm_vocab, K = K_TOPICS,
        prevalence = prevalence_formula, data = cov_df,
        init.type = "Spectral", verbose = FALSE, max.em.its = 75L),
    error = function(e) NULL
  )
  t_stm <- (proc.time() - t0)[3L]

  if (!is.null(fit_stm)) {
    stm_theta <- fit_stm$theta                           # M x K
    V_ilr     <- ilr_contrast(K_TOPICS)                  # K x (K-1)
    stm_Z     <- log(pmax(stm_theta, 1e-10)) %*% V_ilr  # M x (K-1)

    stm_Bz <- tryCatch(
      solve(crossprod(dat$C), crossprod(dat$C, stm_Z)),
      error = function(e) NULL
    )

    if (!is.null(stm_Bz)) {
      pa_stm <- procrustes_align(stm_Bz, Bz0_r)
      metrics <- rbind(metrics, data.frame(
        rep = r, method = "STM",
        mse_Bz = pa_stm$mse, time_s = t_stm
      ))
    }
  }
}

saveRDS(metrics, out_rds)
cat(sprintf("Worker done: M=%d signal=%s reps_ok=%d\n", M_val, sig_name, nrow(metrics)))
