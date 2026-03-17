#!/usr/bin/env Rscript
# ===================================================================
#  Block 4 Worker — Four-Method Comparison
#
#  Runs all four methods on the same simulated data for one
#  (M, scenario) cell, ensuring a perfectly paired comparison:
#
#    1. sgscatm       – ILR-spectral baseline
#    2. sgscatm_ref   – spectral + one-step refine_phi (k-means)
#    3. stm           – variational EM with Spectral init
#    4. stm_warm      – variational EM, warm-started from sgscatm
#
#  Arguments (positional):
#    M         – corpus size (integer)
#    scenario  – "high_sep" | "low_sep" | "weak_sig"
#    N_REP     – number of replicates
#    K         – number of topics
#    P         – number of covariates
#    N_VOCAB   – vocabulary size
#    row_idx   – grid-row index (used for reproducible seeding)
#    out_rds   – output path for the RDS result
#
#  Output: data.frame saved to <out_rds> with columns:
#    rep, method, mse_Bz, mse_phi, time_s, n_iter
#
#  Notes:
#  - mse_Bz  : MSE of Bz after Procrustes alignment to ground truth.
#  - mse_phi : MSE of Phi (topic-word) under optimal row permutation.
#  - time_s  : wall time (for stm_warm this includes sgscatm time).
#  - n_iter  : EM iterations to convergence (STM methods only).
# ===================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 8L)
  stop(paste(
    "Usage: block4_comparison_worker.R",
    "M scenario N_REP K P N_VOCAB row_idx out_rds"
  ))

M_val    <- as.integer(args[1L])
scenario <- args[2L]
N_REP    <- as.integer(args[3L])
K_TOPICS <- as.integer(args[4L])
P_COV    <- as.integer(args[5L])
N_VOCAB  <- as.integer(args[6L])
row_idx  <- as.integer(args[7L])
out_rds  <- args[8L]

# Source all sgscatm package files (no ggplot2 / dplyr dependency)
for (.f in list.files("R", pattern = "[.]R$", full.names = TRUE)) source(.f)
source("scripts/simulation/sim_dgp.R")
source("scripts/simulation/sim_utils.R")

suppressPackageStartupMessages(library(stm))

# ===================================================================
# Scenario definitions
# ===================================================================
# Three scenarios that stress-test different aspects of each method:
#   high_sep  – very distinct topics (sparse Beta), moderate signal
#   low_sep   – overlapping topics  (diffuse Beta), moderate signal
#   weak_sig  – medium topics, weak covariate signal (small Bz)
SCENARIOS <- list(
  high_sep = list(alpha_beta = 0.01, b_max = 0.30, sigma_eps = 0.3),
  low_sep  = list(alpha_beta = 1.00, b_max = 0.30, sigma_eps = 0.3),
  weak_sig = list(alpha_beta = 0.10, b_max = 0.10, sigma_eps = 0.3)
)

if (!scenario %in% names(SCENARIOS))
  stop(sprintf("Unknown scenario '%s'. Choose from: %s",
               scenario, paste(names(SCENARIOS), collapse = ", ")))
sc <- SCENARIOS[[scenario]]

# ===================================================================
# Internal helpers
# ===================================================================

# Recursively generate all permutations of integer vector x.
# For K = 5 this yields 120 permutations.
.gen_perms <- function(x) {
  if (length(x) <= 1L) return(list(x))
  perms <- list()
  for (i in seq_along(x)) {
    sub <- .gen_perms(x[-i])
    for (s in sub) perms <- c(perms, list(c(x[i], s)))
  }
  perms
}

# MSE between Phi_hat (K x N) and Beta_true (K x N) under the
# optimal row permutation.  For K <= 8 exact enumeration is fast.
best_perm_mse_phi <- function(Phi_hat, Beta_true, K) {
  all_p   <- .gen_perms(seq_len(K))
  best_mse <- Inf
  for (p in all_p) {
    m <- mean((Phi_hat[p, ] - Beta_true)^2)
    if (m < best_mse) best_mse <- m
  }
  best_mse
}

# Convert a DTM to STM's list-of-matrices document format.
# Returns a list: $docs (list of 2 x nnz matrices), $vocab, $active_cols.
to_stm_fmt <- function(W) {
  active_cols <- which(colSums(W) > 0L)
  W_act <- W[, active_cols, drop = FALSE]
  docs  <- lapply(seq_len(nrow(W_act)), function(i) {
    idx <- which(W_act[i, ] > 0L)
    rbind(as.integer(idx), as.integer(W_act[i, idx]))
  })
  list(docs = docs, vocab = as.character(active_cols),
       active_cols = active_cols)
}

# Estimate Bz from an STM fit by OLS of ILR-transformed theta on C.
# Uses the same ILR basis as sgscatm for a fair comparison.
.stm_extract_Bz <- function(fit_stm, C, K) {
  V_ilr <- ilr_contrast(K)
  stm_Z <- log(pmax(fit_stm$theta, 1e-10)) %*% V_ilr   # M x (K-1)
  tryCatch(
    solve(crossprod(C), crossprod(C, stm_Z)),
    error = function(e) NULL
  )
}

# Reconstruct a full-vocabulary K x N_vocab Phi matrix from the
# STM beta field (which covers only active vocabulary columns).
.stm_extract_Phi <- function(fit_stm, stm_fmt, N_vocab) {
  tryCatch({
    log_b <- fit_stm$beta$logbeta[[1L]]             # K x |active|
    # Row-normalise in log space for numerical stability
    phi   <- exp(log_b - apply(log_b, 1L, function(x) {
      lse <- max(x); lse + log(sum(exp(x - lse)))
    }))
    Phi_full <- matrix(0, nrow(phi), N_vocab)
    Phi_full[, stm_fmt$active_cols] <- phi
    Phi_full
  }, error = function(e) NULL)
}

# Build the restart object for STM's model= argument.
#
# STM's model= restart path (stm.control else-branch) expects:
#   model$mu         – list(mu = matrix (K-1) x M or (K-1) x 1)
#   model$sigma      – (K-1) x (K-1) covariance matrix
#   model$beta       – list(logbeta = list(K x V_active log-prob matrix))
#   model$eta        – M x (K-1) variational parameters (lambda in STM)
#   model$convergence – list with at least $its (integer)
#
# We convert sgscatm's Pi to STM's logit-normal parameterisation:
#   eta_{ik} = log(Pi_{ik} / Pi_{iK}),  k = 1, ..., K-1
# (topic K as reference category, matching STM's internal convention).
.build_stm_init <- function(fit_sg, stm_fmt, K) {
  Pi  <- pmax(fit_sg$Pi, 1e-10)
  Pi  <- Pi / rowSums(Pi)

  # Per-document logit-normal scores (M x K-1):
  #   eta_{ik} = log(Pi_{ik}/Pi_{iK}), k = 1,...,K-1
  # These become the initial variational parameters (model$eta / lambda).
  eta <- log(Pi[, -K, drop = FALSE]) - log(Pi[, K])   # M x (K-1)

  # Global mu: (K-1) x 1 — STM expects a column mean, not per-document.
  # Use the centroid of sgscatm's eta as starting global mean.
  mu_global <- matrix(colMeans(eta), ncol = 1L)       # (K-1) x 1

  # Log topic-word distributions restricted to active vocabulary
  phi_act <- fit_sg$Phi[, stm_fmt$active_cols, drop = FALSE]
  phi_act <- pmax(phi_act, 1e-10)
  phi_act <- phi_act / rowSums(phi_act)

  list(
    mu    = list(mu = mu_global),                # (K-1) x 1 global mean
    sigma = diag(K - 1L),
    beta  = list(logbeta = list(log(phi_act))),
    eta   = eta,                                 # M x (K-1) per-doc init
    convergence = list(
      bound          = numeric(0),
      its            = 0L,
      stopits        = FALSE,
      converged      = FALSE,
      allow.neg.change = TRUE
    )
  )
}

# ===================================================================
# Main replication loop
# ===================================================================
prevalence_formula <- as.formula(
  paste0("~ ", paste0("V", seq_len(P_COV), collapse = " + "))
)

metrics <- data.frame(
  rep    = integer(0), method  = character(0),
  mse_Bz = numeric(0), mse_phi = numeric(0),
  time_s = numeric(0), n_iter  = integer(0),
  stringsAsFactors = FALSE
)

for (r in seq_len(N_REP)) {

  seed_r <- 40000L + row_idx * 1000L + r
  set.seed(seed_r)

  Bz0_r <- matrix(runif(P_COV * (K_TOPICS - 1L), -sc$b_max, sc$b_max),
                  P_COV, K_TOPICS - 1L)

  dat <- sim_dgp(M = M_val, N = N_VOCAB, K = K_TOPICS, P = P_COV,
                 Bz0        = Bz0_r,
                 sigma_eps  = sc$sigma_eps,
                 alpha_beta = sc$alpha_beta,
                 doc_length = 200L,
                 seed       = seed_r)

  cov_df  <- as.data.frame(dat$C)
  colnames(cov_df) <- paste0("V", seq_len(P_COV))
  stm_fmt <- to_stm_fmt(dat$W)

  # ------------------------------------------------------------------
  # Method 1: sgscatm baseline
  # ------------------------------------------------------------------
  t0     <- proc.time()
  fit_sg <- tryCatch(
    sgscatm(dat$W, dat$C, K = K_TOPICS, lambda = 1, rotate = TRUE),
    error = function(e) NULL
  )
  t_sg <- (proc.time() - t0)[3L]

  if (!is.null(fit_sg)) {
    pa_sg   <- procrustes_align(fit_sg$Bz, Bz0_r)
    phi_sg  <- best_perm_mse_phi(fit_sg$Phi, dat$Beta, K_TOPICS)

    metrics <- rbind(metrics, data.frame(
      rep = r, method = "sgscatm",
      mse_Bz = pa_sg$mse, mse_phi = phi_sg,
      time_s = t_sg, n_iter = NA_integer_
    ))

    # ----------------------------------------------------------------
    # Method 2: sgscatm + refine_phi (k-means one-step)
    # ----------------------------------------------------------------
    t0r     <- proc.time()
    fit_ref <- tryCatch(
      refine_phi(fit_sg, dat$W, method = "kmeans", seed = seed_r),
      error = function(e) NULL
    )
    t_ref <- t_sg + (proc.time() - t0r)[3L]   # cumulative incl. base

    if (!is.null(fit_ref)) {
      phi_ref <- best_perm_mse_phi(fit_ref$Phi, dat$Beta, K_TOPICS)
      metrics <- rbind(metrics, data.frame(
        rep = r, method = "sgscatm_ref",
        mse_Bz = pa_sg$mse,   # Bz is unchanged by refine_phi
        mse_phi = phi_ref,
        time_s = t_ref, n_iter = NA_integer_
      ))
    }
  }

  # ------------------------------------------------------------------
  # Method 3: STM with Spectral initialisation (vanilla)
  # ------------------------------------------------------------------
  t0      <- proc.time()
  fit_stm <- tryCatch(
    stm(documents = stm_fmt$docs, vocab = stm_fmt$vocab,
        K = K_TOPICS, prevalence = prevalence_formula, data = cov_df,
        init.type = "Spectral", verbose = FALSE, max.em.its = 75L),
    error = function(e) NULL
  )
  t_stm <- (proc.time() - t0)[3L]

  if (!is.null(fit_stm)) {
    Bz_stm  <- .stm_extract_Bz(fit_stm, dat$C, K_TOPICS)
    Phi_stm <- .stm_extract_Phi(fit_stm, stm_fmt, N_VOCAB)

    if (!is.null(Bz_stm)) {
      pa_stm  <- procrustes_align(Bz_stm, Bz0_r)
      phi_mse_stm <- if (!is.null(Phi_stm))
        best_perm_mse_phi(Phi_stm, dat$Beta, K_TOPICS) else NA_real_
      n_it_stm <- tryCatch(
        as.integer(fit_stm$convergence$its),
        error = function(e) NA_integer_
      )
      metrics <- rbind(metrics, data.frame(
        rep = r, method = "stm",
        mse_Bz = pa_stm$mse, mse_phi = phi_mse_stm,
        time_s = t_stm, n_iter = n_it_stm
      ))
    }
  }

  # ------------------------------------------------------------------
  # Method 4: STM warm-started from sgscatm exact solution
  # ------------------------------------------------------------------
  if (!is.null(fit_sg)) {
    model_init <- .build_stm_init(fit_sg, stm_fmt, K_TOPICS)

    t0        <- proc.time()
    fit_stm_w <- tryCatch(
      stm(documents = stm_fmt$docs, vocab = stm_fmt$vocab,
          K = K_TOPICS, prevalence = prevalence_formula, data = cov_df,
          model = model_init,          # restart path: skips stm.init entirely
          verbose = FALSE, max.em.its = 75L),
      error = function(e) {
        warning(sprintf("rep %d: STM warm start failed — %s", r, e$message))
        NULL
      }
    )
    # Total time = sgscatm (to build the init) + STM EM
    t_warm <- t_sg + (proc.time() - t0)[3L]

    if (!is.null(fit_stm_w)) {
      Bz_w  <- .stm_extract_Bz(fit_stm_w, dat$C, K_TOPICS)
      Phi_w <- .stm_extract_Phi(fit_stm_w, stm_fmt, N_VOCAB)

      if (!is.null(Bz_w)) {
        pa_w        <- procrustes_align(Bz_w, Bz0_r)
        phi_mse_w   <- if (!is.null(Phi_w))
          best_perm_mse_phi(Phi_w, dat$Beta, K_TOPICS) else NA_real_
        n_it_w <- tryCatch(
          as.integer(fit_stm_w$convergence$its),
          error = function(e) NA_integer_
        )
        metrics <- rbind(metrics, data.frame(
          rep = r, method = "stm_warm",
          mse_Bz = pa_w$mse, mse_phi = phi_mse_w,
          time_s = t_warm, n_iter = n_it_w
        ))
      }
    }
  }
}

saveRDS(metrics, out_rds)
cat(sprintf(
  "Worker done: M=%d scenario=%s reps_ok=%d/%d\n",
  M_val, scenario,
  nrow(metrics[metrics$method == "sgscatm", ]), N_REP
))
