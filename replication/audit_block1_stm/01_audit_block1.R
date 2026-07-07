#!/usr/bin/env Rscript
# ===================================================================
#  Block 1 forensic audit — A2 (mechanical reproduction) and
#  A3 (corrected diagnostics with the two-step estimator)
# ===================================================================
#
#  A2 reproduces the published Block 1 pipeline verbatim (same seeds as
#  run_simulation.R Block 1: 10000*M_index + rep) at 20 reps per M and
#  records every quantity needed to explain H1 (RMSE plateau) and H2
#  (coverage 1.000) numerically.
#
#  A3 runs the two-step estimator (pilot -> oracle GL alignment ->
#  k = 5 Gauss-Newton sweeps at lambda = 0 -> OLS of Z on C) with a
#  sandwich covariance for the OLS step, entrywise coverage after
#  Procrustes alignment (covariance rotated consistently), and the
#  alignment-free row-norm coverage.  Seeds 60000 + M_index*1000 + rep.
#
#  CAVEAT logged in REPORT_AUDIT.md: the GL alignment of the start uses
#  Z_true (oracle).  Feasible truth-free starts were tested and FAIL
#  (raw pilot and simplex-projected-Phi starts reach the same objective
#  value with mse_Bz at the zero-estimator level): the exact objective
#  does not identify the orientation; the aligned start supplies it.
#
#  Usage (from the package root):
#    Rscript replication/audit_block1_stm/01_audit_block1.R [--ncores N]
# ===================================================================

args   <- commandArgs(trailingOnly = TRUE)
NCORES <- if ("--ncores" %in% args) {
  as.integer(args[which(args == "--ncores") + 1L])
} else {
  max(1L, min(10L, parallel::detectCores() - 2L))
}

ROOT <- Sys.getenv("SGSCATM_ROOT", unset = normalizePath(".", mustWork = FALSE))
if (!dir.exists(file.path(ROOT, "R")))
  stop("Run from the package root or set SGSCATM_ROOT.")
RES_DIR <- file.path(ROOT, "replication", "audit_block1_stm", "results")
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)

aud_setup <- function(root) {
  source(file.path(root, "R", "sgscatm_fit.R"))
  source(file.path(root, "R", "ilr_contrast.R"))
  source(file.path(root, "R", "utils.R"))
  source(file.path(root, "R", "ilr_se.R"))
  source(file.path(root, "replication", "simulation", "sim_dgp.R"))
  source(file.path(root, "replication", "simulation", "sim_utils.R"))
  source(file.path(root, "replication", "basin_check", "01_functions.R"))
  invisible(NULL)
}
aud_setup(ROOT)

# Block 1 design constants (run_simulation.R lines 82-93)
N_VOCAB  <- 500L; K_TOPICS <- 5L; P_COV <- 3L
Bz0_TRUE <- matrix(c(
   0.40, -0.20,  0.10,  0.30,
  -0.15,  0.35, -0.25,  0.05,
   0.20,  0.10,  0.40, -0.30
), nrow = P_COV, ncol = K_TOPICS - 1L, byrow = TRUE)
M_VALUES <- c(500L, 1000L, 2000L)

# ===================================================================
#  Sandwich covariance for multivariate OLS  B_hat = (C'C)^{-1} C'Z
#
#  vec(B_hat) - vec(B) = (I_{K-1} ox (C'C)^{-1}) vec(C'E), so with
#  plug-in row residuals r_i,
#    Sigma = (I ox (C'C)^{-1}) [ sum_i (r_i r_i') ox (c_i c_i') ]
#            (I ox (C'C)^{-1}) * M/(M-P)      (HC1)
#  Meat blocks: block (k,k') = sum_i r_ik r_ik' c_i c_i'
#             = crossprod(C, C * (R[,k]*R[,k'])).
#  vec is column-stacking: entry (j,k) of B <-> index (k-1)P + j.
# ===================================================================
aud_sandwich <- function(Z, C, B_hat) {
  M <- nrow(C); P <- ncol(C); Kp <- ncol(Z)
  R    <- Z - C %*% B_hat
  XtXi <- solve(crossprod(C))
  Meat <- matrix(0, P * Kp, P * Kp)
  for (k in seq_len(Kp)) for (kp in seq_len(Kp)) {
    blk <- crossprod(C, C * (R[, k] * R[, kp]))
    Meat[(k - 1L) * P + seq_len(P), (kp - 1L) * P + seq_len(P)] <- blk
  }
  Bread <- kronecker(diag(Kp), XtXi)
  (M / (M - P)) * (Bread %*% Meat %*% Bread)
}

#' Unit test: sandwich vs direct simulation of heteroscedastic OLS.
#' Hard-stops if the mean sandwich estimate misses the empirical
#' covariance of vec(B_hat) by more than 25% in Frobenius norm.
aud_test_sandwich <- function(seed = 7001L, n_sim = 3000L) {
  set.seed(seed)
  M <- 200L; P <- 2L; Kp <- 3L
  C <- cbind(rnorm(M), runif(M, -1, 1))
  B <- matrix(rnorm(P * Kp), P, Kp)
  # heteroscedastic, row-correlated noise: sd depends on |c_i1|
  sdv <- 0.5 * (1 + abs(C[, 1]))
  A   <- matrix(c(1, .3, 0, .3, 1, .2, 0, .2, 1), 3, 3)  # row mixing
  vecs <- matrix(NA_real_, n_sim, P * Kp)
  Sig_acc <- matrix(0, P * Kp, P * Kp)
  for (s in seq_len(n_sim)) {
    E <- (matrix(rnorm(M * Kp), M, Kp) %*% A) * sdv
    Z <- C %*% B + E
    Bh <- solve(crossprod(C), crossprod(C, Z))
    vecs[s, ] <- as.vector(Bh)
    Sig_acc <- Sig_acc + aud_sandwich(Z, C, Bh)
  }
  emp <- cov(vecs)
  avg <- Sig_acc / n_sim
  rel <- norm(avg - emp, "F") / norm(emp, "F")
  if (rel > 0.25)
    stop(sprintf("Sandwich unit test FAILED: rel Frobenius error %.3f", rel))
  rel
}

# ===================================================================
#  A2: mechanical reproduction of the published Block 1 pipeline
# ===================================================================
run_a2_replicate <- function(M, M_idx, rep_id, dummy_cfg) {
  # published seed scheme: run_simulation.R line 128
  seed <- 10000L * M_idx + rep_id
  dat <- sim_dgp(M = M, N = N_VOCAB, K = K_TOPICS, P = P_COV,
                 Bz0 = Bz0_TRUE, sigma_eps = 0.3, alpha_beta = 0.1,
                 doc_length = 200L, seed = seed)
  fit <- sgscatm(dat$W, dat$C, K = K_TOPICS, lambda = 1, rotate = TRUE)
  pa  <- procrustes_align(fit$Bz, dat$Bz0)

  se_res <- tryCatch(ilr_se_analytical(fit), error = function(e) NULL)
  out <- list(M = M, rep = rep_id, seed = seed,
              mse = pa$mse,
              norm_ratio = sqrt(sum(fit$Bz^2)) / sqrt(sum(dat$Bz0^2)),
              center_offset = sqrt(mean((pa$Bz_aligned - dat$Bz0)^2)))
  if (!is.null(se_res)) {
    z_q <- qnorm(0.975)
    # exactly as published (run_simulation.R lines 152-158):
    se_rot <- .rotate_se(se_res, pa$R, P_COV, K_TOPICS - 1L)
    covers <- (dat$Bz0 >= pa$Bz_aligned - z_q * se_rot) &
              (dat$Bz0 <= pa$Bz_aligned + z_q * se_rot)
    # counterfactual: same chain but with the 1/sqrt(M) factor that
    # ilr_se_analytical applies to its own `se` slot
    se_cor <- se_rot / sqrt(M)
    covers_cor <- (dat$Bz0 >= pa$Bz_aligned - z_q * se_cor) &
                  (dat$Bz0 <= pa$Bz_aligned + z_q * se_cor)
    out$coverage_published <- mean(covers)
    out$coverage_scaledSE  <- mean(covers_cor)
    out$med_se_published   <- median(se_rot)
    out$med_se_scaled      <- median(se_cor)
    out$med_abs_offset     <- median(abs(pa$Bz_aligned - dat$Bz0))
  }
  out
}

# ===================================================================
#  A3: corrected diagnostics with the two-step estimator
# ===================================================================
run_a3_replicate <- function(M, M_idx, rep_id, dummy_cfg) {
  seed <- 60000L + M_idx * 1000L + rep_id
  dat <- sim_dgp(M = M, N = N_VOCAB, K = K_TOPICS, P = P_COV,
                 Bz0 = Bz0_TRUE, sigma_eps = 0.3, alpha_beta = 0.1,
                 doc_length = 200L, seed = seed)
  Wf <- dat$W / rowSums(dat$W)
  V  <- dat$V

  t0  <- proc.time()
  fit <- sgscatm(dat$W, dat$C, K = K_TOPICS, lambda = 1, rotate = TRUE)
  pilot_mse <- procrustes_align(fit$Bz, dat$Bz0)$mse   # published convention

  gl <- bc_gl_align(fit$Z, dat$Z_true)                 # ORACLE alignment
  rf <- bc_refine(gl$Z, bc_phi_step(gl$Z, Wf, V), Wf, dat$C, V,
                  lambda = 0, max_sweeps = 5L)         # fixed k = 5
  B_hat <- bc_b_step(rf$Z, dat$C)
  t_all <- (proc.time() - t0)[3]

  # Procrustes alignment + consistently rotated sandwich covariance
  pa    <- procrustes_align(B_hat, dat$Bz0)
  Sigma <- aud_sandwich(rf$Z, dat$C, B_hat)
  Kp    <- K_TOPICS - 1L
  RkI   <- kronecker(t(pa$R), diag(P_COV))
  Sigma_rot <- RkI %*% Sigma %*% t(RkI)
  se    <- matrix(sqrt(pmax(diag(Sigma_rot), 0)), P_COV, Kp)
  z_q   <- qnorm(0.975)
  covers <- (dat$Bz0 >= pa$Bz_aligned - z_q * se) &
            (dat$Bz0 <= pa$Bz_aligned + z_q * se)

  # alignment-free: squared row norms s_j (right-orthogonal invariant)
  s_hat <- rowSums(B_hat^2)
  s_se  <- vapply(seq_len(P_COV), function(j) {
    idx <- j + (seq_len(Kp) - 1L) * P_COV
    sqrt(max(4 * B_hat[j, ] %*% Sigma[idx, idx] %*% B_hat[j, ], 0))
  }, numeric(1))
  s_true   <- rowSums(dat$Bz0^2)
  s_covers <- (s_true >= s_hat - z_q * s_se) & (s_true <= s_hat + z_q * s_se)

  list(M = M, rep = rep_id, seed = seed,
       pilot_mse = pilot_mse,
       mse = pa$mse,                       # corrected two-step, paper metric
       mse_direct = mean((B_hat - dat$Bz0)^2),
       B_tilde = pa$Bz_aligned, se_mat = se,   # for bias/variance split
       norm_ratio = sqrt(sum(B_hat^2)) / sqrt(sum(dat$Bz0^2)),
       coverage = mean(covers),
       cover_mat = covers,
       std_err_mat = (pa$Bz_aligned - dat$Bz0) / se,   # for the QQ plot
       s_hat = s_hat, s_se = s_se, s_covers = s_covers,
       n_fail = sum(rf$trace$n_fail), monotone_ok = !rf$monotone_violation,
       time_s = t_all)
}

# ===================================================================
#  Driver (PSOCK; mclapply is serial on Windows)
# ===================================================================
aud_run_jobs <- function(jobs, fun, ncores,
                         export = c("SIGNAL_LEVELS_DUMMY")) {
  if (ncores <= 1L || length(jobs) <= 1L)
    return(lapply(jobs, function(j) do.call(fun, j)))
  cl <- parallel::makePSOCKcluster(min(ncores, length(jobs)))
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(cl, aud_setup, root = ROOT)
  parallel::clusterExport(cl, c("N_VOCAB", "K_TOPICS", "P_COV", "Bz0_TRUE",
                                "aud_sandwich", "run_a2_replicate",
                                "run_a3_replicate"),
                          envir = globalenv())
  wf <- function(j) do.call(FUN, j)
  environment(wf) <- list2env(list(FUN = fun), parent = globalenv())
  parallel::parLapplyLB(cl, jobs, wf)
}

wrap_try <- function(fun) {
  force(fun)
  function(...) tryCatch(fun(...), error = function(e)
    list(error = conditionMessage(e), args = list(...)))
}

cat("=== Block 1 forensic audit ===\n")
cat("-- Sandwich unit test (M=200, P=2, K-1=3, heteroscedastic) --\n")
t0 <- proc.time()
rel <- aud_test_sandwich()
cat(sprintf("   rel Frobenius error (mean sandwich vs empirical cov): %.3f",
            rel), "\n   PASSED\n")

cat("\n-- A2: published pipeline, 20 reps x M in {500,1000,2000} --\n")
t0 <- proc.time()
jobs <- list()
for (mi in seq_along(M_VALUES))
  for (r in seq_len(20L))
    jobs[[length(jobs) + 1L]] <- list(M = M_VALUES[mi], M_idx = mi,
                                      rep_id = r, dummy_cfg = NULL)
a2 <- aud_run_jobs(jobs, wrap_try(run_a2_replicate), NCORES)
t_a2 <- (proc.time() - t0)[3]
saveRDS(list(results = a2, Bz0 = Bz0_TRUE, time_s = t_a2),
        file.path(RES_DIR, "a2_results.rds"))
cat(sprintf("   [%.1f s] %d/%d replicates OK\n", t_a2,
            sum(vapply(a2, function(x) is.null(x$error), logical(1))),
            length(a2)))

cat("\n-- A3: two-step estimator, 50 reps x M in {500,1000,2000} --\n")
t0 <- proc.time()
jobs <- list()
for (mi in seq_along(M_VALUES))
  for (r in seq_len(50L))
    jobs[[length(jobs) + 1L]] <- list(M = M_VALUES[mi], M_idx = mi,
                                      rep_id = r, dummy_cfg = NULL)
a3 <- aud_run_jobs(jobs, wrap_try(run_a3_replicate), NCORES)
t_a3 <- (proc.time() - t0)[3]
saveRDS(list(results = a3, Bz0 = Bz0_TRUE, time_s = t_a3),
        file.path(RES_DIR, "a3_results.rds"))
cat(sprintf("   [%.1f s] %d/%d replicates OK\n", t_a3,
            sum(vapply(a3, function(x) is.null(x$error), logical(1))),
            length(a3)))
cat("Done. Results in", RES_DIR, "\n")
