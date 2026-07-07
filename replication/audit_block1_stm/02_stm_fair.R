#!/usr/bin/env Rscript
# ===================================================================
#  Fair Block 3 STM comparison — B1 (ALR -> ILR map), B2 (paired
#  rerun on the basin_check replicates), B3 metrics
# ===================================================================
#
#  B1: STM parametrises prevalence in ALR coordinates with reference
#  topic K (eta_K = 0): theta = softmax([eta; 0]), E[eta] = Gamma' x.
#  With L = log(theta), alr = L_{1..K-1} - L_K and ilr = V'L; since
#  V'1 = 0, ilr = V'[alr; 0] = M_map %*% alr with
#      M_map = t(V[1:(K-1), ])          ((K-1) x (K-1)).
#  Row form for score matrices: ILR = ALR %*% t(M_map), hence for the
#  coefficient matrices (dropping STM's intercept row):
#      Gamma_ilr = Gamma_alr %*% t(M_map).
#  Verified to machine precision on synthetic compositions before any
#  STM output is touched (hard stop otherwise), together with the
#  inverse relation alr = A_map %*% ilr, A_map = V[1:(K-1),] - 1 V[K,],
#  M_map %*% A_map = I, and the DGP round-trip Bz0 -> Gamma_alr -> Bz0.
#
#  B2: same data as basin_check E1/E2 (seeds 90000 + regime*1000 + rep,
#  weak = 1, strong = 2; Bz0 drawn inside sim_dgp under the seed).
#  STM fit IN-SESSION (R 4.5.1 native), spectral init, prevalence ~ C,
#  max.em.its = 75 — the old worker's settings.  Two STM extractions:
#    STM-gamma : native prevalence coefficients mapped ALR -> ILR (B1)
#    STM-theta : the old worker's functional (log theta -> ILR -> OLS)
#
#  Usage: Rscript replication/audit_block1_stm/02_stm_fair.R [--ncores N]
# ===================================================================

args   <- commandArgs(trailingOnly = TRUE)
NCORES <- if ("--ncores" %in% args) {
  as.integer(args[which(args == "--ncores") + 1L])
} else {
  max(1L, min(10L, parallel::detectCores() - 2L))
}

ROOT <- Sys.getenv("SGSCATM_ROOT", unset = normalizePath(".", mustWork = FALSE))
RES_DIR <- file.path(ROOT, "replication", "audit_block1_stm", "results")
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)

stm_setup <- function(root) {
  source(file.path(root, "R", "sgscatm_fit.R"))
  source(file.path(root, "R", "ilr_contrast.R"))
  source(file.path(root, "R", "utils.R"))
  source(file.path(root, "replication", "simulation", "sim_dgp.R"))
  source(file.path(root, "replication", "simulation", "sim_utils.R"))
  source(file.path(root, "replication", "basin_check", "01_functions.R"))
  suppressPackageStartupMessages(library(stm))
  invisible(NULL)
}
stm_setup(ROOT)

SIGNAL_LEVELS <- c(weak = 0.15, strong = 0.50)
K_TOPICS <- 5L; P_COV <- 3L; N_VOCAB <- 500L

# ===================================================================
#  B1: coordinate map + hard-stop unit test
# ===================================================================
stm_alr_map <- function(V) t(V[seq_len(nrow(V) - 1L), , drop = FALSE])

stm_test_map <- function(K = 5L, seed = 7101L, tol = 1e-12) {
  set.seed(seed)
  V <- ilr_contrast(K)
  theta <- exp(matrix(rnorm(200 * K, 0, 2), 200, K))
  theta <- theta / rowSums(theta)
  alr <- log(theta[, -K, drop = FALSE] / theta[, K])
  ilr <- log(theta) %*% V
  M_map <- stm_alr_map(V)
  e1 <- max(abs(ilr - alr %*% t(M_map)))
  A_map <- V[seq_len(K - 1L), , drop = FALSE] -
    matrix(V[K, ], K - 1L, K - 1L, byrow = TRUE)
  e2 <- max(abs(alr - ilr %*% t(A_map)))
  e3 <- max(abs(M_map %*% A_map - diag(K - 1L)))
  # DGP round trip: z = Bz0' c  =>  Gamma_alr = Bz0 %*% t(A_map), and
  # mapping back must recover Bz0 exactly
  Bz0 <- matrix(rnorm(3 * (K - 1L)), 3, K - 1L)
  e4 <- max(abs(Bz0 %*% t(A_map) %*% t(M_map) - Bz0))
  errs <- c(ilr_from_alr = e1, alr_from_ilr = e2, inverse = e3,
            roundtrip = e4)
  if (max(errs) > tol)
    stop("ALR->ILR map unit test FAILED: ",
         paste(sprintf("%s=%.2e", names(errs), errs), collapse = ", "))
  errs
}

# ===================================================================
#  B2 replicate: pilot | refined | STM (two extractions), paired data
# ===================================================================
run_b2_replicate <- function(regime, rep_id, M) {
  regime_idx <- match(regime, names(SIGNAL_LEVELS))
  seed <- 90000L + regime_idx * 1000L + rep_id       # basin_check scheme
  dat <- sim_dgp(M = M, N = N_VOCAB, K = K_TOPICS, P = P_COV,
                 b_max = SIGNAL_LEVELS[[regime]], sigma_eps = 0.3,
                 alpha_beta = 0.1, doc_length = 200L, seed = seed)
  Wf <- dat$W / rowSums(dat$W)
  V  <- dat$V
  M_map <- stm_alr_map(V)
  nB0 <- sqrt(sum(dat$Bz0^2))

  out <- list(regime = regime, rep = rep_id, M = M, seed = seed)

  # --- (1) pilot, published pipeline --------------------------------
  t0 <- proc.time()
  fit <- sgscatm(dat$W, dat$C, K = K_TOPICS, lambda = 1, rotate = TRUE)
  out$pilot <- list(
    mse_paper  = procrustes_align(fit$Bz, dat$Bz0)$mse,
    norm_ratio = sqrt(sum(fit$Bz^2)) / nB0,
    mse_gl     = bc_mse_direct(bc_gl_align(fit$Z, dat$Z_true)$Z, dat$C,
                               dat$Bz0),                 # oracle column
    time_s     = (proc.time() - t0)[3])

  # --- (2) pilot + refined (k = 5, lambda = 0) -----------------------
  t0 <- proc.time()
  gl <- bc_gl_align(fit$Z, dat$Z_true)                   # oracle start
  rf <- bc_refine(gl$Z, bc_phi_step(gl$Z, Wf, V), Wf, dat$C, V,
                  lambda = 0, max_sweeps = 5L)
  B_ref <- bc_b_step(rf$Z, dat$C)
  out$refined <- list(
    mse_paper  = procrustes_align(B_ref, dat$Bz0)$mse,
    mse_gl     = mean((B_ref - dat$Bz0)^2),              # model coords
    norm_ratio = sqrt(sum(B_ref^2)) / nB0,
    monotone_ok = !rf$monotone_violation,
    time_s     = out$pilot$time_s + (proc.time() - t0)[3])

  # --- (3) STM, in-session, old worker's settings -------------------
  active_cols <- which(colSums(dat$W) > 0L)
  W_stm <- dat$W[, active_cols, drop = FALSE]
  stm_docs <- lapply(seq_len(nrow(W_stm)), function(i) {
    idx <- which(W_stm[i, ] > 0L)
    rbind(as.integer(idx), as.integer(W_stm[i, idx]))
  })
  cov_df <- as.data.frame(dat$C)
  colnames(cov_df) <- paste0("V", seq_len(P_COV))
  t0 <- proc.time()
  fit_stm <- tryCatch(
    stm::stm(documents = stm_docs, vocab = as.character(active_cols),
             K = K_TOPICS,
             prevalence = ~ V1 + V2 + V3, data = cov_df,
             init.type = "Spectral", verbose = FALSE, max.em.its = 75L),
    error = function(e) NULL)
  t_stm <- (proc.time() - t0)[3]

  if (is.null(fit_stm)) { out$stm_error <- TRUE; return(out) }

  em_its <- tryCatch(length(fit_stm$convergence$bound),
                     error = function(e) NA_integer_)

  # STM-gamma: native prevalence coefficients, ALR -> ILR
  Gamma <- fit_stm$mu$gamma                              # (P+1) x (K-1)
  B_stm_g <- Gamma[-1L, , drop = FALSE] %*% t(M_map)     # P x (K-1), ILR
  out$stm_gamma <- list(
    mse_paper  = procrustes_align(B_stm_g, dat$Bz0)$mse,
    norm_ratio = sqrt(sum(B_stm_g^2)) / nB0,
    time_s     = t_stm, em_its = em_its)

  # STM-theta: the old worker's functional (block3_stm_worker.R 74-80)
  Z_stm <- log(pmax(fit_stm$theta, 1e-10)) %*% V
  B_stm_t <- solve(crossprod(dat$C), crossprod(dat$C, Z_stm))
  out$stm_theta <- list(
    mse_paper  = procrustes_align(B_stm_t, dat$Bz0)$mse,
    norm_ratio = sqrt(sum(B_stm_t^2)) / nB0,
    mse_gl     = bc_mse_direct(bc_gl_align(Z_stm, dat$Z_true)$Z, dat$C,
                               dat$Bz0),                 # oracle column
    time_s     = t_stm, em_its = em_its)
  out
}

# ===================================================================
#  Driver
# ===================================================================
wrap_try <- function(fun) {
  force(fun)
  function(...) tryCatch(fun(...), error = function(e)
    list(error = conditionMessage(e), args = list(...)))
}
stm_run_jobs <- function(jobs, fun, ncores) {
  if (ncores <= 1L || length(jobs) <= 1L)
    return(lapply(jobs, function(j) do.call(fun, j)))
  cl <- parallel::makePSOCKcluster(min(ncores, length(jobs)))
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(cl, stm_setup, root = ROOT)
  parallel::clusterExport(cl, c("SIGNAL_LEVELS", "K_TOPICS", "P_COV",
                                "N_VOCAB", "stm_alr_map",
                                "run_b2_replicate"),
                          envir = globalenv())
  wf <- function(j) do.call(FUN, j)
  environment(wf) <- list2env(list(FUN = fun), parent = globalenv())
  parallel::parLapplyLB(cl, jobs, wf)
}

cat("=== Fair STM comparison ===\n")
cat("-- B1: ALR -> ILR map unit test --\n")
errs <- stm_test_map()
cat(sprintf("   max errors: %s\n   PASSED (tol 1e-12)\n",
            paste(sprintf("%s=%.1e", names(errs), errs), collapse = " ")))
saveRDS(errs, file.path(RES_DIR, "b1_map_test.rds"))

cat("\n-- B2: M = 1000, 20 reps x 2 regimes (paired with basin_check) --\n")
t0 <- proc.time()
jobs <- list()
for (rg in names(SIGNAL_LEVELS))
  for (r in seq_len(20L))
    jobs[[length(jobs) + 1L]] <- list(regime = rg, rep_id = r, M = 1000L)
b2 <- stm_run_jobs(jobs, wrap_try(run_b2_replicate), NCORES)
t_b2 <- (proc.time() - t0)[3]
saveRDS(list(results = b2, time_s = t_b2),
        file.path(RES_DIR, "b2_results.rds"))
n_ok <- sum(vapply(b2, function(x)
  is.null(x$error) && is.null(x$stm_error), logical(1)))
cat(sprintf("   [%.1f s] %d/%d replicates fully OK\n", t_b2, n_ok, length(b2)))

# optional M = 5000 strong if the M = 1000 block was quick enough
stm_times <- unlist(lapply(b2, function(x)
  if (!is.null(x$stm_gamma)) x$stm_gamma$time_s else NULL))
if (length(stm_times) && median(stm_times) < 120) {
  cat("\n-- B2+: M = 5000, strong, 10 reps (optional) --\n")
  t0 <- proc.time()
  jobs <- lapply(seq_len(10L), function(r)
    list(regime = "strong", rep_id = r, M = 5000L))
  b2m5 <- stm_run_jobs(jobs, wrap_try(run_b2_replicate), NCORES)
  t_m5 <- (proc.time() - t0)[3]
  saveRDS(list(results = b2m5, time_s = t_m5),
          file.path(RES_DIR, "b2_m5000_results.rds"))
  cat(sprintf("   [%.1f s] %d/%d replicates OK\n", t_m5,
              sum(vapply(b2m5, function(x)
                is.null(x$error) && is.null(x$stm_error), logical(1))),
              length(b2m5)))
} else {
  cat("\n   Skipping optional M = 5000 (median STM time ",
      sprintf("%.0f", median(stm_times)), "s >= 120s)\n", sep = "")
}
cat("Done. Results in", RES_DIR, "\n")
