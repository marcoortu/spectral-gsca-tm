#!/usr/bin/env Rscript
# ===================================================================
#  Deviance probe — runner for P1 (pathology probe, the decision
#  experiment), P2 (k-curves + inference, conditional), P3 (feasible
#  chain under deviance)
# ===================================================================
#
#  Usage (from the package root):
#    Rscript replication/deviance_probe/02_run.R [--quick] [--ncores N]
#
#  Decision flow: unit tests (hard gates) -> P1 -> gate evaluated
#  in-script -> P2a always, P2b only if gate P1 passes -> P3.
#
#  Seeds: P1/P2a/P3 use the E2 scheme 90000 + regime*1000 + rep
#  (weak = 1, strong = 2; P1 reps 1-5 = the F3/D1 probe data; P3
#  crossover: regime strong, M = 5000, reps 1-10, paired with audit
#  B2+).  P2b: 60000 + M_index*1000 + rep (paired with feasibility F2).
#  Unit tests: 66001.
# ===================================================================

args   <- commandArgs(trailingOnly = TRUE)
QUICK  <- "--quick" %in% args
NCORES <- if ("--ncores" %in% args) {
  as.integer(args[which(args == "--ncores") + 1L])
} else {
  max(1L, min(10L, parallel::detectCores() - 2L))
}

ROOT <- Sys.getenv("SGSCATM_ROOT", unset = normalizePath(".", mustWork = FALSE))
if (!dir.exists(file.path(ROOT, "R")))
  stop("Run from the package root or set SGSCATM_ROOT.")
RES_DIR <- file.path(ROOT, "replication", "deviance_probe", "results")
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)

dv_setup <- function(root) {
  source(file.path(root, "R", "sgscatm_fit.R"))
  source(file.path(root, "R", "ilr_contrast.R"))
  source(file.path(root, "R", "utils.R"))
  source(file.path(root, "replication", "simulation", "sim_dgp.R"))
  source(file.path(root, "replication", "simulation", "sim_utils.R"))
  source(file.path(root, "replication", "basin_check", "01_functions.R"))
  source(file.path(root, "replication", "feasibility", "01_anchors.R"))
  source(file.path(root, "replication", "feasibility",
                   "02_constrained_refine.R"))
  source(file.path(root, "replication", "feasibility", "03_jackknife.R"))
  source(file.path(root, "replication", "deviance_probe",
                   "01_deviance_blocks.R"))
  invisible(NULL)
}
dv_setup(ROOT)

SIGNAL_LEVELS <- c(weak = 0.15, strong = 0.50)
K_TOPICS <- 5L; P_COV <- 3L; N_VOCAB <- 500L
Bz0_TRUE <- matrix(c(
   0.40, -0.20,  0.10,  0.30,
  -0.15,  0.35, -0.25,  0.05,
   0.20,  0.10,  0.40, -0.30
), nrow = P_COV, ncol = K_TOPICS - 1L, byrow = TRUE)

N_REP_P1  <- if (QUICK) 2L else 5L
N_REP_E2  <- if (QUICK) 2L else 20L
N_REP_B1  <- if (QUICK) 3L else 50L
M_GRID    <- if (QUICK) c(500L, 1000L) else c(500L, 1000L, 2000L, 4000L)
CAP       <- 100L

#' oracle-GL pilot start in a deviance-compatible form (Phi rows
#' clipped to 1e-8 and renormalised; the EM step re-grows them)
dv_oracle_start <- function(dat, Wf) {
  fit <- sgscatm(dat$W, dat$C, K = K_TOPICS, lambda = 1, rotate = TRUE)
  gl  <- bc_gl_align(fit$Z, dat$Z_true)
  Phi0 <- bc_phi_step(gl$Z, Wf, dat$V)
  Phi0 <- pmax(Phi0, 1e-8)
  list(Z = gl$Z, Phi = Phi0 / rowSums(Phi0), fit = fit)
}

# ===================================================================
#  Replicate drivers
# ===================================================================

#' P1: truth-start probes (dev + LS, gauge-tracked) + dev oracle start
run_p1_replicate <- function(regime, rep_id) {
  seed <- 90000L + match(regime, names(SIGNAL_LEVELS)) * 1000L + rep_id
  dat <- sim_dgp(M = 1000L, N = N_VOCAB, K = K_TOPICS, P = P_COV,
                 b_max = SIGNAL_LEVELS[[regime]], sigma_eps = 0.3,
                 alpha_beta = 0.1, doc_length = 200L, seed = seed)
  Wn <- dat$W; Wf <- Wn / rowSums(Wn); V <- dat$V
  QT   <- bc_gauge_basis(dat$Z_true, dat$Beta, V)
  eta0 <- bc_pack(dat$Z_true, dat$Beta)
  mse0 <- procrustes_align(bc_b_step(dat$Z_true, dat$C), dat$Bz0)$mse

  t0 <- proc.time()
  r_dev <- dv_refine(dat$Z_true, dat$Beta, Wn, dat$C, V,
                     criterion = "dev", max_sweeps = CAP,
                     apply_rule = FALSE, Bz0 = dat$Bz0,
                     QT = QT, eta0 = eta0)
  t_dev <- (proc.time() - t0)[3]
  r_ls  <- dv_refine(dat$Z_true, dat$Beta, Wn, dat$C, V,
                     criterion = "ls", max_sweeps = CAP,
                     apply_rule = FALSE, Bz0 = dat$Bz0,
                     QT = QT, eta0 = eta0)

  os <- dv_oracle_start(dat, Wf)
  r_dvo <- dv_refine(os$Z, os$Phi, Wn, dat$C, V, criterion = "dev",
                     max_sweeps = CAP, apply_rule = FALSE,
                     Bz0 = dat$Bz0, QT = QT, eta0 = eta0)

  keep <- function(r) r[c("F_path", "mse_path", "nr_path", "gauge_path",
                          "perp_path", "monotone_ok", "rule_stop",
                          "sweeps")]
  list(regime = regime, rep = rep_id, seed = seed, mse0 = mse0,
       dev_truth = keep(r_dev), ls_truth = keep(r_ls),
       dev_oracle = keep(r_dvo), time_dev_s = t_dev)
}

#' P2a: deviance k-curves from the oracle-GL start
run_p2a_replicate <- function(regime, rep_id) {
  seed <- 90000L + match(regime, names(SIGNAL_LEVELS)) * 1000L + rep_id
  dat <- sim_dgp(M = 1000L, N = N_VOCAB, K = K_TOPICS, P = P_COV,
                 b_max = SIGNAL_LEVELS[[regime]], sigma_eps = 0.3,
                 alpha_beta = 0.1, doc_length = 200L, seed = seed)
  Wn <- dat$W; Wf <- Wn / rowSums(Wn)
  os <- dv_oracle_start(dat, Wf)
  r <- dv_refine(os$Z, os$Phi, Wn, dat$C, dat$V, criterion = "dev",
                 max_sweeps = CAP, apply_rule = FALSE, Bz0 = dat$Bz0,
                 track_B_path = TRUE)
  list(regime = regime, rep = rep_id, seed = seed,
       mse_path = r$mse_path, nr_path = r$nr_path,
       rule_stop = fs_rule_from_path(r$B_path, cap = CAP),
       monotone_ok = r$monotone_ok)
}

#' P2b: Block 1 grid — deviance-refined (rule), sandwich coverage
run_p2b_replicate <- function(M, M_idx, rep_id) {
  seed <- 60000L + M_idx * 1000L + rep_id
  dat <- sim_dgp(M = M, N = N_VOCAB, K = K_TOPICS, P = P_COV,
                 Bz0 = Bz0_TRUE, sigma_eps = 0.3, alpha_beta = 0.1,
                 doc_length = 200L, seed = seed)
  Wn <- dat$W; Wf <- Wn / rowSums(Wn)
  t0 <- proc.time()
  os <- dv_oracle_start(dat, Wf)
  r <- dv_refine(os$Z, os$Phi, Wn, dat$C, dat$V, criterion = "dev",
                 max_sweeps = CAP, apply_rule = TRUE)
  Sig <- fs_sandwich(r$Z, dat$C, r$B)
  Kp <- K_TOPICS - 1L
  cv <- fs_coverage_entry(r$B, dat$Bz0, Sig, P_COV, Kp)
  list(M = M, rep = rep_id, seed = seed,
       mse = cv$mse, B_tilde = procrustes_align(r$B, dat$Bz0)$Bz_aligned,
       coverage = mean(cv$covers),
       rn_covers = fs_coverage_rownorm(r$B, dat$Bz0, Sig, P_COV, Kp),
       sweeps = r$sweeps, rule_stop = r$rule_stop,
       monotone_ok = r$monotone_ok, time_s = (proc.time() - t0)[3])
}

#' P3: anchor polish + feasible deviance chain
run_p3_replicate <- function(regime, rep_id, M = 1000L) {
  seed <- 90000L + match(regime, names(SIGNAL_LEVELS)) * 1000L + rep_id
  dat <- sim_dgp(M = M, N = N_VOCAB, K = K_TOPICS, P = P_COV,
                 b_max = SIGNAL_LEVELS[[regime]], sigma_eps = 0.3,
                 alpha_beta = 0.1, doc_length = 200L, seed = seed)
  Wn <- dat$W; V <- dat$V
  nB0 <- sqrt(sum(dat$Bz0^2))

  t0 <- proc.time()
  ap <- fs_anchor_pipeline(Wn, K_TOPICS)
  tv_anchor <- fs_perm_tv(ap$Phi, dat$Beta)

  # polish: likelihood-based anchor recovery (deviance sweeps from the
  # anchored start)
  zi <- dv_z_init(ap$Phi, Wn, V, n_gn = 10L)
  rp <- dv_refine(zi$Z, ap$Phi, Wn, dat$C, V, criterion = "dev",
                  max_sweeps = CAP, apply_rule = TRUE)
  tv_polish <- fs_perm_tv(rp$Phi, dat$Beta)

  # chain: V4-style pilot orientation onto the POLISHED frame
  fit <- sgscatm(Wn, dat$C, K = K_TOPICS, lambda = 1, rotate = TRUE)
  a0 <- colMeans(rp$Z)
  A_hat <- solve(crossprod(fit$Z) + 1e-8 * diag(K_TOPICS - 1L),
                 crossprod(fit$Z, sweep(rp$Z, 2L, a0)))
  Z0 <- sweep(fit$Z %*% A_hat, 2L, a0, "+")
  rc <- dv_refine(Z0, rp$Phi, Wn, dat$C, V, criterion = "dev",
                  max_sweeps = CAP, apply_rule = TRUE, Bz0 = dat$Bz0)
  t_all <- (proc.time() - t0)[3]

  Sig <- fs_sandwich(rc$Z, dat$C, rc$B)
  Kp <- K_TOPICS - 1L
  cv <- fs_coverage_entry(rc$B, dat$Bz0, Sig, P_COV, Kp)
  list(regime = regime, rep = rep_id, M = M, seed = seed,
       tv_anchor = tv_anchor, tv_polish = tv_polish,
       mse_paper = cv$mse,
       mse_perm = fs_perm_mse(rc$B, dat$Bz0, V),
       norm_ratio = sqrt(sum(rc$B^2)) / nB0,
       coverage = mean(cv$covers),
       sweeps_polish = rp$sweeps, sweeps_chain = rc$sweeps,
       monotone_ok = rp$monotone_ok && rc$monotone_ok,
       time_s = t_all)
}

# ===================================================================
#  Parallel helper
# ===================================================================
wrap_try <- function(fun) {
  force(fun)
  function(...) tryCatch(fun(...), error = function(e)
    list(error = conditionMessage(e), args = list(...)))
}
dv_run_jobs <- function(jobs, fun, ncores) {
  if (ncores <= 1L || length(jobs) <= 1L)
    return(lapply(jobs, function(j) do.call(fun, j)))
  cl <- parallel::makePSOCKcluster(min(ncores, length(jobs)))
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(cl, dv_setup, root = ROOT)
  parallel::clusterExport(cl, c("SIGNAL_LEVELS", "K_TOPICS", "P_COV",
                                "N_VOCAB", "Bz0_TRUE", "CAP",
                                "dv_oracle_start", "run_p1_replicate",
                                "run_p2a_replicate", "run_p2b_replicate",
                                "run_p3_replicate"),
                          envir = globalenv())
  wf <- function(j) do.call(FUN, j)
  environment(wf) <- list2env(list(FUN = fun), parent = globalenv())
  parallel::parLapplyLB(cl, jobs, wf)
}
report_ok <- function(res, t_el) {
  n_ok <- sum(vapply(res, function(x) is.null(x$error), logical(1)))
  cat(sprintf("   [%.1f s] %d/%d replicates OK\n", t_el, n_ok,
              length(res)))
  if (n_ok < length(res))
    cat("   ERRORS:", paste(head(unique(unlist(lapply(res, function(x)
      x$error))), 3), collapse = " | "), "\n")
}

cat("=== Deviance probe ===\n")
cat("Mode:", if (QUICK) "QUICK" else "FULL", "| cores:", NCORES, "\n")

# ===================================================================
#  Unit tests (hard gates)
# ===================================================================
cat("\n-- Deviance block unit tests (M=30, N=40, K=3) --\n")
t0 <- proc.time()
ut <- dv_verify()
cat(sprintf("   z-grad rel.err %.2e | EM monotone %d/20 | sweeps monotone %d/5\n",
            ut$grad_relerr, ut$em_ok, ut$sweep_ok))
saveRDS(ut, file.path(RES_DIR, "unit_tests.rds"))
cat(sprintf("   [%.1f s] PASSED\n", (proc.time() - t0)[3]))

# ===================================================================
#  P1: the decision experiment
# ===================================================================
cat("\n-- P1: pathology probe (truth start, dev vs LS, gauge-tracked) --\n")
t0 <- proc.time()
jobs <- list()
for (rg in names(SIGNAL_LEVELS))
  for (r in seq_len(N_REP_P1))
    jobs[[length(jobs) + 1L]] <- list(regime = rg, rep_id = r)
p1 <- dv_run_jobs(jobs, wrap_try(run_p1_replicate), NCORES)
t_p1 <- (proc.time() - t0)[3]
saveRDS(list(results = p1, time_s = t_p1),
        file.path(RES_DIR, "p1_results.rds"))
report_ok(p1, t_p1)

ok1 <- Filter(function(x) is.null(x$error), p1)
gate_p1 <- all(vapply(ok1, function(x) {
  x$dev_truth$mse_path[length(x$dev_truth$mse_path)] <= 2 * x$mse0
}, logical(1)))
rat <- vapply(ok1, function(x)
  x$dev_truth$mse_path[length(x$dev_truth$mse_path)] / x$mse0, numeric(1))
cat(sprintf("   GATE P1 (mse@100 <= 2x mse@0, all reps): %s (ratios: %s)\n",
            if (gate_p1) "PASS" else "FAIL",
            paste(sprintf("%.1f", rat), collapse = " ")))

# ===================================================================
#  P2a: k-curves under deviance (both branches need this)
# ===================================================================
cat("\n-- P2a: deviance k-curves (oracle start, E2 seeds) --\n")
t0 <- proc.time()
jobs <- list()
for (rg in names(SIGNAL_LEVELS))
  for (r in seq_len(N_REP_E2))
    jobs[[length(jobs) + 1L]] <- list(regime = rg, rep_id = r)
p2a <- dv_run_jobs(jobs, wrap_try(run_p2a_replicate), NCORES)
t_p2a <- (proc.time() - t0)[3]
saveRDS(list(results = p2a, time_s = t_p2a),
        file.path(RES_DIR, "p2a_results.rds"))
report_ok(p2a, t_p2a)

# ===================================================================
#  P2b: Block 1 grid (only if gate P1 passed)
# ===================================================================
if (gate_p1) {
  cat("\n-- P2b: Block 1 grid, deviance-refined + sandwich --\n")
  t0 <- proc.time()
  jobs <- list()
  for (mi in seq_along(M_GRID))
    for (r in seq_len(N_REP_B1))
      jobs[[length(jobs) + 1L]] <- list(M = M_GRID[mi], M_idx = mi,
                                        rep_id = r)
  p2b <- dv_run_jobs(jobs, wrap_try(run_p2b_replicate), NCORES)
  t_p2b <- (proc.time() - t0)[3]
  saveRDS(list(results = p2b, time_s = t_p2b),
          file.path(RES_DIR, "p2b_results.rds"))
  report_ok(p2b, t_p2b)
} else {
  cat("\n   Gate P1 FAILED: skipping P2b per protocol.\n")
}

# ===================================================================
#  P3: feasible chain under deviance
# ===================================================================
cat("\n-- P3: anchor polish + feasible deviance chain (M = 1000) --\n")
t0 <- proc.time()
jobs <- list()
for (rg in names(SIGNAL_LEVELS))
  for (r in seq_len(N_REP_E2))
    jobs[[length(jobs) + 1L]] <- list(regime = rg, rep_id = r)
p3 <- dv_run_jobs(jobs, wrap_try(run_p3_replicate), NCORES)
t_p3 <- (proc.time() - t0)[3]
saveRDS(list(results = p3, time_s = t_p3),
        file.path(RES_DIR, "p3_results.rds"))
report_ok(p3, t_p3)

if (!QUICK) {
  cat("\n-- P3+: crossover cell (M = 5000, strong, 10 reps) --\n")
  t0 <- proc.time()
  jobs <- lapply(seq_len(10L), function(r)
    list(regime = "strong", rep_id = r, M = 5000L))
  p3x <- dv_run_jobs(jobs, wrap_try(run_p3_replicate), NCORES)
  t_p3x <- (proc.time() - t0)[3]
  saveRDS(list(results = p3x, time_s = t_p3x),
          file.path(RES_DIR, "p3_m5000_results.rds"))
  report_ok(p3x, t_p3x)
}

cat("\n=== Deviance probe complete ===\n")
cat("Build tables with: Rscript replication/deviance_probe/03_report.R\n")
