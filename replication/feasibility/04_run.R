#!/usr/bin/env Rscript
# ===================================================================
#  Feasibility round — runner for F1 (identified estimator),
#  F2 (split-document jackknife), F3 (k rule), F4 (dress rehearsal)
# ===================================================================
#
#  Usage (from the package root):
#    Rscript replication/feasibility/04_run.R [--quick] [--ncores N]
#
#  Order of execution: unit tests (hard gates) -> F3 (validates the
#  adaptive k rule used everywhere downstream) -> F1c -> F1d -> F2 ->
#  F4(ii).  F4(i) is assembled from the F1c results (each replicate
#  already runs the jackknife on V2 and V3); the V2-vs-V3 choice for
#  F4(ii) is made from the F1c summary and logged.
#
#  Seeds: F3/F1c/F4(i) reuse basin_check E2 (90000 + regime*1000 + rep;
#  weak = 1, strong = 2).  F1d: 30000 + alpha_index*1000 + rep.
#  F2 M-grid and F4(ii): audit A3 scheme 60000 + M_index*1000 + rep
#  (M = 500,1000,2000,4000 -> index 1..4).  F2 L-grid: basin_check E5
#  scheme 70000 + dl_index*1000 + rep (weak regime).  Thinning seed =
#  replicate seed + 500.  Unit tests: 55001 (anchors).
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
RES_DIR <- file.path(ROOT, "replication", "feasibility", "results")
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)

fs_setup <- function(root) {
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
  invisible(NULL)
}
fs_setup(ROOT)

SIGNAL_LEVELS <- c(weak = 0.15, strong = 0.50)
K_TOPICS <- 5L; P_COV <- 3L; N_VOCAB <- 500L
Bz0_TRUE <- matrix(c(
   0.40, -0.20,  0.10,  0.30,
  -0.15,  0.35, -0.25,  0.05,
   0.20,  0.10,  0.40, -0.30
), nrow = P_COV, ncol = K_TOPICS - 1L, byrow = TRUE)

N_REP_E2  <- if (QUICK) 2L else 20L
N_REP_B1  <- if (QUICK) 3L else 50L
N_REP_10  <- if (QUICK) 2L else 10L
M_GRID_F2 <- if (QUICK) c(500L, 1000L) else c(500L, 1000L, 2000L, 4000L)
M_GRID_F4 <- if (QUICK) c(500L, 1000L) else c(500L, 1000L, 2000L)
CAP_RULE  <- 50L

# ===================================================================
#  Replicate drivers
# ===================================================================

#' F3: oracle-start unconstrained refinement, 100 sweeps, full B path
run_f3_replicate <- function(regime, rep_id) {
  seed <- 90000L + match(regime, names(SIGNAL_LEVELS)) * 1000L + rep_id
  dat <- sim_dgp(M = 1000L, N = N_VOCAB, K = K_TOPICS, P = P_COV,
                 b_max = SIGNAL_LEVELS[[regime]], sigma_eps = 0.3,
                 alpha_beta = 0.1, doc_length = 200L, seed = seed)
  Wf <- dat$W / rowSums(dat$W)
  fit <- sgscatm(dat$W, dat$C, K = K_TOPICS, lambda = 1, rotate = TRUE)
  gl  <- bc_gl_align(fit$Z, dat$Z_true)
  rf <- fs_refine_rule(gl$Z, bc_phi_step(gl$Z, Wf, dat$V), Wf, dat$C,
                       dat$V, constrained = FALSE, max_sweeps = 100L,
                       apply_rule = FALSE, track_B_path = TRUE,
                       Bz0 = dat$Bz0)
  # criterion-pathology probe (reps 1-5): the simplex-constrained
  #  descent STARTED AT THE TRUTH — F decreases while mse_Bz degrades,
  #  demonstrating that the exact frequency-LS criterion does not
  #  identify B at finite L (constrained or not)
  path_truth <- NULL
  if (rep_id <= 5L) {
    rt <- fs_refine_rule(dat$Z_true, dat$Beta, Wf, dat$C, dat$V,
                         constrained = TRUE, max_sweeps = 100L,
                         apply_rule = FALSE, Bz0 = dat$Bz0)
    path_truth <- list(F_path = rt$F_path, mse_path = rt$mse_path,
                       monotone_ok = rt$monotone_ok,
                       F_truth = bc_objective(dat$Z_true, dat$Beta, Wf,
                                              dat$V))
  }
  ev_all <- fit$eigenvalues_all
  list(regime = regime, rep = rep_id, seed = seed,
       mse_path = rf$mse_path, F_path = rf$F_path,
       rule_stop = fs_rule_from_path(rf$B_path, cap = CAP_RULE),
       monotone_ok = rf$monotone_ok, path_truth = path_truth,
       relgap = (fit$eigenvalues[K_TOPICS - 1L] - ev_all[K_TOPICS]) /
         fit$eigenvalues[1L],
       rho_glZ = sqrt(sum((gl$Z - dat$Z_true)^2)),
       sat_true = mean(apply(dat$Theta_true, 1L, max) > 0.99),
       nBz0 = sqrt(sum(dat$Bz0^2)))
}

#' F1c: feasible variants + oracle references + jackknife on V2/V3
run_f1c_replicate <- function(regime, rep_id) {
  seed <- 90000L + match(regime, names(SIGNAL_LEVELS)) * 1000L + rep_id
  dat <- sim_dgp(M = 1000L, N = N_VOCAB, K = K_TOPICS, P = P_COV,
                 b_max = SIGNAL_LEVELS[[regime]], sigma_eps = 0.3,
                 alpha_beta = 0.1, doc_length = 200L, seed = seed)
  Wf <- dat$W / rowSums(dat$W)
  V  <- dat$V
  nB0 <- sqrt(sum(dat$Bz0^2))
  out <- list(regime = regime, rep = rep_id, seed = seed)
  ev <- function(B) list(paper = procrustes_align(B, dat$Bz0)$mse,
                         perm = fs_perm_mse(B, dat$Bz0, V),
                         norm_ratio = sqrt(sum(B^2)) / nB0)

  # anchors (truth-free)
  t0 <- proc.time()
  ap <- fs_anchor_pipeline(dat$W, K_TOPICS)
  t_anchor <- (proc.time() - t0)[3]
  out$anchor_tv <- fs_perm_tv(ap$Phi, dat$Beta)      # evaluation only
  out$excl_true <- fs_exclusivity(dat$Beta)

  # V1: anchored Phi fixed, per-document GN from z = 0
  t0 <- proc.time()
  zi <- fs_z_init_gn(ap$Phi, Wf, V, n_gn = 10L)
  B1 <- bc_b_step(zi$Z, dat$C)
  out$V1 <- c(ev(B1), list(time_s = t_anchor + (proc.time() - t0)[3],
                           n_fail = zi$n_fail))

  # V2: joint constrained refinement from (Z_V1, Phi_anchor)
  t0 <- proc.time()
  r2 <- fs_refine_rule(zi$Z, ap$Phi, Wf, dat$C, V, constrained = TRUE,
                       max_sweeps = CAP_RULE)
  t_v2 <- (proc.time() - t0)[3]
  jk2 <- fs_jackknife_B(dat$W, dat$C, V, r2$Z, r2$Phi, seed + 500L)
  out$V2 <- c(ev(r2$B), list(sweeps = r2$sweeps, rule_stop = r2$rule_stop,
                             monotone_ok = r2$monotone_ok,
                             time_s = out$V1$time_s + t_v2))
  out$V2_jk <- c(ev(jk2$B_jk), list(dAB = jk2$B_A - jk2$B_B))

  # V3: pilot-mapped init in the anchored frame
  t0 <- proc.time()
  fit <- sgscatm(dat$W, dat$C, K = K_TOPICS, lambda = 1, rotate = TRUE)
  What <- matrix(fit$w_bar, nrow(Wf), N_VOCAB, byrow = TRUE) +
    fit$Z %*% fit$Psi / K_TOPICS                     # pilot reconstruction
  th3 <- fs_simplex_ls(What, ap$Phi)
  th3 <- pmax(th3, 1e-6); th3 <- th3 / rowSums(th3)
  Z3 <- log(th3) %*% V
  t_init3 <- (proc.time() - t0)[3]
  t0 <- proc.time()
  r3 <- fs_refine_rule(Z3, ap$Phi, Wf, dat$C, V, constrained = TRUE,
                       max_sweeps = CAP_RULE)
  t_v3 <- (proc.time() - t0)[3]
  jk3 <- fs_jackknife_B(dat$W, dat$C, V, r3$Z, r3$Phi, seed + 500L)
  out$V3 <- c(ev(r3$B), list(sweeps = r3$sweeps, rule_stop = r3$rule_stop,
                             monotone_ok = r3$monotone_ok,
                             time_s = t_anchor + t_init3 + t_v3))
  out$V3_jk <- c(ev(jk3$B_jk), list(dAB = jk3$B_A - jk3$B_B))

  # V4 (added variant): anchor-oriented pilot — the anchored V1 scores
  # supply a feasible GL orientation for the pilot's near-perfect score
  # subspace (Ahat = OLS of Z_V1 on fit$Z; per-document noise in Z_V1
  # averages out of the 4x4 fit), then UNCONSTRAINED refinement with
  # the k rule.  Motivated by the criterion pathology: constrained
  # long-run descent degrades B, so the feasible candidate is k-step.
  t0 <- proc.time()
  a0 <- colMeans(zi$Z)
  A_hat <- solve(crossprod(fit$Z) + 1e-8 * diag(K_TOPICS - 1L),
                 crossprod(fit$Z, sweep(zi$Z, 2L, a0)))
  Z4 <- sweep(fit$Z %*% A_hat, 2L, a0, "+")
  r4 <- fs_refine_rule(Z4, bc_phi_step(Z4, Wf, V), Wf, dat$C, V,
                       constrained = FALSE, max_sweeps = CAP_RULE)
  jk4 <- fs_jackknife_B(dat$W, dat$C, V, r4$Z, r4$Phi, seed + 500L)
  out$V4 <- c(ev(r4$B), list(sweeps = r4$sweeps,
                             monotone_ok = r4$monotone_ok,
                             time_s = out$V1$time_s + (proc.time() - t0)[3]))
  out$V4_jk <- c(ev(jk4$B_jk), list(dAB = jk4$B_A - jk4$B_B))

  # sandwich + coverage context for F4(i) (feasible variants, jk-centred)
  for (nm in c("V2", "V3", "V4")) {
    r <- switch(nm, V2 = r2, V3 = r3, V4 = r4)
    jk <- switch(nm, V2 = jk2, V3 = jk3, V4 = jk4)
    Sig <- fs_sandwich(r$Z, dat$C, jk$B_full)
    cv <- fs_coverage_entry(jk$B_jk, dat$Bz0, Sig, P_COV, K_TOPICS - 1L)
    out[[paste0(nm, "_cov")]] <- mean(cv$covers)
  }

  # oracle references (truth used for the start alignment)
  gl <- bc_gl_align(fit$Z, dat$Z_true)
  Phi_gl <- bc_phi_step(gl$Z, Wf, V)
  t0 <- proc.time()
  ro5 <- bc_refine(gl$Z, Phi_gl, Wf, dat$C, V, lambda = 0,
                   max_sweeps = 5L)
  out$oracle_k5 <- c(ev(bc_b_step(ro5$Z, dat$C)),
                     list(time_s = (proc.time() - t0)[3]))
  t0 <- proc.time()
  ror <- fs_refine_rule(gl$Z, Phi_gl, Wf, dat$C, V, constrained = FALSE,
                        max_sweeps = CAP_RULE)
  out$oracle_rule <- c(ev(ror$B), list(sweeps = ror$sweeps,
                                       time_s = (proc.time() - t0)[3]))
  out
}

#' F1d: identification boundary in alpha_beta (strong regime)
run_f1d_replicate <- function(alpha, alpha_idx, rep_id) {
  seed <- 30000L + alpha_idx * 1000L + rep_id
  dat <- sim_dgp(M = 1000L, N = N_VOCAB, K = K_TOPICS, P = P_COV,
                 b_max = 0.50, sigma_eps = 0.3, alpha_beta = alpha,
                 doc_length = 200L, seed = seed)
  Wf <- dat$W / rowSums(dat$W)
  ap <- fs_anchor_pipeline(dat$W, K_TOPICS)
  zi <- fs_z_init_gn(ap$Phi, Wf, dat$V, n_gn = 10L)
  r2 <- fs_refine_rule(zi$Z, ap$Phi, Wf, dat$C, dat$V,
                       constrained = TRUE, max_sweeps = CAP_RULE)
  list(alpha = alpha, rep = rep_id, seed = seed,
       excl_true = fs_exclusivity(dat$Beta),
       anchor_tv = fs_perm_tv(ap$Phi, dat$Beta),
       mse_paper = procrustes_align(r2$B, dat$Bz0)$mse,
       mse_perm = fs_perm_mse(r2$B, dat$Bz0, dat$V),
       sweeps = r2$sweeps, monotone_ok = r2$monotone_ok)
}

#' F2: oracle-start estimator + jackknife on the Block 1 grid / L-grid
run_f2_replicate <- function(M, M_idx, rep_id, dl = 200L,
                             use_block1 = TRUE, seed_base = 60000L) {
  seed <- seed_base + M_idx * 1000L + rep_id
  dat <- if (use_block1) {
    sim_dgp(M = M, N = N_VOCAB, K = K_TOPICS, P = P_COV,
            Bz0 = Bz0_TRUE, sigma_eps = 0.3, alpha_beta = 0.1,
            doc_length = dl, seed = seed)
  } else {                                            # E5 pairing (weak)
    sim_dgp(M = M, N = N_VOCAB, K = K_TOPICS, P = P_COV,
            b_max = SIGNAL_LEVELS[["weak"]], sigma_eps = 0.3,
            alpha_beta = 0.1, doc_length = dl, seed = seed)
  }
  Wf <- dat$W / rowSums(dat$W)
  V  <- dat$V
  t0 <- proc.time()
  fit <- sgscatm(dat$W, dat$C, K = K_TOPICS, lambda = 1, rotate = TRUE)
  gl  <- bc_gl_align(fit$Z, dat$Z_true)
  rf <- fs_refine_rule(gl$Z, bc_phi_step(gl$Z, Wf, V), Wf, dat$C, V,
                       constrained = FALSE, max_sweeps = CAP_RULE)
  jk <- fs_jackknife_B(dat$W, dat$C, V, rf$Z, rf$Phi, seed + 500L)
  t_all <- (proc.time() - t0)[3]

  Sig <- fs_sandwich(rf$Z, dat$C, jk$B_full)
  Kp <- K_TOPICS - 1L
  var_add <- (jk$B_A - jk$B_B)^2 / 4
  cv_full <- fs_coverage_entry(jk$B_full, dat$Bz0, Sig, P_COV, Kp)
  cv_jk   <- fs_coverage_entry(jk$B_jk,   dat$Bz0, Sig, P_COV, Kp)
  cv_jki  <- fs_coverage_entry(jk$B_jk,   dat$Bz0, Sig, P_COV, Kp,
                               var_add = var_add)
  list(M = M, dl = dl, rep = rep_id, seed = seed,
       mse_full = cv_full$mse, mse_jk = cv_jk$mse,
       B_tilde_full = fs_coverage_entry(jk$B_full, dat$Bz0, Sig,
                                        P_COV, Kp)$mse * 0 +
         procrustes_align(jk$B_full, dat$Bz0)$Bz_aligned,
       B_tilde_jk = procrustes_align(jk$B_jk, dat$Bz0)$Bz_aligned,
       cov_full = mean(cv_full$covers), cov_jk = mean(cv_jk$covers),
       cov_jk_infl = mean(cv_jki$covers),
       std_err_jk = cv_jk$std_err,
       rn_full = fs_coverage_rownorm(jk$B_full, dat$Bz0, Sig, P_COV, Kp),
       rn_jk   = fs_coverage_rownorm(jk$B_jk, dat$Bz0, Sig, P_COV, Kp),
       rn_jk_infl = fs_coverage_rownorm(jk$B_jk, dat$Bz0, Sig, P_COV, Kp,
                                        var_add = var_add),
       sweeps = rf$sweeps, monotone_ok = rf$monotone_ok,
       Bz0 = dat$Bz0, time_s = t_all)
}

#' F4(ii): full feasible chain on the Block 1 grid
run_f4_replicate <- function(M, M_idx, rep_id, init_variant = "V3") {
  seed <- 60000L + M_idx * 1000L + rep_id             # audit A3 scheme
  dat <- sim_dgp(M = M, N = N_VOCAB, K = K_TOPICS, P = P_COV,
                 Bz0 = Bz0_TRUE, sigma_eps = 0.3, alpha_beta = 0.1,
                 doc_length = 200L, seed = seed)
  Wf <- dat$W / rowSums(dat$W)
  V  <- dat$V
  t0 <- proc.time()
  ap <- fs_anchor_pipeline(dat$W, K_TOPICS)
  if (init_variant == "V3") {
    fit <- sgscatm(dat$W, dat$C, K = K_TOPICS, lambda = 1, rotate = TRUE)
    What <- matrix(fit$w_bar, nrow(Wf), N_VOCAB, byrow = TRUE) +
      fit$Z %*% fit$Psi / K_TOPICS
    th <- fs_simplex_ls(What, ap$Phi)
    th <- pmax(th, 1e-6); th <- th / rowSums(th)
    Z0 <- log(th) %*% V
  } else {
    Z0 <- fs_z_init_gn(ap$Phi, Wf, V, n_gn = 10L)$Z
  }
  if (init_variant == "V4") {
    # anchor-oriented pilot + UNCONSTRAINED k-rule refinement (the
    # constrained criterion demonstrably degrades B; see F3 pathology)
    fit <- sgscatm(dat$W, dat$C, K = K_TOPICS, lambda = 1, rotate = TRUE)
    a0 <- colMeans(Z0)
    A_hat <- solve(crossprod(fit$Z) + 1e-8 * diag(K_TOPICS - 1L),
                   crossprod(fit$Z, sweep(Z0, 2L, a0)))
    Z0 <- sweep(fit$Z %*% A_hat, 2L, a0, "+")
    rf <- fs_refine_rule(Z0, bc_phi_step(Z0, Wf, V), Wf, dat$C, V,
                         constrained = FALSE, max_sweeps = CAP_RULE)
  } else {
    rf <- fs_refine_rule(Z0, ap$Phi, Wf, dat$C, V, constrained = TRUE,
                         max_sweeps = CAP_RULE)
  }
  jk <- fs_jackknife_B(dat$W, dat$C, V, rf$Z, rf$Phi, seed + 500L)
  t_all <- (proc.time() - t0)[3]

  Sig <- fs_sandwich(rf$Z, dat$C, jk$B_full)
  Kp <- K_TOPICS - 1L
  var_add <- (jk$B_A - jk$B_B)^2 / 4
  cv_jk  <- fs_coverage_entry(jk$B_jk, dat$Bz0, Sig, P_COV, Kp)
  cv_jki <- fs_coverage_entry(jk$B_jk, dat$Bz0, Sig, P_COV, Kp,
                              var_add = var_add)
  list(M = M, rep = rep_id, seed = seed,
       mse_full = procrustes_align(jk$B_full, dat$Bz0)$mse,
       mse_jk = cv_jk$mse,
       mse_perm_jk = fs_perm_mse(jk$B_jk, dat$Bz0, V),
       cov_jk = mean(cv_jk$covers), cov_jk_infl = mean(cv_jki$covers),
       rn_jk = fs_coverage_rownorm(jk$B_jk, dat$Bz0, Sig, P_COV, Kp),
       rn_jk_infl = fs_coverage_rownorm(jk$B_jk, dat$Bz0, Sig, P_COV, Kp,
                                        var_add = var_add),
       sweeps = rf$sweeps, monotone_ok = rf$monotone_ok,
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
fs_run_jobs <- function(jobs, fun, ncores) {
  if (ncores <= 1L || length(jobs) <= 1L)
    return(lapply(jobs, function(j) do.call(fun, j)))
  cl <- parallel::makePSOCKcluster(min(ncores, length(jobs)))
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(cl, fs_setup, root = ROOT)
  parallel::clusterExport(cl, c("SIGNAL_LEVELS", "K_TOPICS", "P_COV",
                                "N_VOCAB", "Bz0_TRUE", "CAP_RULE",
                                "run_f3_replicate", "run_f1c_replicate",
                                "run_f1d_replicate", "run_f2_replicate",
                                "run_f4_replicate"),
                          envir = globalenv())
  wf <- function(j) do.call(FUN, j)
  environment(wf) <- list2env(list(FUN = fun), parent = globalenv())
  parallel::parLapplyLB(cl, jobs, wf)
}
report_ok <- function(res, label, t_el) {
  n_ok <- sum(vapply(res, function(x) is.null(x$error), logical(1)))
  cat(sprintf("   [%.1f s] %d/%d replicates OK\n", t_el, n_ok, length(res)))
  if (n_ok < length(res)) {
    msgs <- unique(unlist(lapply(res, function(x) x$error)))
    cat("   ERRORS:", paste(head(msgs, 3), collapse = " | "), "\n")
  }
}

cat("=== Feasibility round ===\n")
cat("Mode:", if (QUICK) "QUICK" else "FULL", "| cores:", NCORES, "\n")

# ===================================================================
#  Unit tests (hard gates)
# ===================================================================
cat("\n-- F1a unit test: anchor recovery at alpha_beta = 0.05 --\n")
t0 <- proc.time()
at <- fs_test_anchors()
cat(sprintf(paste0("   mean row TV = %.4f (gate <= 0.15) | perm-metric",
                   " self-test = %.2e | true exclusivity = %.3f\n"),
            at$tv, at$perm_zero, at$exclusivity_true))
saveRDS(at, file.path(RES_DIR, "f1a_unit_test.rds"))
if (!at$pass) stop("F1a anchor unit test FAILED (TV > 0.15).")
if (at$perm_zero > 1e-20) stop("Permutation-representation self-test FAILED.")
cat(sprintf("   [%.1f s] PASSED\n", (proc.time() - t0)[3]))

# ===================================================================
#  F3: k rule (first — the rule feeds everything downstream)
# ===================================================================
cat("\n-- F3: k curves, oracle-start, E2 seeds, 100 sweeps --\n")
t0 <- proc.time()
jobs <- list()
for (rg in names(SIGNAL_LEVELS))
  for (r in seq_len(N_REP_E2))
    jobs[[length(jobs) + 1L]] <- list(regime = rg, rep_id = r)
f3 <- fs_run_jobs(jobs, wrap_try(run_f3_replicate), NCORES)
t_f3 <- (proc.time() - t0)[3]
saveRDS(list(results = f3, time_s = t_f3),
        file.path(RES_DIR, "f3_results.rds"))
report_ok(f3, "F3", t_f3)

# ===================================================================
#  F1c: feasible variants vs oracle reference (E2 seeds)
# ===================================================================
cat("\n-- F1c: V1/V2/V3 vs oracle, E2 seeds --\n")
t0 <- proc.time()
f1c <- fs_run_jobs(jobs, wrap_try(run_f1c_replicate), NCORES)
t_f1c <- (proc.time() - t0)[3]
saveRDS(list(results = f1c, time_s = t_f1c),
        file.path(RES_DIR, "f1c_results.rds"))
report_ok(f1c, "F1c", t_f1c)

# choose the F4 init variant from the F1c summary (lower mean Procrustes
# MSE across both regimes, jackknifed)
okc <- Filter(function(x) is.null(x$error), f1c)
mjk <- vapply(c("V2_jk", "V3_jk", "V4_jk"), function(nm)
  mean(vapply(okc, function(x) x[[nm]]$paper, numeric(1))), numeric(1))
INIT_F4 <- c("V2", "V3", "V4")[which.min(mjk)]
cat(sprintf("   F4 init variant: %s (jk mse V2 %.5f | V3 %.5f | V4 %.5f)\n",
            INIT_F4, mjk[1], mjk[2], mjk[3]))

# ===================================================================
#  F1d: identification boundary (strong regime)
# ===================================================================
if (!QUICK) {
  cat("\n-- F1d: alpha_beta boundary --\n")
  ALPHAS <- c(0.05, 0.1, 0.3, 1.0)
  t0 <- proc.time()
  jobs_d <- list()
  for (ai in seq_along(ALPHAS))
    for (r in seq_len(N_REP_10))
      jobs_d[[length(jobs_d) + 1L]] <- list(alpha = ALPHAS[ai],
                                            alpha_idx = ai, rep_id = r)
  f1d <- fs_run_jobs(jobs_d, wrap_try(run_f1d_replicate), NCORES)
  t_f1d <- (proc.time() - t0)[3]
  saveRDS(list(results = f1d, time_s = t_f1d),
          file.path(RES_DIR, "f1d_results.rds"))
  report_ok(f1d, "F1d", t_f1d)
}

# ===================================================================
#  F2: jackknife on the Block 1 grid + L-grid
# ===================================================================
cat("\n-- F2: Block 1 grid, oracle-start + jackknife --\n")
t0 <- proc.time()
jobs_2 <- list()
for (mi in seq_along(M_GRID_F2))
  for (r in seq_len(N_REP_B1))
    jobs_2[[length(jobs_2) + 1L]] <- list(M = M_GRID_F2[mi], M_idx = mi,
                                          rep_id = r)
f2 <- fs_run_jobs(jobs_2, wrap_try(run_f2_replicate), NCORES)
t_f2 <- (proc.time() - t0)[3]
saveRDS(list(results = f2, time_s = t_f2),
        file.path(RES_DIR, "f2_results.rds"))
report_ok(f2, "F2", t_f2)

cat("\n-- F2: L-grid (E5 pairing, weak regime) --\n")
t0 <- proc.time()
DLS <- c(50L, 200L, 1000L)
jobs_l <- list()
for (di in seq_along(DLS))
  for (r in seq_len(N_REP_10))
    jobs_l[[length(jobs_l) + 1L]] <- list(M = 1000L, M_idx = di,
                                          rep_id = r, dl = DLS[di],
                                          use_block1 = FALSE,
                                          seed_base = 70000L)
f2l <- fs_run_jobs(jobs_l, wrap_try(run_f2_replicate), NCORES)
t_f2l <- (proc.time() - t0)[3]
saveRDS(list(results = f2l, time_s = t_f2l),
        file.path(RES_DIR, "f2_lgrid_results.rds"))
report_ok(f2l, "F2-L", t_f2l)

# ===================================================================
#  F4(ii): full feasible chain on the Block 1 grid
# ===================================================================
cat(sprintf("\n-- F4(ii): feasible chain (init %s) on the Block 1 grid --\n",
            INIT_F4))
t0 <- proc.time()
jobs_4 <- list()
for (mi in seq_along(M_GRID_F4))
  for (r in seq_len(N_REP_B1))
    jobs_4[[length(jobs_4) + 1L]] <- list(M = M_GRID_F4[mi], M_idx = mi,
                                          rep_id = r,
                                          init_variant = INIT_F4)
f4 <- fs_run_jobs(jobs_4, wrap_try(run_f4_replicate), NCORES)
t_f4 <- (proc.time() - t0)[3]
saveRDS(list(results = f4, init_variant = INIT_F4, time_s = t_f4),
        file.path(RES_DIR, "f4_results.rds"))
report_ok(f4, "F4", t_f4)

cat("\n=== Feasibility runs complete ===\n")
cat("Build tables with: Rscript replication/feasibility/05_report.R\n")
