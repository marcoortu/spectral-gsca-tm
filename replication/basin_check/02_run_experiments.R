#!/usr/bin/env Rscript
# ===================================================================
#  Basin-condition verification for the Newton refinement (Block 3)
# ===================================================================
#
#  Experiments (see replication/basin_check/REPORT.md):
#    E0  smoke test        — tiny config, full pipeline end-to-end
#    E1  operational basin — refine from pilot vs from truth, endpoints
#    E2  accuracy table    — mse_Bz: pilot | refined | oracle floor
#    E3  Hessian at truth  — spectrum, gauge subspace, gamma_perp
#    E4  Kantorovich       — L_H, rho_perp, basin ratio r
#    E5  bias floor        — refined-from-truth vs doc_length (optional)
#
#  Usage (from the package root):
#    Rscript replication/basin_check/02_run_experiments.R            # full
#    Rscript replication/basin_check/02_run_experiments.R --quick    # E0 only
#    Rscript replication/basin_check/02_run_experiments.R --no-e5
#
#  Seed scheme: seed = 90000 + regime_index*1000 + rep
#  (regime_index: weak = 1, strong = 2); E0 uses 80000 + ..., E5 uses
#  70000 + dl_index*1000 + rep.  Bz0 is drawn INSIDE sim_dgp() under
#  the replicate seed (entries ~ U(-b_max, b_max)), mirroring Block 3's
#  per-replicate random Bz0.
# ===================================================================

args    <- commandArgs(trailingOnly = TRUE)
QUICK   <- "--quick" %in% args
RUN_E5  <- !("--no-e5" %in% args) && !QUICK
NCORES  <- if ("--ncores" %in% args) {
  as.integer(args[which(args == "--ncores") + 1L])
} else {
  max(1L, min(10L, parallel::detectCores() - 2L))
}

ROOT <- Sys.getenv("SGSCATM_ROOT", unset = normalizePath(".", mustWork = FALSE))
if (!dir.exists(file.path(ROOT, "R")))
  stop("Run from the package root or set SGSCATM_ROOT.")
BC_DIR  <- file.path(ROOT, "replication", "basin_check")
RES_DIR <- file.path(BC_DIR, "results")
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)

# --- source package + simulation infrastructure (no modifications) --
bc_setup <- function(root) {
  source(file.path(root, "R", "sgscatm_fit.R"))
  source(file.path(root, "R", "ilr_contrast.R"))
  source(file.path(root, "R", "utils.R"))
  source(file.path(root, "replication", "simulation", "sim_dgp.R"))
  source(file.path(root, "replication", "simulation", "sim_utils.R"))
  source(file.path(root, "replication", "basin_check", "01_functions.R"))
  invisible(NULL)
}
bc_setup(ROOT)

SIGNAL_LEVELS <- c(weak = 0.15, strong = 0.50)   # Block 3 b_max values

# --- design configs -------------------------------------------------
CFG_FULL <- list(
  M = 1000L, N = 500L, K = 5L, P = 3L,
  sigma_eps = 0.3, alpha_beta = 0.1, doc_length = 200L,
  n_rep_e12 = 20L, n_rep_e34 = 5L, n_rep_e5 = 10L,
  k_eigs = 30L, seed_base = 90000L
)
CFG_SMOKE <- list(
  M = 200L, N = 200L, K = 5L, P = 3L,
  sigma_eps = 0.3, alpha_beta = 0.1, doc_length = 200L,
  n_rep_e12 = 2L, n_rep_e34 = 1L, n_rep_e5 = 2L,
  k_eigs = 10L, seed_base = 80000L
)

cat("=== Basin-condition verification ===\n")
cat("Mode:", if (QUICK) "QUICK (smoke only)" else "FULL", "| cores:", NCORES, "\n")


# ===================================================================
#  Derivative verification gate (hard stop on failure)
# ===================================================================
cat("\n-- Derivative verification (tiny instance, M=30 N=40 K=3) --\n")
t0 <- proc.time()
chk <- bc_verify_derivatives()
cat(sprintf("   grad rel.err: lambda=0 %.2e | lambda=1 %.2e\n",
            chk$grad_relerr_lambda0, chk$grad_relerr_lambda1))
cat(sprintf("   HVP  rel.err: lambda=0 %.2e | lambda=1 %.2e\n",
            chk$hvp_relerr_lambda0, chk$hvp_relerr_lambda1))
cat(sprintf("   gauge tangency (rel): %.2e\n", chk$gauge_tangency))
cat(sprintf("   [%.1f s] PASSED\n", (proc.time() - t0)[3]))
saveRDS(chk, file.path(RES_DIR, "derivative_checks.rds"))


# ===================================================================
#  Per-replicate drivers
# ===================================================================

#' E1 + E2 replicate: pilot, alignments, four refinements, endpoints
run_e12_replicate <- function(regime, rep_id, cfg) {
  regime_idx <- match(regime, names(SIGNAL_LEVELS))
  seed <- cfg$seed_base + regime_idx * 1000L + rep_id
  dat <- sim_dgp(M = cfg$M, N = cfg$N, K = cfg$K, P = cfg$P,
                 b_max = SIGNAL_LEVELS[[regime]],
                 sigma_eps = cfg$sigma_eps, alpha_beta = cfg$alpha_beta,
                 doc_length = cfg$doc_length, seed = seed)
  Wf <- dat$W / rowSums(dat$W)
  V  <- dat$V

  t0  <- proc.time()
  fit <- sgscatm(dat$W, dat$C, K = cfg$K, lambda = 1, rotate = TRUE)
  t_pilot <- (proc.time() - t0)[3]

  # pilot metrics: paper metric exactly as Block 3 computes it
  pilot_paper <- procrustes_align(fit$Bz, dat$Bz0)$mse

  # alignments of the pilot scores into model coordinates
  gl <- bc_gl_align(fit$Z, dat$Z_true)
  op <- bc_op_align(fit$Z, dat$Z_true)
  Phi_gl <- bc_phi_step(gl$Z, Wf, V)
  Phi_op <- bc_phi_step(op$Z, Wf, V)

  pilot <- list(
    paper_mse   = pilot_paper,
    gl_mse      = bc_mse_direct(gl$Z, dat$C, dat$Bz0),
    op_mse      = bc_mse_direct(op$Z, dat$C, dat$Bz0),
    op_mse_proc = bc_mse_paper(op$Z, dat$C, dat$Bz0),
    rho_gl      = bc_rho(gl$Z, Phi_gl, dat$Z_true, dat$Beta),
    rho_op      = bc_rho(op$Z, Phi_op, dat$Z_true, dat$Beta),
    A_gl        = gl$A,
    Bz_norm     = sqrt(sum(fit$Bz^2)),
    Bz0_norm    = sqrt(sum(dat$Bz0^2)),
    time_s      = t_pilot
  )

  refits <- list()
  for (lam in c(0, 1)) {
    # lambda = 1 uses the profiled two-block solver and a 600-sweep cap:
    # the three-block variant zig-zags along the Z ~ C B valley (rate
    # ~0.989/sweep) and even the profiled tail decays at ~0.98/sweep,
    # so the brief's 100-sweep cap would leave endpoints meaningless.
    msw <- if (lam > 0) 600L else 100L
    t1 <- proc.time()
    rp <- bc_refine(gl$Z, Phi_gl, Wf, dat$C, V, lambda = lam,
                    Bz0 = dat$Bz0, max_sweeps = msw)
    t_rp <- (proc.time() - t1)[3]
    t1 <- proc.time()
    rt <- bc_refine(dat$Z_true, dat$Beta, Wf, dat$C, V, lambda = lam,
                    Bz0 = dat$Bz0, max_sweeps = msw)
    t_rt <- (proc.time() - t1)[3]

    dZ <- sqrt(sum((rp$Z - rt$Z)^2)) / sqrt(sum(rt$Z^2))
    dF <- abs(rp$F - rt$F) / (1 + abs(rt$F))

    # gauge-invariant endpoint diagnostics: at lambda = 0 the minimiser
    # is a 20-dim gauge orbit, so endpoints can coincide as models while
    # differing in Z.  Compare the fitted matrices Theta(Z)Phi (exactly
    # gauge-invariant) and the difference projected off the gauge
    # tangent at the truth-start endpoint.
    Fit_p <- bc_theta(rp$Z, V) %*% rp$Phi
    Fit_t <- bc_theta(rt$Z, V) %*% rt$Phi
    dFit  <- sqrt(sum((Fit_p - Fit_t)^2)) / sqrt(sum(Fit_t^2))
    QTb   <- bc_gauge_basis(rt$Z, rt$Phi, V)
    dv    <- bc_pack(rp$Z, rp$Phi) - bc_pack(rt$Z, rt$Phi)
    dvp   <- dv - QTb %*% crossprod(QTb, dv)
    dEta_perp <- sqrt(sum(dvp^2)) /
      sqrt(sum(bc_pack(rt$Z, rt$Phi)^2))

    refits[[sprintf("lambda%g", lam)]] <- list(
      lambda = lam,
      # E1: endpoint coincidence (strict pre-registered criterion)
      dZ_rel = dZ, dF_rel = dF,
      same_basin = (dZ < 1e-3) && (dF < 1e-8),
      # E1': gauge-aware coincidence
      dFit_rel = dFit, dEta_perp_rel = dEta_perp,
      same_basin_gauge = (dFit < 1e-3) && (dF < 1e-8),
      sweeps_pilot = rp$sweeps, sweeps_truth = rt$sweeps,
      converged_pilot = rp$converged, converged_truth = rt$converged,
      monotone_ok = !rp$monotone_violation && !rt$monotone_violation,
      # E2: accuracy under both metrics
      refined_paper = bc_mse_paper(rp$Z, dat$C, dat$Bz0),
      refined_gl    = bc_mse_direct(rp$Z, dat$C, dat$Bz0),
      truthref_paper = bc_mse_paper(rt$Z, dat$C, dat$Bz0),
      truthref_gl    = bc_mse_direct(rt$Z, dat$C, dat$Bz0),
      rho_refined = bc_rho(rp$Z, rp$Phi, dat$Z_true, dat$Beta),
      F_pilot_end = rp$F, F_truth_end = rt$F,
      trace_pilot = rp$trace, trace_truth = rt$trace,
      time_pilot_s = t_rp, time_truth_s = t_rt
    )
  }

  list(regime = regime, rep = rep_id, seed = seed,
       pilot = pilot, refits = refits)
}

#' Spectral analysis of the (profiled) Hessian at a point eta
#'
#' Returns lambda_max (power iteration), the k numerically smallest raw
#' eigenvalues, gauge diagnostics (||H q||/lambda_max per gauge column,
#' principal cosines between the gauge span and the smallest raw
#' eigenvectors) and, at lambda = 0, gamma_perp from deflated Lanczos.
bc_analyze_hessian <- function(eta, dims, Wf, V, lam, C, QT, d, k, seed) {
  hvp <- bc_hvp_factory(eta, dims, Wf, V, lam, C)
  pw  <- bc_power_iter(hvp, d, iters = 50L, seed = seed)
  shift <- 1.01 * pw$norm
  raw <- bc_smallest_eigs(hvp, d, k = k, shift = shift)
  res <- list(lambda_max = pw$value, raw_smallest = raw$values,
              raw_nconv = raw$nconv, raw_niter = raw$niter)
  # direct gauge-block diagnostics: ARPACK cannot resolve the ~20-fold
  # degenerate near-null cluster from a single Krylov sequence, so the
  # trustworthy checks are ||H q|| per gauge column and the spectrum of
  # the 20 x 20 compression Q_T' H Q_T.
  HQ <- vapply(seq_len(ncol(QT)), function(j) hvp(QT[, j]), numeric(d))
  res$HQt_relnorms <- sqrt(colSums(HQ^2)) / pw$norm
  QtHQ <- crossprod(QT, HQ)
  res$QtHQ_eigs <- eigen((QtHQ + t(QtHQ)) / 2, symmetric = TRUE,
                         only.values = TRUE)$values
  kk <- min(ncol(QT), ncol(raw$vectors))
  res$principal_cosines <-
    bc_principal_cosines(QT, raw$vectors[, seq_len(kk), drop = FALSE])
  if (lam == 0) {
    defl <- bc_smallest_eigs(hvp, d, k = k, shift = shift, QT = QT)
    res$deflated_smallest <- defl$values
    res$gamma_perp <- defl$values[1L]
  } else {
    res$gamma_raw <- raw$values[1L]
  }
  res
}

#' E3 + E4 replicate: Hessian spectrum at the truth eta0 AND at the
#' exact M-estimator eta_star (refined-from-truth endpoint), plus the
#' Kantorovich diagnostic.  eta_star analysis is an addition to the
#' brief: at finite samples eta0 is not a stationary point of F, so
#' H(eta0) can be slightly indefinite; the Newton-Kantorovich gamma
#' lives at eta_star.
run_e34_replicate <- function(regime, rep_id, cfg) {
  regime_idx <- match(regime, names(SIGNAL_LEVELS))
  seed <- cfg$seed_base + regime_idx * 1000L + rep_id
  dat <- sim_dgp(M = cfg$M, N = cfg$N, K = cfg$K, P = cfg$P,
                 b_max = SIGNAL_LEVELS[[regime]],
                 sigma_eps = cfg$sigma_eps, alpha_beta = cfg$alpha_beta,
                 doc_length = cfg$doc_length, seed = seed)
  Wf <- dat$W / rowSums(dat$W)
  V  <- dat$V
  dims <- list(M = cfg$M, Km1 = cfg$K - 1L, K = cfg$K, N = cfg$N)
  d    <- cfg$M * (cfg$K - 1L) + cfg$K * cfg$N

  eta0 <- bc_pack(dat$Z_true, dat$Beta)
  QT   <- bc_gauge_basis(dat$Z_true, dat$Beta, V)

  # pilot in model coordinates (GL alignment), for the segment in E4
  fit <- sgscatm(dat$W, dat$C, K = cfg$K, lambda = 1, rotate = TRUE)
  gl  <- bc_gl_align(fit$Z, dat$Z_true)
  Phi_gl <- bc_phi_step(gl$Z, Wf, V)
  eta_pil <- bc_pack(gl$Z, Phi_gl)

  diffv    <- eta_pil - eta0
  rho_full <- sqrt(sum(diffv^2))
  dperp    <- diffv - QT %*% crossprod(QT, diffv)
  rho_perp <- sqrt(sum(dperp^2))

  out <- list(regime = regime, rep = rep_id, seed = seed, d = d,
              rho_full = rho_full, rho_perp = rho_perp)

  for (lam in c(0, 1)) {
    key <- sprintf("lambda%g", lam)
    t0  <- proc.time()
    res <- list(lambda = lam)

    # --- E3a: Hessian at the truth eta0 (pre-registered) -------------
    res$at_truth <- bc_analyze_hessian(eta0, dims, Wf, V, lam, dat$C,
                                       QT = QT, d = d, k = cfg$k_eigs,
                                       seed = seed + round(lam))

    # --- E3b: Hessian at the exact M-estimator eta_star ---------------
    rt <- bc_refine(dat$Z_true, dat$Beta, Wf, dat$C, V, lambda = lam,
                    max_sweeps = if (lam > 0) 600L else 200L)
    eta_star <- bc_pack(rt$Z, rt$Phi)
    QTs <- bc_gauge_basis(rt$Z, rt$Phi, V)
    res$at_star <- bc_analyze_hessian(eta_star, dims, Wf, V, lam, dat$C,
                                      QT = QTs, d = d, k = cfg$k_eigs,
                                      seed = seed + round(lam) + 100L)
    res$star_F <- rt$F
    res$star_sweeps <- rt$sweeps
    res$star_converged <- rt$converged
    res$rho_truth_to_star <- sqrt(sum((eta_star - eta0)^2))

    # pilot displacement measured from eta_star (basin of the estimator)
    dstar <- eta_pil - eta_star
    res$rho_star_full <- sqrt(sum(dstar^2))
    dsp <- dstar - QTs %*% crossprod(QTs, dstar)
    res$rho_star_perp <- sqrt(sum(dsp^2))

    # --- E4: local Hessian Lipschitz constant along eta0 -> eta_pil ---
    pts <- seq(0, 1, length.out = 5L)
    seg_len <- rho_full / 4
    ratios <- numeric(4L)
    for (j in seq_len(4L)) {
      eta_s <- eta0 + pts[j]     * diffv
      eta_t <- eta0 + pts[j + 1] * diffv
      hs <- bc_hvp_factory(eta_s, dims, Wf, V, lam, dat$C)
      ht <- bc_hvp_factory(eta_t, dims, Wf, V, lam, dat$C)
      Dop <- function(v) ht(v) - hs(v)
      ratios[j] <- bc_power_iter(Dop, d, iters = 30L,
                                 seed = seed + j)$norm / seg_len
    }
    res$L_H_ratios <- ratios
    res$L_H <- max(ratios)

    # basin ratios (Kantorovich is sufficient, not necessary — reported
    # as is): anchored at eta0 (pre-registered) and at eta_star
    gam0 <- if (lam == 0) res$at_truth$gamma_perp else res$at_truth$gamma_raw
    gams <- if (lam == 0) res$at_star$gamma_perp  else res$at_star$gamma_raw
    res$gamma_truth <- gam0
    res$gamma_star  <- gams
    res$r_truth <- 2 * res$L_H * rho_perp / gam0
    res$r_star  <- 2 * res$L_H * res$rho_star_perp / gams
    res$time_s <- (proc.time() - t0)[3]
    out[[key]] <- res
  }
  out
}

#' E5 replicate: refined-from-truth accuracy vs document length
run_e5_replicate <- function(dl, dl_idx, rep_id, cfg) {
  seed <- 70000L + dl_idx * 1000L + rep_id
  dat <- sim_dgp(M = cfg$M, N = cfg$N, K = cfg$K, P = cfg$P,
                 b_max = SIGNAL_LEVELS[["weak"]],
                 sigma_eps = cfg$sigma_eps, alpha_beta = cfg$alpha_beta,
                 doc_length = dl, seed = seed)
  Wf <- dat$W / rowSums(dat$W)
  rt <- bc_refine(dat$Z_true, dat$Beta, Wf, dat$C, dat$V, lambda = 0,
                  Bz0 = dat$Bz0)
  list(doc_length = dl, rep = rep_id, seed = seed,
       mse_paper = bc_mse_paper(rt$Z, dat$C, dat$Bz0),
       mse_gl    = bc_mse_direct(rt$Z, dat$C, dat$Bz0),
       rho_end   = bc_rho(rt$Z, rt$Phi, dat$Z_true, dat$Beta),
       sweeps = rt$sweeps, converged = rt$converged,
       monotone_ok = !rt$monotone_violation)
}


# ===================================================================
#  Parallel execution helper (PSOCK: mclapply is serial on Windows)
# ===================================================================
bc_run_jobs <- function(jobs, fun, cfg, ncores) {
  if (ncores <= 1L || length(jobs) <= 1L) {
    return(lapply(jobs, function(j) do.call(fun, c(j, list(cfg = cfg)))))
  }
  cl <- parallel::makePSOCKcluster(min(ncores, length(jobs)))
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(cl, bc_setup, root = ROOT)
  parallel::clusterExport(cl, c("SIGNAL_LEVELS", "bc_analyze_hessian",
                                "run_e12_replicate", "run_e34_replicate",
                                "run_e5_replicate"),
                          envir = globalenv())
  # worker fn gets a minimal environment (the frame of bc_run_jobs holds
  # the cluster connection object, which must not be serialised)
  wf <- function(j) do.call(FUN, c(j, list(cfg = CFG)))
  environment(wf) <- list2env(list(FUN = fun, CFG = cfg),
                              parent = globalenv())
  parallel::parLapplyLB(cl, jobs, wf)
}

wrap_try <- function(fun) {
  force(fun)   # force the promise: it must serialise to PSOCK workers
  function(...) {
    tryCatch(fun(...), error = function(e)
      list(error = conditionMessage(e), args = list(...)))
  }
}


# ===================================================================
#  E0: smoke test (always runs first; gate for the full runs)
# ===================================================================
cat("\n-- E0: smoke test (M=200, N=200, 2 reps, both regimes) --\n")
t0 <- proc.time()
jobs_e0 <- do.call(rbind, lapply(names(SIGNAL_LEVELS), function(rg)
  data.frame(regime = rg, rep_id = seq_len(CFG_SMOKE$n_rep_e12),
             stringsAsFactors = FALSE)))
e0_e12 <- lapply(seq_len(nrow(jobs_e0)), function(i)
  wrap_try(run_e12_replicate)(jobs_e0$regime[i], jobs_e0$rep_id[i],
                              cfg = CFG_SMOKE))
e0_e34 <- lapply(names(SIGNAL_LEVELS), function(rg)
  wrap_try(run_e34_replicate)(rg, 1L, cfg = CFG_SMOKE))
t_e0 <- (proc.time() - t0)[3]
saveRDS(list(e12 = e0_e12, e34 = e0_e34, cfg = CFG_SMOKE, time_s = t_e0),
        file.path(RES_DIR, "e0_smoke.rds"))

e0_err <- c(vapply(e0_e12, function(x) !is.null(x$error), logical(1)),
            vapply(e0_e34, function(x) !is.null(x$error), logical(1)))
e0_mono <- vapply(e0_e12, function(x) {
  if (!is.null(x$error)) return(FALSE)
  all(vapply(x$refits, `[[`, logical(1), "monotone_ok"))
}, logical(1))
cat(sprintf("   [%.1f s] errors: %d/%d | monotone: %s\n", t_e0,
            sum(e0_err), length(e0_err),
            if (all(e0_mono)) "all OK" else "VIOLATION"))
if (any(e0_err)) {
  msgs <- unlist(lapply(c(e0_e12, e0_e34), function(x) x$error))
  stop("E0 smoke test FAILED:\n", paste(msgs, collapse = "\n"))
}
if (!all(e0_mono)) stop("E0 smoke test FAILED: monotonicity violation.")
cat("   E0 PASSED\n")

if (QUICK) {
  cat("\nQUICK mode: stopping after E0.\n")
  quit(save = "no", status = 0L)
}


# ===================================================================
#  E1 + E2: full-size refinement study
# ===================================================================
cfg <- CFG_FULL
cat(sprintf("\n-- E1/E2: M=%d, %d reps x 2 regimes --\n",
            cfg$M, cfg$n_rep_e12))
t0 <- proc.time()
jobs <- list()
for (rg in names(SIGNAL_LEVELS))
  for (r in seq_len(cfg$n_rep_e12))
    jobs[[length(jobs) + 1L]] <- list(regime = rg, rep_id = r)
e12 <- bc_run_jobs(jobs, wrap_try(run_e12_replicate), cfg, NCORES)
t_e12 <- (proc.time() - t0)[3]
saveRDS(list(results = e12, cfg = cfg, time_s = t_e12),
        file.path(RES_DIR, "e12_results.rds"))
cat(sprintf("   [%.1f s] %d/%d replicates OK\n", t_e12,
            sum(vapply(e12, function(x) is.null(x$error), logical(1))),
            length(e12)))

# ===================================================================
#  E3 + E4: Hessian analysis + Kantorovich diagnostic
# ===================================================================
cat(sprintf("\n-- E3/E4: %d reps x 2 regimes (d = %d) --\n",
            cfg$n_rep_e34, cfg$M * (cfg$K - 1L) + cfg$K * cfg$N))
t0 <- proc.time()
jobs <- list()
for (rg in names(SIGNAL_LEVELS))
  for (r in seq_len(cfg$n_rep_e34))
    jobs[[length(jobs) + 1L]] <- list(regime = rg, rep_id = r)
e34 <- bc_run_jobs(jobs, wrap_try(run_e34_replicate), cfg, NCORES)
t_e34 <- (proc.time() - t0)[3]
saveRDS(list(results = e34, cfg = cfg, time_s = t_e34),
        file.path(RES_DIR, "e34_results.rds"))
cat(sprintf("   [%.1f s] %d/%d replicates OK\n", t_e34,
            sum(vapply(e34, function(x) is.null(x$error), logical(1))),
            length(e34)))

# ===================================================================
#  E5: finite-L bias floor (optional)
# ===================================================================
if (RUN_E5) {
  DL_VALUES <- c(50L, 200L, 1000L)
  cat(sprintf("\n-- E5: doc_length in {%s}, weak regime, %d reps --\n",
              paste(DL_VALUES, collapse = ", "), cfg$n_rep_e5))
  t0 <- proc.time()
  jobs <- list()
  for (di in seq_along(DL_VALUES))
    for (r in seq_len(cfg$n_rep_e5))
      jobs[[length(jobs) + 1L]] <- list(dl = DL_VALUES[di], dl_idx = di,
                                        rep_id = r)
  e5 <- bc_run_jobs(jobs, wrap_try(run_e5_replicate), cfg, NCORES)
  t_e5 <- (proc.time() - t0)[3]
  saveRDS(list(results = e5, cfg = cfg, time_s = t_e5),
          file.path(RES_DIR, "e5_results.rds"))
  cat(sprintf("   [%.1f s] %d/%d replicates OK\n", t_e5,
              sum(vapply(e5, function(x) is.null(x$error), logical(1))),
              length(e5)))
}

cat("\n=== All experiments complete ===\n")
cat("Results in", RES_DIR, "\n")
cat("Build tables/figures with: Rscript replication/basin_check/03_report.R\n")
