#!/usr/bin/env Rscript
# ===================================================================
#  Task C — operating-regime map (exploratory, no registered
#  prediction): does the pilot subspace or the per-document Newton
#  step break as compositions approach the simplex boundary?
#
#  b_max in {0.15, 0.50, 1.00, 1.50}, M = 1000, 10 replicates,
#  seeds 40000 + bmax_index*1000 + rep.  Per cell: GL-aligned pilot
#  MSE (oracle subspace diagnostic), refined (k = 5, lambda = 0) MSE
#  under the paper metric, and GN health diagnostics (max Levenberg
#  damping, Armijo failure count, monotonicity, share of documents
#  with max_k theta_ik > 0.99 at the truth and at the endpoint).
#
#  Usage: Rscript replication/audit_block1_stm/03_regime_map.R [--ncores N]
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

map_setup <- function(root) {
  source(file.path(root, "R", "sgscatm_fit.R"))
  source(file.path(root, "R", "ilr_contrast.R"))
  source(file.path(root, "R", "utils.R"))
  source(file.path(root, "replication", "simulation", "sim_dgp.R"))
  source(file.path(root, "replication", "simulation", "sim_utils.R"))
  source(file.path(root, "replication", "basin_check", "01_functions.R"))
  invisible(NULL)
}
map_setup(ROOT)

BMAX_VALUES <- c(0.15, 0.50, 1.00, 1.50)
K_TOPICS <- 5L; P_COV <- 3L; N_VOCAB <- 500L

run_c_replicate <- function(b_max, b_idx, rep_id) {
  seed <- 40000L + b_idx * 1000L + rep_id
  dat <- sim_dgp(M = 1000L, N = N_VOCAB, K = K_TOPICS, P = P_COV,
                 b_max = b_max, sigma_eps = 0.3, alpha_beta = 0.1,
                 doc_length = 200L, seed = seed)
  Wf <- dat$W / rowSums(dat$W)
  V  <- dat$V
  Mn <- nrow(Wf)

  sat_true <- mean(apply(dat$Theta_true, 1L, max) > 0.99)

  fit <- sgscatm(dat$W, dat$C, K = K_TOPICS, lambda = 1, rotate = TRUE)
  gl  <- bc_gl_align(fit$Z, dat$Z_true)
  pilot_gl_mse    <- bc_mse_direct(gl$Z, dat$C, dat$Bz0)
  pilot_paper_mse <- procrustes_align(fit$Bz, dat$Bz0)$mse

  # 5-sweep refinement via the exported blocks, tracking Levenberg nu
  Z <- gl$Z; Phi <- bc_phi_step(Z, Wf, V)
  nu <- rep(1e-6, Mn)
  F_prev <- bc_objective(Z, Phi, Wf, V)
  n_fail <- 0L; nu_max <- 1e-6; mono_ok <- TRUE
  for (s in seq_len(5L)) {
    zs <- bc_z_step(Z, Phi, Wf, V, lambda = 0, CB = NULL, nu = nu,
                    n_gn = 2L)
    Z <- zs$Z; nu <- zs$nu
    n_fail <- n_fail + zs$n_fail
    nu_max <- max(nu_max, nu)
    Phi <- bc_phi_step(Z, Wf, V)
    F_cur <- bc_objective(Z, Phi, Wf, V)
    if (F_cur > F_prev + 1e-12 * (1 + abs(F_prev))) mono_ok <- FALSE
    F_prev <- F_cur
  }
  B_ref <- bc_b_step(Z, dat$C)
  sat_end <- mean(apply(bc_theta(Z, V), 1L, max) > 0.99)

  list(b_max = b_max, rep = rep_id, seed = seed,
       sat_true = sat_true, sat_end = sat_end,
       pilot_gl_mse = pilot_gl_mse, pilot_paper_mse = pilot_paper_mse,
       refined_paper_mse = procrustes_align(B_ref, dat$Bz0)$mse,
       n_fail = n_fail, nu_max = nu_max, monotone_ok = mono_ok)
}

wrap_try <- function(fun) {
  force(fun)
  function(...) tryCatch(fun(...), error = function(e)
    list(error = conditionMessage(e), args = list(...)))
}

cat("=== Task C: operating-regime map ===\n")
t0 <- proc.time()
jobs <- list()
for (bi in seq_along(BMAX_VALUES))
  for (r in seq_len(10L))
    jobs[[length(jobs) + 1L]] <- list(b_max = BMAX_VALUES[bi], b_idx = bi,
                                      rep_id = r)
if (NCORES > 1L) {
  cl <- parallel::makePSOCKcluster(min(NCORES, length(jobs)))
  parallel::clusterCall(cl, map_setup, root = ROOT)
  parallel::clusterExport(cl, c("K_TOPICS", "P_COV", "N_VOCAB",
                                "run_c_replicate"), envir = globalenv())
  FUN <- wrap_try(run_c_replicate)
  wf <- function(j) do.call(FUN, j)
  environment(wf) <- list2env(list(FUN = FUN), parent = globalenv())
  cres <- parallel::parLapplyLB(cl, jobs, wf)
  parallel::stopCluster(cl)
} else {
  cres <- lapply(jobs, function(j) do.call(wrap_try(run_c_replicate), j))
}
t_c <- (proc.time() - t0)[3]
saveRDS(list(results = cres, time_s = t_c),
        file.path(RES_DIR, "c_results.rds"))
cat(sprintf("[%.1f s] %d/%d replicates OK\n", t_c,
            sum(vapply(cres, function(x) is.null(x$error), logical(1))),
            length(cres)))
