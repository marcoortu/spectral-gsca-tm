#!/usr/bin/env Rscript
# ===================================================================
#  Spot-check runner — SC0 (unit gates), SC-A (bias field b(z)),
#  SC-B (gradient direction G0, two routes), SC-C (optional cosine)
# ===================================================================
#
#  Usage (from the package root):
#    Rscript replication/spotcheck/02_run.R [--quick] [--ncores N]
#
#  Seeds: Phi0-A 77001, Phi0-B 77002, test directions 77010,
#  SC-A multinomial draws 77000 + point_idx*100 + L_idx, brute cell
#  77777 (+ chunk), SC-B z-law/multinomial base draws 78200 + L,
#  SC-B directions 78100, SC-B route-1 78000 + regime index, SC0 77000.
#  SC-C reuses the feasibility E2 seeds (90000 + 2000 + rep).
# ===================================================================

args   <- commandArgs(trailingOnly = TRUE)
QUICK  <- "--quick" %in% args
ONLY_A <- "--only-a" %in% args   # rerun SC-A/brute only (certificate fix)
NCORES <- if ("--ncores" %in% args) {
  as.integer(args[which(args == "--ncores") + 1L])
} else {
  max(1L, min(10L, parallel::detectCores() - 2L))
}

ROOT <- Sys.getenv("SGSCATM_ROOT", unset = normalizePath(".", mustWork = FALSE))
if (!dir.exists(file.path(ROOT, "R")))
  stop("Run from the package root or set SGSCATM_ROOT.")
RES_DIR <- file.path(ROOT, "replication", "spotcheck", "results")
dir.create(RES_DIR, recursive = TRUE, showWarnings = FALSE)

sc_setup <- function(root) {
  source(file.path(root, "R", "sgscatm_fit.R"))
  source(file.path(root, "R", "ilr_contrast.R"))
  source(file.path(root, "R", "utils.R"))
  source(file.path(root, "replication", "simulation", "sim_dgp.R"))
  source(file.path(root, "replication", "simulation", "sim_utils.R"))
  source(file.path(root, "replication", "basin_check", "01_functions.R"))
  source(file.path(root, "replication", "feasibility", "01_anchors.R"))
  source(file.path(root, "replication", "feasibility",
                   "02_constrained_refine.R"))
  source(file.path(root, "replication", "spotcheck", "01_formulas.R"))
  invisible(NULL)
}
sc_setup(ROOT)

R_REP  <- if (QUICK) 2000L else 20000L
L_GRID <- c(100L, 200L, 400L, 800L)
K_W <- 5L; N_W <- 500L                       # working design

sc_phi0 <- function(seed, K = K_W, N = N_W, alpha = 0.1) {
  set.seed(seed)
  .rdirichlet_matrix(K, N, alpha)
}

#' deterministic test-point table (9 points)
sc_points <- function() {
  set.seed(77010L)
  dirs <- lapply(1:6, function(i) {
    u <- rnorm(K_W - 1L); u / sqrt(sum(u^2))
  })
  norms <- c(0.5, 0.5, 1.0, 1.0, 1.5, 1.5)
  pts <- list(list(label = "A_z0", phi_seed = 77001L,
                   z = rep(0, K_W - 1L)))
  for (i in 1:6)
    pts[[length(pts) + 1L]] <- list(
      label = sprintf("A_u%d_n%.1f", i, norms[i]), phi_seed = 77001L,
      z = dirs[[i]] * norms[i])
  pts[[length(pts) + 1L]] <- list(label = "B_z0", phi_seed = 77002L,
                                  z = rep(0, K_W - 1L))
  pts[[length(pts) + 1L]] <- list(label = "B_u3_n1.0", phi_seed = 77002L,
                                  z = dirs[[3]] * 1.0)
  pts
}

# ===================================================================
#  SC-A replicate job: one (test point, L) cell
# ===================================================================
run_a_cell <- function(point_idx, L_idx) {
  pts <- sc_points()
  pt  <- pts[[point_idx]]
  L   <- L_GRID[L_idx]
  Phi0 <- sc_phi0(pt$phi_seed)
  V <- ilr_contrast(K_W)
  ob <- sc_objects(pt$z, Phi0, V, want_G0 = FALSE)

  set.seed(77000L + point_idx * 100L + L_idx)
  n <- t(rmultinom(R_REP, L, ob$p))
  Wf <- n / L
  gn <- sc_gn_batch(Wf, Phi0, V, pt$z, tol = 1e-8)
  ok <- gn$ok
  E  <- Wf[ok, , drop = FALSE] -
    matrix(ob$p, sum(ok), N_W, byrow = TRUE)
  Mlin <- E %*% ob$J %*% ob$Hi                # rows = m(e_r)'
  Dev  <- gn$Z[ok, , drop = FALSE] -
    matrix(pt$z, sum(ok), K_W - 1L, byrow = TRUE) - Mlin
  bhat <- L * colMeans(Dev)
  bse  <- L * apply(Dev, 2L, sd) / sqrt(sum(ok))
  list(point = pt$label, point_idx = point_idx, L = L,
       bhat = bhat, bse = bse, b_closed = ob$b,
       b1_closed = ob$b1, a2half = as.vector(-0.5 * ob$Hi %*% ob$a2),
       n_fail = gn$n_fail, R_used = sum(ok))
}

#' brute-force cell chunk (K = 3, N = 30, L = 50, no control variate)
run_brute_chunk <- function(chunk_id, R_chunk = 20000L) {
  K <- 3L; N <- 30L; L <- 50L
  V <- ilr_contrast(K)
  set.seed(77777L)
  Phi0 <- .rdirichlet_matrix(K, N, 0.3)
  u <- rnorm(K - 1L); z <- u / sqrt(sum(u^2))
  ob <- sc_objects(z, Phi0, V, want_G0 = FALSE)
  set.seed(77700L + chunk_id)
  n <- t(rmultinom(R_chunk, L, ob$p))
  Wf <- n / L
  gn <- sc_gn_batch(Wf, Phi0, V, z, tol = 1e-8)
  ok <- gn$ok
  E <- Wf[ok, , drop = FALSE] - matrix(ob$p, sum(ok), N, byrow = TRUE)
  Mlin <- E %*% ob$J %*% ob$Hi
  Dz <- gn$Z[ok, , drop = FALSE] - matrix(z, sum(ok), K - 1L, byrow = TRUE)
  list(chunk = chunk_id, R_used = sum(ok), n_fail = gn$n_fail,
       sum_raw = colSums(Dz), sum_cv = colSums(Dz - Mlin),
       ss_raw = colSums(Dz^2), ss_cv = colSums((Dz - Mlin)^2),
       b_closed = ob$b, L = L, z = z)
}

# ===================================================================
#  SC-B route 1: closed-form average of G0 under the Block-3 z law
# ===================================================================
sc_draw_z <- function(S, b_max, seed) {
  set.seed(seed)
  t(vapply(seq_len(S), function(s) {
    B <- matrix(runif(3L * (K_W - 1L), -b_max, b_max), 3L, K_W - 1L)
    as.vector(crossprod(B, rnorm(3L))) + rnorm(K_W - 1L, 0, 0.3)
  }, numeric(K_W - 1L)))
}

run_b1_regime <- function(regime, S = 5000L) {
  b_max <- c(weak = 0.15, strong = 0.50)[[regime]]
  Phi0 <- sc_phi0(77001L)
  V <- ilr_contrast(K_W)
  Zs <- sc_draw_z(S, b_max, 78000L + match(regime, c("weak", "strong")))
  nb <- 10L
  batch <- rep(seq_len(nb), length.out = S)
  Gsum <- matrix(0, K_W, N_W)
  Gb <- lapply(seq_len(nb), function(i) matrix(0, K_W, N_W))
  cnt <- integer(nb)
  psum <- numeric(N_W)
  for (s in seq_len(S)) {
    ob <- sc_objects(Zs[s, ], Phi0, V, want_G0 = TRUE)
    Gsum <- Gsum + ob$G0
    Gb[[batch[s]]] <- Gb[[batch[s]]] + ob$G0
    cnt[batch[s]] <- cnt[batch[s]] + 1L
    psum <- psum + ob$p
  }
  Gbar <- Gsum / S
  Ustar <- Gbar / sqrt(sum(Gbar^2))
  proj <- vapply(seq_len(nb), function(i)
    sum((Gb[[i]] / cnt[i]) * Ustar), numeric(1))
  list(regime = regime, Gbar = Gbar,
       fro = sqrt(sum(Gbar^2)), fro_se = sd(proj) / sqrt(nb),
       pbar = psum / S)
}

# ===================================================================
#  SC-B route 2: finite differences with common random numbers
# ===================================================================
sc_b2_dirs <- function(pbar) {
  set.seed(78100L)
  U <- lapply(1:4, function(i) {
    u <- matrix(rnorm(K_W * N_W), K_W, N_W); u / sqrt(sum(u^2))
  })
  low <- order(pbar)[1:50]
  for (i in 1:2) {
    u <- matrix(0, K_W, N_W)
    u[, low] <- rnorm(K_W * 50L)
    U[[4L + i]] <- u / sqrt(sum(u^2))
  }
  U
}

sc_b2_base <- function(L, R = R_REP) {
  Phi0 <- sc_phi0(77001L)
  V <- ilr_contrast(K_W)
  Zs <- sc_draw_z(R, 0.50, 78200L + L)             # strong regime
  Th <- bc_theta(Zs, V)
  P0 <- Th %*% Phi0
  set.seed(78300L + L)
  n <- matrix(0L, R, N_W)
  for (r in seq_len(R)) n[r, ] <- rmultinom(1L, L, P0[r, ])
  list(Zs = Zs, Wf = n / L, Phi0 = Phi0, V = V)
}

run_b2_cell <- function(L, dir_idx, h, sgn, pbar) {
  bs <- sc_b2_base(L)
  U  <- sc_b2_dirs(pbar)[[dir_idx]]
  Phi <- bs$Phi0 + sgn * h * U
  gn <- sc_gn_batch_z0(bs$Wf, Phi, bs$V, bs$Zs, tol = 1e-10)
  Th <- bc_theta(gn$Z, bs$V)
  loss <- rowSums((bs$Wf - Th %*% Phi)^2)
  list(L = L, dir = dir_idx, h = h, sgn = sgn, loss = loss,
       ok = gn$ok, n_fail = gn$n_fail)
}

run_b2_g0 <- function(L) {
  bs <- sc_b2_base(L)
  R <- nrow(bs$Zs)
  Gsum <- matrix(0, K_W, N_W)
  for (r in seq_len(R))
    Gsum <- Gsum + sc_objects(bs$Zs[r, ], bs$Phi0, bs$V,
                              want_G0 = TRUE)$G0
  list(L = L, Gbar = Gsum / R)
}

#' batch GN with per-row start matrix (SC-B: each row starts at its z_r)
sc_gn_batch_z0 <- function(Wf, Phi, V, Z0, tol = 1e-10, max_calls = 60L) {
  Z <- Z0; nu <- rep(1e-6, nrow(Wf))
  for (it in seq_len(max_calls)) {
    gm <- apply(abs(sc_grad_batch(Z, Phi, Wf, V)), 1L, max)
    if (max(gm) < tol) break
    zs <- bc_z_step(Z, Phi, Wf, V, lambda = 0, CB = NULL, nu = nu,
                    n_gn = 2L)
    Z <- zs$Z; nu <- zs$nu
  }
  gm <- apply(abs(sc_grad_batch(Z, Phi, Wf, V)), 1L, max)
  list(Z = Z, gmax = gm, ok = gm < tol, n_fail = sum(gm >= tol))
}

# ===================================================================
#  SC-C: cosine of -G0 vs the early LS-pathology Phi displacement
#  (Phi paths were not stored in the feasibility RDS; recomputed on
#  the same seeds, 10 constrained-LS sweeps from the truth)
# ===================================================================
run_c_replicate <- function(rep_id) {
  seed <- 90000L + 2L * 1000L + rep_id                 # strong regime
  dat <- sim_dgp(M = 1000L, N = N_W, K = K_W, P = 3L, b_max = 0.50,
                 sigma_eps = 0.3, alpha_beta = 0.1, doc_length = 200L,
                 seed = seed)
  Wf <- dat$W / rowSums(dat$W)
  Z <- dat$Z_true; Phi <- dat$Beta
  nu <- rep(1e-6, 1000L)
  for (s in 1:10) {
    zs <- bc_z_step(Z, Phi, Wf, dat$V, lambda = 0, CB = NULL, nu = nu,
                    n_gn = 2L)
    Z <- zs$Z; nu <- zs$nu
    Phi <- fs_phi_step_proj(Z, Phi, Wf, dat$V)$Phi
  }
  list(rep = rep_id, dPhi = Phi - dat$Beta)
}

# ===================================================================
#  Driver
# ===================================================================
wrap_try <- function(fun) {
  force(fun)
  function(...) tryCatch(fun(...), error = function(e)
    list(error = conditionMessage(e), args = list(...)))
}
sc_run_jobs <- function(jobs, fun, ncores) {
  if (ncores <= 1L || length(jobs) <= 1L)
    return(lapply(jobs, function(j) do.call(fun, j)))
  cl <- parallel::makePSOCKcluster(min(ncores, length(jobs)))
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(cl, sc_setup, root = ROOT)
  parallel::clusterExport(cl, c("R_REP", "L_GRID", "K_W", "N_W",
                                "sc_phi0", "sc_points", "sc_draw_z",
                                "sc_b2_dirs", "sc_b2_base",
                                "sc_gn_batch_z0",
                                "run_a_cell", "run_brute_chunk",
                                "run_b1_regime", "run_b2_cell",
                                "run_b2_g0", "run_c_replicate"),
                          envir = globalenv())
  wf <- function(j) do.call(FUN, j)
  environment(wf) <- list2env(list(FUN = fun), parent = globalenv())
  parallel::parLapplyLB(cl, jobs, wf)
}
report_ok <- function(res, t_el) {
  n_ok <- sum(vapply(res, function(x) is.null(x$error), logical(1)))
  cat(sprintf("   [%.1f s] %d/%d jobs OK\n", t_el, n_ok, length(res)))
  if (n_ok < length(res))
    cat("   ERRORS:", paste(head(unique(unlist(lapply(res, function(x)
      x$error))), 3), collapse = " | "), "\n")
}

cat("=== Section-4 spot-check ===\n")
cat("Mode:", if (QUICK) "QUICK" else "FULL", "| cores:", NCORES, "\n")

cat("\n-- SC0 unit tests (K=3, N=20) --\n")
t0 <- proc.time()
ut <- sc_verify()
cat(sprintf("   H2theta %.1e | H2f %.1e | opt-vs-naive %.1e | GN gmax %.1e\n",
            ut$H2theta_relerr, ut$H2f_relerr, ut$opt_vs_naive, ut$gn_gmax))
saveRDS(ut, file.path(RES_DIR, "sc0_unit_tests.rds"))
cat(sprintf("   [%.1f s] PASSED\n", (proc.time() - t0)[3]))

cat("\n-- SC-A: bias field, 9 points x 4 L, R =", R_REP, "--\n")
t0 <- proc.time()
jobs <- list()
for (p in seq_along(sc_points()))
  for (li in seq_along(L_GRID))
    jobs[[length(jobs) + 1L]] <- list(point_idx = p, L_idx = li)
scA <- sc_run_jobs(jobs, wrap_try(run_a_cell), NCORES)
t_a <- (proc.time() - t0)[3]
saveRDS(list(results = scA, time_s = t_a),
        file.path(RES_DIR, "sca_results.rds"))
report_ok(scA, t_a)

cat("\n-- SC-A brute cell (K=3, N=30, L=50, R = 300k, no CV) --\n")
t0 <- proc.time()
n_chunks <- if (QUICK) 2L else 15L
scAb <- sc_run_jobs(lapply(seq_len(n_chunks), function(i)
  list(chunk_id = i)), wrap_try(run_brute_chunk), NCORES)
t_ab <- (proc.time() - t0)[3]
saveRDS(list(results = scAb, time_s = t_ab),
        file.path(RES_DIR, "sca_brute_results.rds"))
report_ok(scAb, t_ab)

if (ONLY_A) {
  cat("\n--only-a: skipping SC-B and SC-C (results already on disk)\n")
  quit(save = "no", status = 0L)
}

cat("\n-- SC-B route 1: E[G0] under the Block-3 z law --\n")
t0 <- proc.time()
scB1 <- sc_run_jobs(list(list(regime = "weak"), list(regime = "strong")),
                    wrap_try(run_b1_regime), NCORES)
t_b1 <- (proc.time() - t0)[3]
saveRDS(list(results = scB1, time_s = t_b1),
        file.path(RES_DIR, "scb1_results.rds"))
report_ok(scB1, t_b1)
pbar <- scB1[[which(vapply(scB1, function(x)
  identical(x$regime, "strong"), logical(1)))]]$pbar

cat("\n-- SC-B route 2: FD with CRN (strong regime) --\n")
t0 <- proc.time()
jobs <- list()
for (L in c(200L, 400L))
  for (d in 1:6)
    for (h in c(1e-3, 1e-4))
      for (sgn in c(1, -1))
        jobs[[length(jobs) + 1L]] <- list(L = L, dir_idx = d, h = h,
                                          sgn = sgn, pbar = pbar)
scB2 <- sc_run_jobs(jobs, wrap_try(run_b2_cell), NCORES)
scB2g <- sc_run_jobs(list(list(L = 200L), list(L = 400L)),
                     wrap_try(run_b2_g0), NCORES)
t_b2 <- (proc.time() - t0)[3]
saveRDS(list(results = scB2, g0 = scB2g, time_s = t_b2),
        file.path(RES_DIR, "scb2_results.rds"))
report_ok(scB2, t_b2)

cat("\n-- SC-C: cosine vs LS-pathology Phi displacement (5 reps) --\n")
t0 <- proc.time()
scC <- sc_run_jobs(lapply(1:5, function(r) list(rep_id = r)),
                   wrap_try(run_c_replicate), NCORES)
t_c <- (proc.time() - t0)[3]
saveRDS(list(results = scC, time_s = t_c),
        file.path(RES_DIR, "scc_results.rds"))
report_ok(scC, t_c)

cat("\n=== Spot-check runs complete ===\n")
cat("Build tables with: Rscript replication/spotcheck/03_report.R\n")
