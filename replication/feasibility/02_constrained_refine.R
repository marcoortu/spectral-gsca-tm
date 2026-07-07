#' ===================================================================
#'  Feasibility round — F1b: simplex-constrained refinement and the
#'  shared sweep loop with the F3 adaptive k rule
#' ===================================================================
#'
#'  The Phi-step becomes projected gradient descent on
#'  ||Wf - Theta(Z) Phi||_F^2 with each row of Phi constrained to the
#'  probability simplex (Duchi projection, step 1/Lipschitz, halving
#'  safeguard) — this removes the Z-vs-Phi scale/gauge freedom that the
#'  audit showed was silently supplied by the oracle alignment.
#'  Z-steps are unchanged (bc_z_step, lambda = 0).
#'
#'  fs_refine_rule() is the single sweep loop used by F1c/F2/F3/F4:
#'  constrained or unconstrained Phi-step, optional per-sweep B path,
#'  and the adaptive k rule
#'      stop when  max|B_s - B_{s-1}| / rms(B_{s-1}) < rule_tol
#'      for `patience` consecutive sweeps; cap max_sweeps.
#' ===================================================================

#' Simplex-constrained Phi-step: FISTA (accelerated projected gradient
#' with function-value restart) on the compressed quadratic form
#' (Theta enters only through G = Theta'Theta and TW = Theta'Wf, so an
#' iteration costs O(K^2 N) — essentially free).  Plain PGD with 10
#' iterations was tried first and freezes the Phi block (cond(G) is
#' 1e2-1e3, so 10 steps of 1/Lip move Phi microscopically; measured:
#' F stalls and even the oracle-start constrained run degraded to
#' mse 0.067).  Monotonicity is enforced by the restart: the returned
#' Phi never has a larger blockwise objective than the input.
fs_phi_step_proj <- function(Z, Phi, Wf, V, iters = 300L) {
  Th <- bc_theta(Z, V)
  G  <- crossprod(Th)
  TW <- crossprod(Th, Wf)
  lip <- 2 * max(eigen(G, symmetric = TRUE, only.values = TRUE)$values)
  qf <- function(P) sum(P * (G %*% P)) - 2 * sum(P * TW)
  f_best <- qf(Phi); P_best <- Phi
  Y <- Phi; P_prev <- Phi; t_prev <- 1
  for (it in seq_len(iters)) {
    Gr <- 2 * (G %*% Y - TW)
    Pn <- fs_proj_simplex_rows(Y - Gr / lip)
    fn <- qf(Pn)
    if (fn < f_best) { f_best <- fn; P_best <- Pn }
    if (fn > qf(P_prev)) {                   # restart momentum
      Y <- Pn; t_prev <- 1
    } else {
      t_new <- (1 + sqrt(1 + 4 * t_prev^2)) / 2
      Y <- Pn + ((t_prev - 1) / t_new) * (Pn - P_prev)
      t_prev <- t_new
    }
    P_prev <- Pn
  }
  list(Phi = P_best, n_halved = 0L)
}

#' Per-document z from scratch (or a start) with Phi held fixed:
#' damped GN via bc_z_step (per-document Armijo + Levenberg).
fs_z_init_gn <- function(Phi, Wf, V, Z0 = NULL, n_gn = 10L) {
  M <- nrow(Wf)
  if (is.null(Z0)) Z0 <- matrix(0, M, ncol(V))
  zs <- bc_z_step(Z0, Phi, Wf, V, lambda = 0, CB = NULL,
                  nu = rep(1e-6, M), n_gn = n_gn)
  list(Z = zs$Z, n_fail = zs$n_fail)
}

#' Shared sweep loop (Z-step + Phi-step), with the adaptive k rule.
#'
#' @param constrained TRUE: fs_phi_step_proj (rows on the simplex);
#'   FALSE: bc_phi_step (unconstrained ridge LS — the audit estimator).
#' @param apply_rule FALSE runs all max_sweeps (F3 uses this with
#'   track_B_path = TRUE and evaluates every k and the rule post hoc).
#' @param Bz0 evaluation-only: if given, mse_Bz (paper metric) traced.
fs_refine_rule <- function(Z0, Phi0, Wf, C, V,
                           constrained = TRUE,
                           max_sweeps = 50L, rule_tol = 1e-3,
                           patience = 2L, n_gn = 2L, phi_iters = 10L,
                           apply_rule = TRUE, track_B_path = FALSE,
                           Bz0 = NULL) {
  Z <- Z0; Phi <- Phi0
  M <- nrow(Z); Km1 <- ncol(Z)
  nu <- rep(1e-6, M)
  F_cur <- bc_objective(Z, Phi, Wf, V)
  B_prev <- bc_b_step(Z, C)
  monotone_ok <- TRUE; n_fail <- 0L; hits <- 0L
  rule_stop <- NA_integer_
  B_path <- if (track_B_path)
    array(NA_real_, c(nrow(B_prev), Km1, max_sweeps)) else NULL
  F_path <- numeric(max_sweeps)
  mse_path <- numeric(max_sweeps)

  s_used <- 0L
  for (s in seq_len(max_sweeps)) {
    zs <- bc_z_step(Z, Phi, Wf, V, lambda = 0, CB = NULL, nu = nu,
                    n_gn = n_gn)
    Z <- zs$Z; nu <- zs$nu; n_fail <- n_fail + zs$n_fail
    if (constrained) {
      Phi <- fs_phi_step_proj(Z, Phi, Wf, V, iters = phi_iters)$Phi
    } else {
      Phi <- bc_phi_step(Z, Wf, V)
    }
    F_new <- bc_objective(Z, Phi, Wf, V)
    if (F_new > F_cur + 1e-12 * (1 + abs(F_cur))) monotone_ok <- FALSE
    F_cur <- F_new
    F_path[s] <- F_new
    s_used <- s

    B_new <- bc_b_step(Z, C)
    if (track_B_path) B_path[, , s] <- B_new
    if (!is.null(Bz0)) mse_path[s] <- procrustes_align(B_new, Bz0)$mse
    relch <- max(abs(B_new - B_prev)) / (1e-8 + sqrt(mean(B_prev^2)))
    B_prev <- B_new
    hits <- if (relch < rule_tol) hits + 1L else 0L
    if (hits >= patience && is.na(rule_stop)) rule_stop <- s
    if (apply_rule && !is.na(rule_stop)) break
  }
  list(Z = Z, Phi = Phi, B = B_prev, sweeps = s_used,
       rule_stop = rule_stop, monotone_ok = monotone_ok,
       n_fail = n_fail, F_path = F_path[seq_len(s_used)],
       mse_path = if (!is.null(Bz0)) mse_path[seq_len(s_used)] else NULL,
       B_path = if (track_B_path) B_path[, , seq_len(s_used), drop = FALSE]
                else NULL)
}

#' Post-hoc simulation of the adaptive rule from a stored B path
#' (same definition as in fs_refine_rule).
fs_rule_from_path <- function(B_path, rule_tol = 1e-3, patience = 2L,
                              cap = 50L) {
  S <- dim(B_path)[3]
  hits <- 0L
  for (s in 2:min(S, cap)) {
    relch <- max(abs(B_path[, , s] - B_path[, , s - 1])) /
      (1e-8 + sqrt(mean(B_path[, , s - 1]^2)))
    hits <- if (relch < rule_tol) hits + 1L else 0L
    if (hits >= patience) return(s)
  }
  min(S, cap)
}
