#' Exact-objective engine and k-step Gauss-Newton refinement (internal)
#'
#' Internal numerical core promoted verbatim (up to renaming) from the
#' validated replication code in `replication/basin_check/01_functions.R`
#' (derivative-checked against numDeriv there) and the read-out /
#' B-stationarity refinement of `replication/feasibility/`. These routines
#' implement the exact (non-linearised) reconstruction objective
#' \deqn{F_0(Z,\Phi)=\lVert \tilde W - \Theta(Z)\,\Phi\rVert_F^2,\qquad
#'       \Theta(Z)=\mathrm{softmax\_rows}(ZV^\top),}
#' its analytic gradient, a damped per-document Gauss-Newton z-step, and the
#' least-squares \eqn{\Phi}- and B-steps used by [sgscatm_chain()].
#'
#' All functions are non-exported (prefixed `.sg_`) and touch no other package
#' internals. The refinement z-step criterion is **least squares**, never the
#' multinomial deviance (a settled design choice; see the manuscript Lemma 20
#' and `replication/deviance_probe/`).
#' @keywords internal
#' @name sgscatm-refine-internal
NULL

# Row-wise softmax of Z V' (numerically stable)
.sg_theta <- function(Z, V) {
  S <- tcrossprod(Z, V)
  S <- S - apply(S, 1L, max)
  E <- exp(S)
  E / rowSums(E)
}

# Residual of the OLS projection of Z on C (profiled penalty term)
.sg_z_resid <- function(Z, C, ridge = 1e-8) {
  B <- solve(crossprod(C) + ridge * diag(ncol(C)), crossprod(C, Z))
  Z - C %*% B
}

# Profiled exact objective value
.sg_objective <- function(Z, Phi, Wf, V, lambda = 0, C = NULL) {
  R <- .sg_theta(Z, V) %*% Phi - Wf
  f <- sum(R * R)
  if (lambda > 0) {
    Zr <- .sg_z_resid(Z, C)
    f  <- f + lambda * sum(Zr * Zr)
  }
  f
}

# Profiled objective + analytic gradient (Gz per-document, Gphi pooled)
.sg_grad <- function(Z, Phi, Wf, V, lambda = 0, C = NULL) {
  Theta <- .sg_theta(Z, V)
  R  <- Theta %*% Phi - Wf
  f  <- sum(R * R)
  U  <- R %*% t(Phi)
  TU <- Theta * U
  Gz   <- 2 * (TU - Theta * rowSums(TU)) %*% V
  Gphi <- 2 * crossprod(Theta, R)
  if (lambda > 0) {
    Zr <- .sg_z_resid(Z, C)
    f  <- f + lambda * sum(Zr * Zr)
    Gz <- Gz + 2 * lambda * Zr
  }
  list(F = f, Gz = Gz, Gphi = Gphi)
}

# Phi-step: ridge-stabilised least squares given Z
.sg_phi_step <- function(Z, Wf, V, ridge = 1e-8) {
  Th <- .sg_theta(Z, V)
  solve(crossprod(Th) + ridge * diag(ncol(Th)), crossprod(Th, Wf))
}

# B-step: closed-form OLS of Z on C
.sg_b_step <- function(Z, C, ridge = 1e-8) {
  solve(crossprod(C) + ridge * diag(ncol(C)), crossprod(C, Z))
}

# Z-step: n_gn damped Gauss-Newton iterations per document (least squares),
# Armijo backtracking + Levenberg damping, optional trust-region cap on ||dz||.
# lambda = 0 for the refinement (CB = NULL); lambda > 0 supported for the pilot
# penalty via CB = C %*% B. Verbatim from bc_z_step with an added dz_cap guard.
.sg_z_step <- function(Z, Phi, Wf, V, lambda = 0, CB = NULL, nu,
                       n_gn = 2L, max_bt = 30L, ridge = 1e-8, dz_cap = Inf) {
  M   <- nrow(Z); Km1 <- ncol(Z)
  A   <- tcrossprod(Phi)
  Pw  <- Wf %*% t(Phi)
  w2  <- rowSums(Wf * Wf)

  fdoc <- function(Zm, rows) {
    Th <- .sg_theta(Zm, V)
    f  <- rowSums((Th %*% A) * Th) - 2 * rowSums(Th * Pw[rows, , drop = FALSE]) +
          w2[rows]
    if (lambda > 0)
      f <- f + lambda * rowSums((Zm - CB[rows, , drop = FALSE])^2)
    f
  }

  n_fail_tot <- 0L
  for (it in seq_len(n_gn)) {
    Th  <- .sg_theta(Z, V)
    ThA <- Th %*% A
    Uu  <- ThA - Pw
    TU  <- Th * Uu
    G   <- 2 * (TU - Th * rowSums(TU)) %*% V
    f0  <- rowSums(ThA * Th) - 2 * rowSums(Th * Pw) + w2
    if (lambda > 0) {
      G  <- G + 2 * lambda * (Z - CB)
      f0 <- f0 + lambda * rowSums((Z - CB)^2)
    }

    S  <- Th %*% V
    Ba <- vector("list", Km1)
    for (a in seq_len(Km1))
      Ba[[a]] <- Th * (matrix(V[, a], M, ncol(Th), byrow = TRUE) - S[, a])
    Ca <- lapply(Ba, function(B) B %*% A)
    Harr <- array(0, c(M, Km1, Km1))
    for (a in seq_len(Km1)) for (b in seq_len(Km1))
      Harr[, a, b] <- rowSums(Ba[[a]] * Ca[[b]])

    Delta <- matrix(0, M, Km1)
    for (i in seq_len(M)) {
      Hi <- Harr[i, , ] + diag(lambda + nu[i] + ridge, Km1)
      Delta[i, ] <- tryCatch(-solve(Hi, G[i, ] / 2),
                             error = function(e) rep(0, Km1))
    }
    # trust-region cap: rescale any step longer than dz_cap
    if (is.finite(dz_cap)) {
      nrm <- sqrt(rowSums(Delta^2))
      big <- nrm > dz_cap
      if (any(big)) Delta[big, ] <- Delta[big, , drop = FALSE] *
        (dz_cap / nrm[big])
    }

    gd   <- rowSums(G * Delta)
    Delta[gd >= 0, ] <- 0
    gd   <- pmin(gd, 0)
    step <- rep(1, M)
    Znew <- Z + Delta
    fnew <- fdoc(Znew, seq_len(M))
    acc  <- fnew <= f0 + 1e-4 * step * gd
    bt   <- 0L
    while (any(!acc) && bt < max_bt) {
      bt  <- bt + 1L
      idx <- which(!acc)
      step[idx] <- step[idx] / 2
      Znew[idx, ] <- Z[idx, , drop = FALSE] +
        step[idx] * Delta[idx, , drop = FALSE]
      fnew[idx] <- fdoc(Znew[idx, , drop = FALSE], idx)
      acc[idx]  <- fnew[idx] <= f0[idx] + 1e-4 * step[idx] * gd[idx]
    }
    fail <- !acc & rowSums(Delta != 0) > 0
    if (any(fail)) {
      Znew[fail, ] <- Z[fail, , drop = FALSE]
      nu[fail]     <- nu[fail] * 10
      n_fail_tot   <- n_fail_tot + sum(fail)
    }
    nu[acc] <- pmax(nu[acc] / 3, 1e-6)
    Z <- Znew
  }
  list(Z = Z, nu = nu, n_fail = n_fail_tot)
}

# Per-document read-out from scratch (or a warm start) with Phi held fixed:
# damped GN via .sg_z_step (trust-region capped). Promoted from fs_z_init_gn.
.sg_readout_gn <- function(Phi, Wf, V, Z0 = NULL, n_gn = 10L, dz_cap = 1) {
  M <- nrow(Wf)
  if (is.null(Z0)) Z0 <- matrix(0, M, ncol(V))
  zs <- .sg_z_step(Z0, Phi, Wf, V, lambda = 0, CB = NULL,
                   nu = rep(1e-6, M), n_gn = n_gn, dz_cap = dz_cap)
  list(Z = zs$Z, n_fail = zs$n_fail)
}

# k-step refinement with the B-stationarity stop.
#   mode = "frozen_phi": Z-only sweeps, Phi held at Phi0 (theory baseline,
#          Section 3.7 / Lemma 17 regime).
#   mode = "joint":      Z-step then unconstrained ridge-LS Phi-step each
#          sweep (the empirically-best V4 variant; Prop 19 alternating).
# Stop: max_{j,a}|dB|/rms(B) < rule_tol for `patience` consecutive sweeps.
.sg_refine <- function(Z0, Phi0, Wf, C, V,
                       mode = c("frozen_phi", "joint"),
                       max_sweeps = 100L, rule_tol = 1e-3, patience = 2L,
                       n_gn = 2L, dz_cap = 1) {
  mode <- match.arg(mode)
  Z <- Z0; Phi <- Phi0
  M <- nrow(Z)
  nu <- rep(1e-6, M)
  F_cur <- .sg_objective(Z, Phi, Wf, V)
  B_prev <- .sg_b_step(Z, C)
  monotone_ok <- TRUE; n_fail <- 0L; hits <- 0L; f_hits <- 0L
  rule_stop <- NA_integer_; s_used <- 0L
  F_path <- numeric(max_sweeps)
  for (s in seq_len(max_sweeps)) {
    zs <- .sg_z_step(Z, Phi, Wf, V, lambda = 0, CB = NULL, nu = nu,
                     n_gn = n_gn, dz_cap = dz_cap)
    Z <- zs$Z; nu <- zs$nu; n_fail <- n_fail + zs$n_fail
    if (mode == "joint") Phi <- .sg_phi_step(Z, Wf, V)
    F_new <- .sg_objective(Z, Phi, Wf, V)
    if (F_new > F_cur + 1e-12 * (1 + abs(F_cur))) monotone_ok <- FALSE
    rel_f <- (F_cur - F_new) / (1 + abs(F_cur))     # F-convergence backstop
    F_cur <- F_new; F_path[s] <- F_new; s_used <- s
    B_new <- .sg_b_step(Z, C)
    relch <- max(abs(B_new - B_prev)) / (1e-8 + sqrt(mean(B_prev^2)))
    B_prev <- B_new
    # primary stop: B-stationarity. Backstop (fires cleanly at small M where B
    # drifts within numerical noise but the objective has converged): relative
    # objective decrease below tol_f for `patience` consecutive sweeps.
    hits   <- if (relch < rule_tol) hits + 1L else 0L
    f_hits <- if (rel_f < 1e-9) f_hits + 1L else 0L
    if ((hits >= patience || f_hits >= patience) && is.na(rule_stop))
      rule_stop <- s
    if (!is.na(rule_stop)) break
  }
  list(Z = Z, Phi = Phi, B = B_prev, sweeps = s_used, rule_stop = rule_stop,
       monotone_ok = monotone_ok, n_fail = n_fail,
       F_path = F_path[seq_len(s_used)])
}
