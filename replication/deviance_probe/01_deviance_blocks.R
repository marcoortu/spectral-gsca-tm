#' ===================================================================
#'  Deviance probe — multinomial-deviance refinement blocks
#' ===================================================================
#'
#'  Criterion (counts n_ij, p_i = Phi' theta(z_i), rows of Phi on the
#'  simplex, theta = softmax(V z)):
#'      F_dev(Z, Phi) = - sum_ij n_ij log p_ij .
#'  Zero-count cells contribute nothing — the deviance weighs each cell
#'  by its information instead of fitting its noise (the repair for D1).
#'
#'  Blocks (derived from scratch; unit-tested below against numDeriv):
#'    z-step  : damped Fisher scoring per document.  With
#'              D = diag(theta) - theta theta',  A = Phi' D V (N x K-1):
#'                grad_z F_dev = - V' D Phi (n/p)
#'                I(z)         = L * A' diag(1/p) A
#'              (Fisher information of the multinomial with L draws),
#'              Armijo backtracking + per-document Levenberg damping,
#'              mirroring bc_z_step.
#'    Phi-step: pLSI EM M-step given Theta —
#'                Phi_kj <- Phi_kj * (Theta' (n/p))_kj, rows renormal-
#'              ised — monotone in F_dev for fixed Theta.
#'
#'  dv_refine() is the shared sweep loop for BOTH criteria ("dev" and
#'  the constrained LS of the feasibility round) so that P1 overlays
#'  run on identical data with identical recording, including the
#'  gauge-drift decomposition against the basin_check E3 basis.
#'  All functions are prefixed dv_.
#' ===================================================================

# -------------------------------------------------------------------
#  Objective
# -------------------------------------------------------------------
dv_F <- function(Z, Phi, Wn, V) {
  P <- pmax(bc_theta(Z, V) %*% Phi, 1e-12)
  -sum(Wn * log(P))
}

# -------------------------------------------------------------------
#  z-step: damped Fisher scoring, vectorised over documents
# -------------------------------------------------------------------
dv_z_step <- function(Z, Phi, Wn, V, nu, n_gn = 2L, max_bt = 30L,
                      ridge = 1e-8) {
  M <- nrow(Z); Km1 <- ncol(Z); K <- nrow(Phi)
  L <- rowSums(Wn)

  fdoc <- function(Zm, rows) {
    P <- pmax(bc_theta(Zm, V) %*% Phi, 1e-12)
    -rowSums(Wn[rows, , drop = FALSE] * log(P))
  }

  n_fail_tot <- 0L
  for (it in seq_len(n_gn)) {
    Th <- bc_theta(Z, V)
    P  <- pmax(Th %*% Phi, 1e-12)
    Ratio <- Wn / P                          # 0 where n = 0
    U  <- Ratio %*% t(Phi)                   # rows = (Phi (n/p))'
    TU <- Th * U
    G  <- -(TU - Th * rowSums(TU)) %*% V     # M x (K-1) gradient
    f0 <- -rowSums(Wn * log(P))

    # Fisher middle matrices Mid_i = Phi diag(1/p_i) Phi'  (K x K)
    invP <- 1 / P
    Mid <- array(0, c(M, K, K))
    for (a in seq_len(K)) for (b in a:K) {
      v <- invP %*% (Phi[a, ] * Phi[b, ])
      Mid[, a, b] <- v
      if (b > a) Mid[, b, a] <- v
    }
    # B_i = D_i V, assembled as in bc_z_step
    S  <- Th %*% V
    Ba <- vector("list", Km1)
    for (a in seq_len(Km1))
      Ba[[a]] <- Th * (matrix(V[, a], M, K, byrow = TRUE) - S[, a])
    # I_i = L_i * B_i' Mid_i B_i
    Ca <- array(0, c(M, K, Km1))
    for (a in seq_len(Km1)) for (k in seq_len(K)) {
      acc <- 0
      for (l in seq_len(K)) acc <- acc + Mid[, k, l] * Ba[[a]][, l]
      Ca[, k, a] <- acc
    }
    Harr <- array(0, c(M, Km1, Km1))
    for (a in seq_len(Km1)) for (b in seq_len(Km1)) {
      acc <- 0
      for (k in seq_len(K)) acc <- acc + Ba[[a]][, k] * Ca[, k, b]
      Harr[, a, b] <- L * acc
    }

    # per-document damped Fisher solves, with a trust-region cap on the
    # step norm (radius 1 in ILR units; z-components have sd ~ 0.35).
    # Without the cap, an accepted Fisher step from a blended-Phi start
    # jumps a document straight to softmax saturation (the per-document
    # deviance MLE under a blended Phi is extreme), where gradients die
    # and the document is stranded (measured: nr 12, mse 8.7).
    Delta <- matrix(0, M, Km1)
    for (i in seq_len(M)) {
      Fi <- Harr[i, , ] + diag(nu[i] + ridge, Km1)
      Delta[i, ] <- tryCatch(-solve(Fi, G[i, ]),
                             error = function(e) rep(0, Km1))
    }
    dn <- sqrt(rowSums(Delta^2))
    over <- dn > 1
    if (any(over)) Delta[over, ] <- Delta[over, , drop = FALSE] / dn[over]

    # vectorised Armijo backtracking (c1 = 1e-4)
    gd <- rowSums(G * Delta)
    Delta[gd >= 0, ] <- 0
    gd <- pmin(gd, 0)
    step <- rep(1, M)
    Znew <- Z + Delta
    fnew <- fdoc(Znew, seq_len(M))
    acc  <- fnew <= f0 + 1e-4 * step * gd
    bt <- 0L
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
      nu[fail] <- nu[fail] * 10
      n_fail_tot <- n_fail_tot + sum(fail)
    }
    nu[acc] <- pmax(nu[acc] / 3, 1e-6)
    Z <- Znew
  }
  list(Z = Z, nu = nu, n_fail = n_fail_tot)
}

# -------------------------------------------------------------------
#  Phi-step: pLSI EM (monotone in F_dev for fixed Theta)
# -------------------------------------------------------------------
dv_phi_em <- function(Z, Phi, Wn, V, iters = 5L, floor = 1e-12) {
  Th <- bc_theta(Z, V)
  for (it in seq_len(iters)) {
    P <- pmax(Th %*% Phi, 1e-12)
    Phi <- Phi * crossprod(Th, Wn / P)
    Phi <- pmax(Phi, floor)
    Phi <- Phi / rowSums(Phi)
  }
  Phi
}

#' Per-document z from a start (or z = 0) with Phi fixed (deviance)
dv_z_init <- function(Phi, Wn, V, Z0 = NULL, n_gn = 10L) {
  M <- nrow(Wn)
  if (is.null(Z0)) Z0 <- matrix(0, M, ncol(V))
  zs <- dv_z_step(Z0, Phi, Wn, V, nu = rep(1e-6, M), n_gn = n_gn)
  list(Z = zs$Z, n_fail = zs$n_fail)
}

# -------------------------------------------------------------------
#  Shared sweep loop for both criteria, with rule + gauge tracking
# -------------------------------------------------------------------
#' @param criterion "dev": dv_z_step + dv_phi_em (Phi rows stay on the
#'   simplex by construction); "ls": bc_z_step + fs_phi_step_proj (the
#'   feasibility round's simplex-constrained LS blocks).  F_path is the
#'   criterion's own objective.
#' @param QT,eta0 optional gauge tracking: per sweep, the displacement
#'   pack(Z,Phi) - eta0 is split into ||projection on span(QT)|| and
#'   the orthogonal remainder (the "statistical tilt").
#' @param phi_first for "dev" (default TRUE): apply one EM Phi-step
#'   before the first z-step.  A start whose Phi has floored cells
#'   (e.g. the clipped oracle-GL Phi: 44% of LS entries are negative)
#'   otherwise produces p ~ 1e-8 on counted cells, and the 1/p terms
#'   catapult per-document z into softmax saturation before EM can
#'   regrow those cells (measured: mse 2-8 instead of ~0.007).  One
#'   multiplicative EM update regrows them from the data; at the truth
#'   this is a harmless first block move.
dv_refine <- function(Z0, Phi0, Wn, C, V, criterion = c("dev", "ls"),
                      max_sweeps = 100L, rule_tol = 1e-3, patience = 2L,
                      apply_rule = TRUE, n_gn = 2L, em_iters = 5L,
                      Bz0 = NULL, QT = NULL, eta0 = NULL,
                      track_B_path = FALSE, phi_first = TRUE) {
  criterion <- match.arg(criterion)
  Z <- Z0; Phi <- Phi0
  if (criterion == "dev" && phi_first)
    Phi <- dv_phi_em(Z, Phi, Wn, V, iters = 1L)
  M <- nrow(Z); Km1 <- ncol(Z)
  nu <- rep(1e-6, M)
  Wf <- if (criterion == "ls") Wn / rowSums(Wn) else NULL
  fobj <- function() if (criterion == "dev") dv_F(Z, Phi, Wn, V) else
    bc_objective(Z, Phi, Wf, V)

  F_cur <- fobj()
  B_prev <- bc_b_step(Z, C)
  monotone_ok <- TRUE; n_fail <- 0L; hits <- 0L
  rule_stop <- NA_integer_
  F_path <- mse_path <- nr_path <- gauge_path <- perp_path <-
    numeric(max_sweeps)
  B_path <- if (track_B_path)
    array(NA_real_, c(nrow(B_prev), Km1, max_sweeps)) else NULL
  nB0 <- if (!is.null(Bz0)) sqrt(sum(Bz0^2)) else NA_real_

  s_used <- 0L
  for (s in seq_len(max_sweeps)) {
    if (criterion == "dev") {
      zs <- dv_z_step(Z, Phi, Wn, V, nu = nu, n_gn = n_gn)
      Z <- zs$Z; nu <- zs$nu; n_fail <- n_fail + zs$n_fail
      Phi <- dv_phi_em(Z, Phi, Wn, V, iters = em_iters)
    } else {
      zs <- bc_z_step(Z, Phi, Wf, V, lambda = 0, CB = NULL, nu = nu,
                      n_gn = n_gn)
      Z <- zs$Z; nu <- zs$nu; n_fail <- n_fail + zs$n_fail
      Phi <- fs_phi_step_proj(Z, Phi, Wf, V)$Phi
    }
    F_new <- fobj()
    if (F_new > F_cur + 1e-10 * (1 + abs(F_cur))) monotone_ok <- FALSE
    F_cur <- F_new
    F_path[s] <- F_new
    s_used <- s

    B_new <- bc_b_step(Z, C)
    if (track_B_path) B_path[, , s] <- B_new
    if (!is.null(Bz0)) {
      mse_path[s] <- procrustes_align(B_new, Bz0)$mse
      nr_path[s]  <- sqrt(sum(B_new^2)) / nB0
    }
    if (!is.null(QT)) {
      d <- bc_pack(Z, Phi) - eta0
      g2 <- sum(crossprod(QT, d)^2)
      gauge_path[s] <- sqrt(g2)
      perp_path[s]  <- sqrt(max(sum(d^2) - g2, 0))
    }
    relch <- max(abs(B_new - B_prev)) / (1e-8 + sqrt(mean(B_prev^2)))
    B_prev <- B_new
    hits <- if (relch < rule_tol) hits + 1L else 0L
    if (hits >= patience && is.na(rule_stop)) rule_stop <- s
    if (apply_rule && !is.na(rule_stop)) break
  }
  keep <- seq_len(s_used)
  list(Z = Z, Phi = Phi, B = B_prev, sweeps = s_used,
       rule_stop = rule_stop, monotone_ok = monotone_ok,
       n_fail = n_fail, F_path = F_path[keep],
       mse_path = if (!is.null(Bz0)) mse_path[keep] else NULL,
       nr_path = if (!is.null(Bz0)) nr_path[keep] else NULL,
       gauge_path = if (!is.null(QT)) gauge_path[keep] else NULL,
       perp_path = if (!is.null(QT)) perp_path[keep] else NULL,
       B_path = if (track_B_path) B_path[, , keep, drop = FALSE]
                else NULL)
}

# -------------------------------------------------------------------
#  Unit tests (hard gates; tiny instance M = 30, N = 40, K = 3)
# -------------------------------------------------------------------
dv_verify <- function(seed = 66001L) {
  dat <- sim_dgp(M = 30L, N = 40L, K = 3L, P = 2L, b_max = 0.5,
                 sigma_eps = 0.3, alpha_beta = 0.1, doc_length = 100L,
                 seed = seed)
  Wn <- dat$W; V <- dat$V
  set.seed(seed + 1L)
  Z   <- dat$Z_true + matrix(rnorm(60, 0, 0.1), 30, 2)
  Phi <- dat$Beta

  # (i) analytic z-gradient vs numDeriv (Phi fixed)
  fZ <- function(zv) dv_F(matrix(zv, 30, 2), Phi, Wn, V)
  Th <- bc_theta(Z, V)
  P  <- pmax(Th %*% Phi, 1e-12)
  U  <- (Wn / P) %*% t(Phi)
  TU <- Th * U
  G  <- -(TU - Th * rowSums(TU)) %*% V
  g_nd <- numDeriv::grad(fZ, as.vector(Z))
  rel <- sqrt(sum((as.vector(G) - g_nd)^2)) / sqrt(sum(g_nd^2))
  if (rel > 1e-6)
    stop(sprintf("dv z-gradient check FAILED: rel err %.2e", rel))

  # (ii) EM Phi-step decreases F_dev on 20 random starts
  set.seed(seed + 2L)
  em_ok <- 0L
  for (t in 1:20) {
    P0 <- matrix(rgamma(3 * 40, 0.5), 3, 40)
    P0 <- P0 / rowSums(P0)
    Zr <- matrix(rnorm(60, 0, 0.5), 30, 2)
    f_before <- dv_F(Zr, P0, Wn, V)
    f_after  <- dv_F(Zr, dv_phi_em(Zr, P0, Wn, V, iters = 5L), Wn, V)
    if (f_after <= f_before + 1e-8 * (1 + abs(f_before)))
      em_ok <- em_ok + 1L
  }
  if (em_ok < 20L)
    stop(sprintf("dv EM monotonicity FAILED: %d/20", em_ok))

  # (iii) full sweeps monotone from 5 random starts
  set.seed(seed + 3L)
  sweep_ok <- 0L
  for (t in 1:5) {
    P0 <- matrix(rgamma(3 * 40, 0.5), 3, 40); P0 <- P0 / rowSums(P0)
    Zr <- matrix(rnorm(60, 0, 0.5), 30, 2)
    r <- dv_refine(Zr, P0, Wn, dat$C, V, criterion = "dev",
                   max_sweeps = 10L, apply_rule = FALSE)
    if (r$monotone_ok) sweep_ok <- sweep_ok + 1L
  }
  if (sweep_ok < 5L)
    stop(sprintf("dv sweep monotonicity FAILED: %d/5", sweep_ok))

  list(grad_relerr = rel, em_ok = em_ok, sweep_ok = sweep_ok)
}
