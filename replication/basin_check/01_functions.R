#' ===================================================================
#' Basin-condition verification — core numerical routines
#' ===================================================================
#'
#' Implements the EXACT (non-linearised) refinement objective
#'
#'   F0(Z, Phi)      = || Wf - Theta(Z) %*% Phi ||_F^2
#'   F1(Z, Phi, B)   = F0 + lambda * || Z - C %*% B ||_F^2
#'
#' with Theta(Z) = softmax_rows(Z V'), Wf = row-normalised counts
#' (the same convention sgscatm() applies internally via scale_W = TRUE),
#' plus:
#'   - analytic gradients (verified against numDeriv on a tiny instance)
#'   - a damped block Gauss-Newton refinement (Z-step / Phi-step / B-step)
#'   - Hessian-vector products by central differences of the gradient
#'   - the gauge (general-linear reparametrisation) tangent basis
#'   - power iteration / Lanczos helpers for the spectral analysis
#'
#' For lambda > 0 the objective is always handled in PROFILED form:
#' B is the closed-form OLS minimiser given Z, so
#'   F1_prof(Z, Phi) = F0 + lambda * || (I - P_C) Z ||_F^2 ,
#' whose gradient in Z is grad(F0) + 2*lambda*(Z - C B(Z)) by the
#' envelope theorem.  The block descent's B-step realises the profiling.
#'
#' All functions are prefixed bc_ and touch no package internals.
#' ===================================================================


# -------------------------------------------------------------------
#  Objective and analytic gradient
# -------------------------------------------------------------------

#' Row-wise softmax of Z V'  (numerically stable)
bc_theta <- function(Z, V) {
  S <- tcrossprod(Z, V)                 # M x K
  S <- S - apply(S, 1L, max)
  E <- exp(S)
  E / rowSums(E)
}

#' Residual of the OLS projection of Z on C (profiled penalty term)
bc_z_resid <- function(Z, C, ridge = 1e-8) {
  B <- solve(crossprod(C) + ridge * diag(ncol(C)), crossprod(C, Z))
  Z - C %*% B
}

#' Profiled objective value
bc_objective <- function(Z, Phi, Wf, V, lambda = 0, C = NULL) {
  R <- bc_theta(Z, V) %*% Phi - Wf
  f <- sum(R * R)
  if (lambda > 0) {
    Zr <- bc_z_resid(Z, C)
    f  <- f + lambda * sum(Zr * Zr)
  }
  f
}

#' Profiled objective + analytic gradient
#'
#' Gradients (per document i, theta_i = softmax(V z_i)):
#'   r_i      = Phi' theta_i - w_i
#'   grad_zi  = 2 V' (diag(theta_i) - theta_i theta_i') Phi r_i
#'              (+ 2*lambda*(z_i - B(Z)' c_i) at lambda > 0, profiled)
#'   grad_Phi = 2 Theta' (Theta Phi - Wf)
#'
#' Vectorised over documents: with U = R Phi' (rows = (Phi r_i)'),
#'   Gz = 2 * (Theta*U - Theta * rowSums(Theta*U)) %*% V .
bc_grad <- function(Z, Phi, Wf, V, lambda = 0, C = NULL) {
  Theta <- bc_theta(Z, V)
  R  <- Theta %*% Phi - Wf              # M x N
  f  <- sum(R * R)
  U  <- R %*% t(Phi)                    # M x K
  TU <- Theta * U
  Gz   <- 2 * (TU - Theta * rowSums(TU)) %*% V
  Gphi <- 2 * crossprod(Theta, R)
  if (lambda > 0) {
    Zr <- bc_z_resid(Z, C)
    f  <- f + lambda * sum(Zr * Zr)
    Gz <- Gz + 2 * lambda * Zr
  }
  list(F = f, Gz = Gz, Gphi = Gphi)
}


# -------------------------------------------------------------------
#  eta = (vec(Z), vec(Phi)) packing (column-major)
# -------------------------------------------------------------------

bc_pack <- function(Z, Phi) c(as.vector(Z), as.vector(Phi))

bc_unpack <- function(eta, M, Km1, K, N) {
  nz <- M * Km1
  list(Z   = matrix(eta[seq_len(nz)], M, Km1),
       Phi = matrix(eta[nz + seq_len(K * N)], K, N))
}

#' Gradient as a flat vector of eta (for numDeriv checks and HVPs)
bc_grad_eta <- function(eta, dims, Wf, V, lambda = 0, C = NULL) {
  p <- bc_unpack(eta, dims$M, dims$Km1, dims$K, dims$N)
  g <- bc_grad(p$Z, p$Phi, Wf, V, lambda, C)
  bc_pack(g$Gz, g$Gphi)
}

bc_obj_eta <- function(eta, dims, Wf, V, lambda = 0, C = NULL) {
  p <- bc_unpack(eta, dims$M, dims$Km1, dims$K, dims$N)
  bc_objective(p$Z, p$Phi, Wf, V, lambda, C)
}

#' Hessian-vector product operator at a fixed point eta
#'
#' Central finite differences of the analytic gradient:
#'   Hv = (g(eta + h v) - g(eta - h v)) / (2h),
#' h = h_scale * (1 + ||eta||) / ||v||   (h_scale tuned on the tiny
#' instance against the full numDeriv Hessian; default 1e-5).
bc_hvp_factory <- function(eta, dims, Wf, V, lambda = 0, C = NULL,
                           h_scale = 1e-5) {
  n_eta <- sqrt(sum(eta^2))
  function(v) {
    nv <- sqrt(sum(v^2))
    if (nv == 0) return(v * 0)
    h <- h_scale * (1 + n_eta) / nv
    (bc_grad_eta(eta + h * v, dims, Wf, V, lambda, C) -
     bc_grad_eta(eta - h * v, dims, Wf, V, lambda, C)) / (2 * h)
  }
}


# -------------------------------------------------------------------
#  Block Gauss-Newton refinement
# -------------------------------------------------------------------

#' Phi-step: ridge-stabilised least squares given Z
bc_phi_step <- function(Z, Wf, V, ridge = 1e-8) {
  Th <- bc_theta(Z, V)
  solve(crossprod(Th) + ridge * diag(ncol(Th)), crossprod(Th, Wf))
}

#' B-step: closed-form OLS of Z on C
bc_b_step <- function(Z, C, ridge = 1e-8) {
  solve(crossprod(C) + ridge * diag(ncol(C)), crossprod(C, Z))
}

#' Z-step: n_gn damped Gauss-Newton iterations per document, vectorised
#'
#' Per document the GN system is
#'   (J_i'J_i + lambda I + nu_i I) delta = -(J_i'r_i + lambda(z_i - B'c_i))
#' with J_i = Phi'(diag(theta_i) - theta_i theta_i')V, followed by
#' Armijo backtracking (halving, max_bt halvings).  Levenberg damping
#' nu_i starts at nu0 and is multiplied by 10 when a document fails to
#' find a decrease (the document then keeps its current z for this
#' iteration); on success it relaxes back towards nu0.
#'
#' All document-level quantities are computed in vectorised form; only
#' the (K-1)x(K-1) solves loop over documents.
bc_z_step <- function(Z, Phi, Wf, V, lambda = 0, CB = NULL, nu,
                      n_gn = 2L, max_bt = 30L, ridge = 1e-8) {
  M   <- nrow(Z); Km1 <- ncol(Z)
  A   <- tcrossprod(Phi)                # K x K
  Pw  <- Wf %*% t(Phi)                  # M x K
  w2  <- rowSums(Wf * Wf)               # M

  # per-document objective f_i(z_i) (constant w2 included so f_i >= 0)
  fdoc <- function(Zm, rows) {
    Th <- bc_theta(Zm, V)
    f  <- rowSums((Th %*% A) * Th) - 2 * rowSums(Th * Pw[rows, , drop = FALSE]) +
          w2[rows]
    if (lambda > 0)
      f <- f + lambda * rowSums((Zm - CB[rows, , drop = FALSE])^2)
    f
  }

  n_fail_tot <- 0L
  for (it in seq_len(n_gn)) {
    Th  <- bc_theta(Z, V)
    ThA <- Th %*% A
    Uu  <- ThA - Pw                     # M x K, rows = (Phi r_i)'
    TU  <- Th * Uu
    G   <- 2 * (TU - Th * rowSums(TU)) %*% V
    f0  <- rowSums(ThA * Th) - 2 * rowSums(Th * Pw) + w2
    if (lambda > 0) {
      G  <- G + 2 * lambda * (Z - CB)
      f0 <- f0 + lambda * rowSums((Z - CB)^2)
    }

    # GN normal matrices: with B_i = (diag(theta_i)-theta_i theta_i')V,
    # J_i'J_i = B_i' A B_i.  Entries assembled by Km1^2 vectorised
    # contractions; B_i[k,a] = theta_ik * (V[k,a] - (theta_i'V)[a]).
    S  <- Th %*% V                      # M x Km1
    Ba <- vector("list", Km1)
    for (a in seq_len(Km1))
      Ba[[a]] <- Th * (matrix(V[, a], M, ncol(Th), byrow = TRUE) - S[, a])
    Ca <- lapply(Ba, function(B) B %*% A)     # A symmetric
    Harr <- array(0, c(M, Km1, Km1))
    for (a in seq_len(Km1)) for (b in seq_len(Km1))
      Harr[, a, b] <- rowSums(Ba[[a]] * Ca[[b]])

    # per-document solves (small (K-1)x(K-1) systems)
    Delta <- matrix(0, M, Km1)
    for (i in seq_len(M)) {
      Hi <- Harr[i, , ] + diag(lambda + nu[i] + ridge, Km1)
      Delta[i, ] <- tryCatch(-solve(Hi, G[i, ] / 2),
                             error = function(e) rep(0, Km1))
    }

    # vectorised Armijo backtracking (c1 = 1e-4)
    gd   <- rowSums(G * Delta)
    Delta[gd >= 0, ] <- 0               # not a descent direction: skip doc
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
      Znew[fail, ] <- Z[fail, , drop = FALSE]   # keep old z
      nu[fail]     <- nu[fail] * 10
      n_fail_tot   <- n_fail_tot + sum(fail)
    }
    nu[acc] <- pmax(nu[acc] / 3, 1e-6)
    Z <- Znew
  }
  list(Z = Z, nu = nu, n_fail = n_fail_tot)
}

#' Z-step with the B-block profiled out exactly (lambda > 0)
#'
#' The three-block descent (Z-step | Phi-step | B-step) has a slow
#' linear tail at lambda = 1: the penalty term dominates the
#' reconstruction term, so (Z, B) zig-zag along the valley Z ~ C B
#' (empirical rate ~0.989/sweep).  Profiling B out removes that mode:
#' the penalty becomes lambda ||(I - P_C) Z||^2 with Hessian
#' 2 lambda (I - Q Q') (x) I_{K-1}, Q = orth(C) — a per-document
#' diagonal plus a rank-P(K-1) coupling handled in closed form by the
#' Woodbury identity.  Writing A_i = J_i'J_i + (lambda + nu_i + ridge)I
#' and S = Q' Delta, the GN system
#'   [blockdiag(A_i) - lambda (QQ') (x) I] Delta = -g/2
#' reduces to per-document solves plus one P(K-1) x P(K-1) solve:
#'   delta_i = A_i^{-1}(-g_i/2) + lambda A_i^{-1} S' q_i,
#'   (I - lambda sum_i W_i (x) q_i q_i') vec(S) = vec(Q' X),
#' with W_i = A_i^{-1}, X the per-document solutions.  The step is
#' damped by a single global Armijo backtracking on the profiled F
#' (the direction couples documents, so per-document line searches no
#' longer apply); Levenberg damping is global on failure.
bc_z_step_prof <- function(Z, Phi, Wf, V, lambda, C, nu,
                           n_gn = 2L, max_bt = 30L, ridge = 1e-8) {
  M <- nrow(Z); Km1 <- ncol(Z)
  A  <- tcrossprod(Phi)
  Pw <- Wf %*% t(Phi)
  w2 <- rowSums(Wf * Wf)
  Qc <- qr.Q(qr(C))                   # M x P, P_C = Qc Qc'
  P  <- ncol(Qc)

  Fglob <- function(Zm) {
    Th <- bc_theta(Zm, V)
    Zr <- Zm - Qc %*% crossprod(Qc, Zm)
    sum(rowSums((Th %*% A) * Th) - 2 * rowSums(Th * Pw) + w2) +
      lambda * sum(Zr * Zr)
  }

  n_fail_tot <- 0L
  for (it in seq_len(n_gn)) {
    Th  <- bc_theta(Z, V)
    ThA <- Th %*% A
    TU  <- Th * (ThA - Pw)
    Zr  <- Z - Qc %*% crossprod(Qc, Z)
    G   <- 2 * (TU - Th * rowSums(TU)) %*% V + 2 * lambda * Zr
    f0  <- sum(rowSums(ThA * Th) - 2 * rowSums(Th * Pw) + w2) +
      lambda * sum(Zr * Zr)

    S  <- Th %*% V
    Ba <- vector("list", Km1)
    for (a in seq_len(Km1))
      Ba[[a]] <- Th * (matrix(V[, a], M, ncol(Th), byrow = TRUE) - S[, a])
    Ca <- lapply(Ba, function(B) B %*% A)
    Harr <- array(0, c(M, Km1, Km1))
    for (a in seq_len(Km1)) for (b in seq_len(Km1))
      Harr[, a, b] <- rowSums(Ba[[a]] * Ca[[b]])

    # per-document solves + Woodbury accumulation
    X    <- matrix(0, M, Km1)
    Wk   <- array(0, c(M, Km1, Km1))
    Kacc <- matrix(0, P * Km1, P * Km1)
    ok   <- TRUE
    for (i in seq_len(M)) {
      Ai <- Harr[i, , ] + diag(lambda + nu[i] + ridge, Km1)
      Wi <- tryCatch(solve(Ai), error = function(e) NULL)
      if (is.null(Wi)) { ok <- FALSE; break }
      Wk[i, , ] <- Wi
      X[i, ]    <- Wi %*% (-G[i, ] / 2)
      Kacc <- Kacc + kronecker(Wi, tcrossprod(Qc[i, ]))
    }
    if (!ok) { nu <- nu * 10; n_fail_tot <- n_fail_tot + 1L; next }
    Svec <- solve(diag(P * Km1) - lambda * Kacc,
                  as.vector(crossprod(Qc, X)))
    Smat <- matrix(Svec, P, Km1)
    Vv   <- Qc %*% Smat                # rows = (S' q_i)'
    Corr <- matrix(0, M, Km1)
    for (a in seq_len(Km1)) for (b in seq_len(Km1))
      Corr[, a] <- Corr[, a] + Wk[, a, b] * Vv[, b]
    Delta <- X + lambda * Corr

    # global Armijo backtracking (c1 = 1e-4), Levenberg on failure
    gd <- sum(G * Delta)
    if (gd >= 0) { nu <- nu * 10; n_fail_tot <- n_fail_tot + 1L; next }
    step <- 1; accepted <- FALSE
    for (bt in 0:max_bt) {
      fnew <- Fglob(Z + step * Delta)
      if (fnew <= f0 + 1e-4 * step * gd) { accepted <- TRUE; break }
      step <- step / 2
    }
    if (accepted) {
      Z  <- Z + step * Delta
      nu <- pmax(nu / 3, 1e-6)
    } else {
      nu <- nu * 10
      n_fail_tot <- n_fail_tot + 1L
    }
  }
  list(Z = Z, nu = nu, n_fail = n_fail_tot)
}

#' Full block Gauss-Newton refinement
#'
#' One sweep = Z-step (n_gn damped GN iterations per document),
#' Phi-step (ridge LS), B-step (OLS, lambda > 0 only).
#'
#' Convergence: relative F decrease over a sweep < tol_f  AND
#' max_i ||grad_zi||_inf < tol_g_scale * (1 + F), or max_sweeps.
#' Monotone decrease of the (profiled) objective is checked at every
#' sweep with tolerance 1e-12*(1+|F|); a violation sets
#' `monotone_violation = TRUE` and stops the run with diagnostics kept
#' in the trace (flag-and-record rather than stop(), so parallel jobs
#' retain their results — deviation logged in REPORT.md).
#'
#' @param Bz0 optional true coefficients: if given, mse_Bz (paper
#'   metric: OLS of current Z on C, then procrustes_align) is traced.
#' @param profile_B at lambda > 0: TRUE (default) profiles B into the
#'   Z-step (bc_z_step_prof, two-block descent); FALSE keeps the
#'   three-block variant with per-document decoupled Z-steps.
bc_refine <- function(Z0, Phi0, Wf, C, V, lambda = 0, Bz0 = NULL,
                      max_sweeps = 100L, n_gn = 2L,
                      tol_f = 1e-10, tol_g_scale = 1e-7,
                      profile_B = TRUE) {
  Z <- Z0; Phi <- Phi0
  M <- nrow(Z)
  B  <- if (lambda > 0) bc_b_step(Z, C) else NULL
  nu <- rep(1e-6, M)

  g0    <- bc_grad(Z, Phi, Wf, V, lambda, C)
  F_cur <- g0$F
  trace <- data.frame(sweep = 0L, F = F_cur,
                      gmax = max(abs(g0$Gz)),
                      mse_Bz = if (!is.null(Bz0))
                        procrustes_align(bc_b_step(Z, C), Bz0)$mse else NA_real_,
                      n_fail = 0L)
  monotone_violation <- FALSE
  converged <- FALSE

  for (s in seq_len(max_sweeps)) {
    zs <- if (lambda > 0 && profile_B) {
      bc_z_step_prof(Z, Phi, Wf, V, lambda, C, nu = nu, n_gn = n_gn)
    } else {
      bc_z_step(Z, Phi, Wf, V, lambda,
                CB = if (lambda > 0) C %*% B else NULL,
                nu = nu, n_gn = n_gn)
    }
    Z   <- zs$Z; nu <- zs$nu
    Phi <- bc_phi_step(Z, Wf, V)
    if (lambda > 0) B <- bc_b_step(Z, C)

    g     <- bc_grad(Z, Phi, Wf, V, lambda, C)
    F_new <- g$F
    gmax  <- max(abs(g$Gz))
    trace <- rbind(trace, data.frame(
      sweep = s, F = F_new, gmax = gmax,
      mse_Bz = if (!is.null(Bz0))
        procrustes_align(bc_b_step(Z, C), Bz0)$mse else NA_real_,
      n_fail = zs$n_fail))

    if (F_new > F_cur + 1e-12 * (1 + abs(F_cur))) {
      monotone_violation <- TRUE
      warning(sprintf(
        "bc_refine: monotonicity violated at sweep %d (F %.15e -> %.15e)",
        s, F_cur, F_new))
      break
    }
    rel_dec <- (F_cur - F_new) / (1 + abs(F_cur))
    F_cur   <- F_new
    if (rel_dec < tol_f && gmax < tol_g_scale * (1 + abs(F_cur))) {
      converged <- TRUE
      break
    }
  }

  list(Z = Z, Phi = Phi, B = if (lambda > 0) B else bc_b_step(Z, C),
       F = F_cur, sweeps = max(trace$sweep), converged = converged,
       monotone_violation = monotone_violation, trace = trace)
}


# -------------------------------------------------------------------
#  Gauge subspace (lambda = 0 reparametrisation invariance)
# -------------------------------------------------------------------

#' Orthonormal basis of the first-order gauge directions at (Z, Phi)
#'
#' F0 is invariant under (Theta, Phi) -> (Theta G^{-1}, G Phi) with
#' G = I + tE, E 1_K = 0.  First-order directions:
#'   dPhi     = E Phi
#'   dtheta_i = -E' theta_i
#'   dz_i     = V' (dtheta_i / theta_i)      (z = V' log theta)
#' Basis for {E : E 1 = 0}: E_ab = e_a (e_b - e_K)', a in 1..K,
#' b in 1..(K-1) — dimension K*(K-1).
#'
#' Returns the d x K(K-1) orthonormalised matrix Q_T (d = M(K-1)+KN).
bc_gauge_basis <- function(Z, Phi, V) {
  Theta <- bc_theta(Z, V)
  K <- ncol(Theta); Km1 <- K - 1L
  M <- nrow(Z); N <- ncol(Phi)
  Tmat <- matrix(0, M * Km1 + K * N, K * Km1)
  col <- 0L
  for (a in seq_len(K)) for (b in seq_len(Km1)) {
    col <- col + 1L
    E <- matrix(0, K, K)
    E[a, b] <- 1; E[a, K] <- E[a, K] - 1
    dPhi   <- E %*% Phi
    dTheta <- -Theta %*% E              # row i = (-E' theta_i)'
    dZ     <- (dTheta / Theta) %*% V
    Tmat[, col] <- c(as.vector(dZ), as.vector(dPhi))
  }
  qr.Q(qr(Tmat))
}


# -------------------------------------------------------------------
#  Spectral helpers
# -------------------------------------------------------------------

#' Power iteration for a symmetric operator; returns the dominant
#' Rayleigh quotient (signed) and its absolute value (operator norm).
bc_power_iter <- function(Afun, d, iters = 50L, seed = 1L) {
  set.seed(seed)
  v <- rnorm(d); v <- v / sqrt(sum(v^2))
  lam <- 0
  for (i in seq_len(iters)) {
    w  <- Afun(v)
    lam <- sum(v * w)
    nw <- sqrt(sum(w^2))
    if (nw == 0) break
    v <- w / nw
  }
  list(value = lam, norm = abs(lam), vector = v)
}

#' Smallest-k eigenpairs of a symmetric operator via Lanczos on the
#' shifted operator  B = shift*I - H  (top of B = bottom of H).
#' With deflation (QT given), B = shift*P - P H P, P = I - QT QT':
#' gauge directions map to 0, the orthocomplement spectrum maps to
#' shift - lambda, so the top-k of B are the smallest-k of P H P
#' restricted to span(QT)^perp.
bc_smallest_eigs <- function(hvp, d, k = 30L, shift, QT = NULL,
                             ncv = NULL, tol = 1e-7, maxitr = 1000L) {
  if (is.null(ncv)) ncv <- min(d, max(4L * k, 80L))
  proj <- if (is.null(QT)) identity else
    function(x) x - QT %*% crossprod(QT, x)
  Bop <- function(x, args) {
    if (is.null(QT)) shift * x - hvp(x)
    else { xp <- proj(x); as.vector(proj(shift * xp - hvp(xp))) }
  }
  res <- RSpectra::eigs_sym(Bop, k = k, which = "LM", n = d,
                            opts = list(ncv = ncv, tol = tol,
                                        maxitr = maxitr, retvec = TRUE))
  ord <- order(res$values, decreasing = TRUE)     # largest B = smallest H
  list(values  = shift - res$values[ord],
       vectors = res$vectors[, ord, drop = FALSE],
       nconv   = res$nconv, niter = res$niter)
}

#' Principal cosines between span(Q1) and span(Q2) (orthonormal inputs)
bc_principal_cosines <- function(Q1, Q2) {
  sv <- svd(crossprod(Q1, Q2))
  pmin(sv$d, 1)
}


# -------------------------------------------------------------------
#  Pilot alignment and Bz metrics
# -------------------------------------------------------------------

#' Oracle general-linear (GL) alignment of pilot scores to Z_true:
#'   A_hat = (Zhat'Zhat)^{-1} Zhat' Zc,   Zc = centred Z_true,
#'   Z_al  = Zhat A_hat + 1 colMeans(Z_true)'.
bc_gl_align <- function(Zhat, Z_true, ridge = 1e-8) {
  Zc <- scale(Z_true, center = TRUE, scale = FALSE)
  A  <- solve(crossprod(Zhat) + ridge * diag(ncol(Zhat)),
              crossprod(Zhat, Zc))
  Zal <- Zhat %*% A +
    matrix(colMeans(Z_true), nrow(Zhat), ncol(Zhat), byrow = TRUE)
  list(Z = Zal, A = A)
}

#' Orthogonal-Procrustes-only alignment (rotation, no scaling), same
#' centring convention as bc_gl_align for comparability.
bc_op_align <- function(Zhat, Z_true) {
  Zc <- scale(Z_true, center = TRUE, scale = FALSE)
  sv <- svd(crossprod(Zhat, Zc))
  R  <- sv$u %*% t(sv$v)
  Zal <- Zhat %*% R +
    matrix(colMeans(Z_true), nrow(Zhat), ncol(Zhat), byrow = TRUE)
  list(Z = Zal, R = R)
}

#' Paper metric: OLS of Z on C, then orthogonal Procrustes to Bz0
bc_mse_paper <- function(Z, C, Bz0) {
  procrustes_align(bc_b_step(Z, C), Bz0)$mse
}

#' Direct metric: OLS of Z on C, entry-wise MSE without any rotation
#' (used with GL-aligned or model-coordinate Z, where the alignment
#' already absorbed the rotation)
bc_mse_direct <- function(Z, C, Bz0) {
  mean((bc_b_step(Z, C) - Bz0)^2)
}

#' rho components between (Z, Phi) and the truth
bc_rho <- function(Z, Phi, Z_true, Phi_true) {
  dz <- sqrt(sum((Z - Z_true)^2))
  dp <- sqrt(sum((Phi - Phi_true)^2))
  c(rho = sqrt(dz^2 + dp^2), rho_Z = dz, rho_Phi = dp)
}


# -------------------------------------------------------------------
#  Derivative verification on a tiny instance (hard-stop gate)
# -------------------------------------------------------------------

#' Verify analytic gradient, HVP step size, and gauge tangency on a
#' tiny instance (M = 30, N = 40, K = 3).  Stops on failure.
bc_verify_derivatives <- function(seed = 4242L, h_scale = 1e-5,
                                  tol_grad = 1e-6) {
  dat <- sim_dgp(M = 30L, N = 40L, K = 3L, P = 2L, b_max = 0.5,
                 sigma_eps = 0.3, alpha_beta = 0.1, doc_length = 100L,
                 seed = seed)
  Wf <- dat$W / rowSums(dat$W)
  V  <- dat$V
  dims <- list(M = 30L, Km1 = 2L, K = 3L, N = 40L)
  set.seed(seed + 1L)
  Z   <- dat$Z_true + matrix(rnorm(30 * 2, 0, 0.1), 30, 2)
  Phi <- dat$Beta   + matrix(rnorm(3 * 40, 0, 0.01), 3, 40)
  eta <- bc_pack(Z, Phi)

  out <- list()
  for (lam in c(0, 1)) {
    g_an <- bc_grad_eta(eta, dims, Wf, V, lam, dat$C)
    g_nd <- numDeriv::grad(bc_obj_eta, eta, dims = dims, Wf = Wf, V = V,
                           lambda = lam, C = dat$C)
    rel <- sqrt(sum((g_an - g_nd)^2)) / sqrt(sum(g_nd^2))
    out[[sprintf("grad_relerr_lambda%g", lam)]] <- rel
    if (rel > tol_grad)
      stop(sprintf("Gradient check FAILED (lambda=%g): rel err %.3e > %.0e",
                   lam, rel, tol_grad))

    # HVP vs numDeriv Jacobian of the analytic gradient
    Jac <- numDeriv::jacobian(bc_grad_eta, eta, dims = dims, Wf = Wf,
                              V = V, lambda = lam, C = dat$C)
    Hnd <- (Jac + t(Jac)) / 2
    hvp <- bc_hvp_factory(eta, dims, Wf, V, lam, dat$C, h_scale)
    set.seed(seed + 2L)
    relh <- vapply(1:5, function(j) {
      v <- rnorm(length(eta))
      sqrt(sum((hvp(v) - Hnd %*% v)^2)) / sqrt(sum((Hnd %*% v)^2))
    }, numeric(1))
    out[[sprintf("hvp_relerr_lambda%g", lam)]] <- max(relh)
    if (max(relh) > 1e-4)
      stop(sprintf("HVP check FAILED (lambda=%g): rel err %.3e > 1e-4",
                   lam, max(relh)))
  }

  # gauge tangency: g . t = 0 exactly for every gauge direction (lambda=0)
  QT <- bc_gauge_basis(Z, Phi, V)
  g0 <- bc_grad_eta(eta, dims, Wf, V, 0, NULL)
  gt <- max(abs(crossprod(QT, g0)))
  out$gauge_tangency <- gt / (1 + sqrt(sum(g0^2)))
  if (out$gauge_tangency > 1e-8)
    stop(sprintf("Gauge tangency check FAILED: %.3e", out$gauge_tangency))

  out
}
