#' Asymptotic Covariance of the ILR Path Coefficients (Standardized Scale)
#'
#' Computes the asymptotic variance-covariance matrix and standard errors of
#' \eqn{\mathrm{vec}(\hat{\mathbf{B}}_z)} on the **standardized** score scale,
#' via the empirical influence function of the spectral estimator.
#'
#' @details
#' The estimator standardizes the ILR scores as
#' \eqn{\tilde{\mathbf{Z}} = \sqrt{M}\,\mathbf{Z}^\ast} so that
#' \eqn{M^{-1}\tilde{\mathbf{Z}}^\top\tilde{\mathbf{Z}} = \mathbf{I}}, and
#' \eqn{\hat{\mathbf{B}}_z = (\mathbf{C}^\top\mathbf{C})^{-1}\mathbf{C}^\top
#' \tilde{\mathbf{Z}}}. Writing the standardized score of document \eqn{i} on
#' component \eqn{k} as \eqn{\tilde z_{ik} = \rho_k^{-1/2}\,\mathbf{v}_k^\top
#' \tilde{\mathbf{w}}_i}, where \eqn{(\mathbf{v}_k,\rho_k)} are the eigenvectors
#' and eigenvalues **on the O(1) scale** of the word Gram
#' \eqn{\boldsymbol{\Sigma}_W = M^{-1}\tilde{\mathbf{W}}^\top\tilde{\mathbf{W}}}
#' (i.e. \eqn{\rho_k} equals the eigenvalue of
#' \eqn{M^{-1}\tilde{\mathbf W}\tilde{\mathbf W}^\top}, never the O(M) raw
#' eigenvalue), the per-document influence of column \eqn{k} is the sum of
#' three terms:
#' \deqn{\psi_{ik} = \underbrace{\boldsymbol{\Sigma}_C^{-1}\mathbf{c}_i
#'   (\tilde z_{ik} - \mathbf{c}_i^\top \hat{\mathbf b}_k)}_{\text{(A) regression}}
#'   + \underbrace{\sum_{l\neq k} (\mathbf{v}_l^\top\tilde{\mathbf w}_i)
#'   (\mathbf{v}_k^\top\tilde{\mathbf w}_i)\,\mathbf{g}_{kl}}_{\text{(B) eigenvector}}
#'   - \underbrace{\tfrac12\,\hat{\mathbf b}_k (\tilde z_{ik}^2 - 1)}_{\text{(C) eigenvalue}},}
#' with \eqn{\mathbf{g}_{kl} = \rho_k^{-1/2}\boldsymbol{\Sigma}_C^{-1}
#' \boldsymbol{\Sigma}_{CW}\mathbf{v}_l / (\rho_k - \rho_l)}. Term (B) is the
#' fluctuation of the eigenvector and term (C) — **new** relative to the
#' previous implementation — is the fluctuation of the eigenvalue used to
#' standardize the score. The covariance is a *mean* over documents and is
#' divided by \eqn{M} at the end:
#' \deqn{\hat{\mathbf V} = M^{-1}\sum_i \psi_i\psi_i^\top, \qquad
#'   \widehat{\mathrm{SE}}(\mathrm{vec}\,\hat{\mathbf B}_z)
#'   = \sqrt{\mathrm{diag}(\hat{\mathbf V})/M}.}
#' The final \eqn{1/M} (missing before) makes the SE \eqn{O(M^{-1/2})} rather
#' than \eqn{O(1)}.
#'
#' The regularisation weight \eqn{\lambda} does **not** enter at first order:
#' with \eqn{\lambda} fixed the covariate block \eqn{\sqrt{\lambda}\mathbf{Q}_C}
#' perturbs the leading word eigenspace by \eqn{O(\lambda/M)}, which is
#' asymptotically negligible; the influence above therefore depends on the data
#' only through the word Gram and the covariate cross-moment, not on
#' \eqn{\lambda}.
#'
#' **Rotational identification.** The spectral estimand is identified only up
#' to an orthogonal rotation of the score basis. When the estimate is aligned
#' to a reference (Procrustes in simulation, varimax in practice), the aligned
#' estimator has no variance along the rotational tangent
#' \eqn{\{\mathrm{vec}(\hat{\mathbf B}_z\boldsymbol{\Omega}):
#' \boldsymbol{\Omega}^\top=-\boldsymbol{\Omega}\}}. With `identified = TRUE`
#' (default) that tangent is projected out of \eqn{\hat{\mathbf V}}, which is
#' required for the SEs to match the sampling variability of the aligned
#' estimator (empirical SD / jackknife).
#'
#' @param fit An `"sgscatm"` object from [sgscatm()]. The centred DTM
#'   `fit$W_tilde` and centred covariates `fit$C_centred` are used.
#' @param W,C Optional; if supplied their dimensions are checked against `fit`.
#'   The computation uses the fit's stored centred matrices.
#' @param rotation Optional \eqn{(K-1)\times(K-1)} orthogonal matrix
#'   \eqn{\mathbf R} (e.g. Procrustes or varimax). If supplied, SEs are returned
#'   for the rotated estimate \eqn{\hat{\mathbf B}_z\mathbf R} via
#'   \eqn{(\mathbf R^\top\!\otimes\mathbf I_P)\hat{\mathbf V}
#'   (\mathbf R\otimes\mathbf I_P)}.
#' @param identified Logical; if TRUE (default) project out the rotational
#'   tangent (see Details).
#' @param r Integer; retained rank for the eigenvector sum (default
#'   `min(100, N, M-1)`). Bulk directions contribute negligibly.
#' @param tol Numeric tolerance for dropping (near-)degenerate eigenpairs from
#'   the term-(B) sum. Default 1e-8 on the O(1) scale.
#'
#' @return A list with `vcov` (\eqn{\hat{\mathbf V}/M}, on the standardized
#'   scale), `se` (a \eqn{P\times(K-1)} matrix), `B` (the standardized, possibly
#'   rotated, \eqn{\hat{\mathbf B}_z}), `rho` (the O(1) eigenvalues
#'   \eqn{\rho_1,\dots,\rho_{K-1}}), and `scale = "standardized"`.
#' @seealso [sgscatm()], [ilr_se()]
#' @export
sgscatm_vcov <- function(fit, W = NULL, C = NULL, rotation = NULL,
                         identified = TRUE, r = NULL, tol = 1e-8) {
  stopifnot(inherits(fit, "sgscatm"))
  if (is.null(fit$W_tilde) || is.null(fit$C_centred))
    stop("fit lacks stored W_tilde/C_centred; refit with current sgscatm().")

  Wt <- fit$W_tilde                        # M x N centred (scaled) DTM
  Cc <- fit$C_centred                      # M x P centred covariates
  M  <- nrow(Wt); N <- ncol(Wt); P <- ncol(Cc)
  K  <- fit$K; Km1 <- K - 1L

  if (!is.null(W) && any(dim(as.matrix(W)) != c(M, N)))
    warning("dim(W) does not match fit; using fit$W_tilde.")
  if (!is.null(C) && nrow(as.matrix(C)) != M)
    warning("nrow(C) does not match fit; using fit$C_centred.")

  # --- word Gram eigen on the O(1) scale: rho = eig( (1/M) W~' W~ ) --------
  SigmaW <- crossprod(Wt) / M                       # N x N (never M x M)
  eg  <- eigen(SigmaW, symmetric = TRUE)
  rho <- eg$values; Vv <- eg$vectors
  if (is.null(r)) r <- min(100L, N, M - 1L)
  r   <- min(r, sum(rho > tol * max(rho)), N)
  r   <- max(r, Km1)
  rho <- rho[seq_len(r)]; Vv <- Vv[, seq_len(r), drop = FALSE]

  # --- covariate moments ---------------------------------------------------
  SigC    <- crossprod(Cc) / M                      # P x P
  SigCinv <- .safe_chol_solve(SigC)
  SigCW   <- crossprod(Cc, Wt) / M                  # P x N

  A    <- Wt %*% Vv                                  # M x r : a_il = v_l' w~_i
  SWV  <- SigCinv %*% (SigCW %*% Vv)                 # P x r
  rk   <- sqrt(rho[seq_len(Km1)])
  ztil <- sweep(A[, seq_len(Km1), drop = FALSE], 2L, rk, "/")   # M x (K-1)
  Bhat <- sweep(SWV[, seq_len(Km1), drop = FALSE], 2L, rk, "/") # P x (K-1)
  CSinv <- Cc %*% SigCinv                            # M x P

  psi <- array(0, dim = c(M, P, Km1))
  for (k in seq_len(Km1)) {
    # (A) regression term
    resid <- as.numeric(ztil[, k] - Cc %*% Bhat[, k])
    tA <- CSinv * resid
    # (B) eigenvector-fluctuation term
    dk <- 1 / (rho[k] - rho); dk[k] <- 0
    dk[!is.finite(dk)] <- 0
    dk[abs(rho[k] - rho) < tol * max(rho)] <- 0      # guard degeneracy
    Gk <- sweep(SWV, 2L, dk / sqrt(rho[k]), "*")     # P x r  (col l = g_kl)
    tB <- A[, k] * (A %*% t(Gk))
    # (C) eigenvalue-fluctuation term
    tC <- outer(ztil[, k]^2 - 1, -0.5 * Bhat[, k])
    psi[, , k] <- tA + tB + tC
  }

  Vhat <- crossprod(matrix(psi, M, P * Km1)) / M     # O(1), P(K-1) square
  B    <- Bhat

  if (!is.null(rotation)) {
    R <- as.matrix(rotation)
    if (any(dim(R) != c(Km1, Km1))) stop("rotation must be (K-1)x(K-1).")
    Kk   <- kronecker(t(R), diag(P))
    Vhat <- Kk %*% Vhat %*% t(Kk)
    B    <- Bhat %*% R
  }
  if (isTRUE(identified)) {
    Pr   <- .rot_tangent_projector(B)
    Vhat <- Pr %*% Vhat %*% Pr
  }

  se <- matrix(sqrt(pmax(diag(Vhat), 0) / M), P, Km1)
  list(vcov = Vhat / M, se = se, B = B, rho = rho[seq_len(Km1)],
       scale = "standardized")
}

# Orthogonal projector removing the rotational tangent of a P x (K-1) matrix B,
# i.e. span{ vec(B (E_ij - E_ji)) : i<j }.
.rot_tangent_projector <- function(B) {
  P <- nrow(B); Km1 <- ncol(B)
  if (Km1 < 2L) return(diag(P * Km1))
  cols <- vector("list", Km1 * (Km1 - 1L) / 2L); m <- 0L
  for (i in seq_len(Km1 - 1L)) for (j in (i + 1L):Km1) {
    O <- matrix(0, Km1, Km1); O[i, j] <- 1; O[j, i] <- -1
    m <- m + 1L; cols[[m]] <- as.vector(B %*% O)
  }
  D <- do.call(cbind, cols)
  qrD <- qr(D)
  Q <- qr.Q(qrD)[, seq_len(qrD$rank), drop = FALSE]
  diag(P * Km1) - tcrossprod(Q)
}

# Guarded solve for small SPD matrices (Cholesky, QR fallback).
.safe_chol_solve <- function(A, tol = 1e-12) {
  ch <- tryCatch(chol(A), error = function(e) NULL)
  if (!is.null(ch)) return(chol2inv(ch))
  qr.solve(qr(A, tol = tol), diag(nrow(A)))
}
