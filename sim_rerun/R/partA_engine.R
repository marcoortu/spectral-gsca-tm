# ===================================================================
#  partA_engine.R  —  corrected estimator, dual-form analytic SE,
#                     and leave-one-out jackknife (Part A)
# ===================================================================
#  Reuses the package spectral solver sgscatm() and closure map; adds
#  the corrected standardized scale Z_tilde = sqrt(M) Z*, the dual
#  fixed-dimension influence-function SE, and an independent LOO
#  jackknife cross-check.
# ===================================================================

# --- Build the identity-normalized true coefficient matrix B0 -------
# P x (K-1) with orthonormal ROWS (R R' = I_P) scaled by sqrt(rho2):
#   B0 B0' = rho2 I_P  (exactly).  Isotropic residual sigma_eps=sqrt(1-rho2).
# With P < K-1 an exactly-identity (K-1)x(K-1) score covariance is
# unattainable (signal is rank P); diag(Cov(z)) is reported honestly.
build_B0 <- function(K = 5L, P = 3L, rho2 = 0.5, seed = 20260703L) {
  Km1 <- K - 1L
  stopifnot(P <= Km1)                      # rows orthonormal needs P <= K-1
  set.seed(seed)
  Rraw <- matrix(rnorm(Km1 * P), Km1, P)   # (K-1) x P
  Q    <- qr.Q(qr(Rraw))                   # (K-1) x P, orthonormal columns
  Rmat <- t(Q)                             # P x (K-1), R R' = I_P
  B0   <- Rmat * sqrt(rho2)                # P x (K-1)
  err  <- max(abs(tcrossprod(B0) - rho2 * diag(P)))
  stopifnot(err < 1e-10)                   # ||B0 B0' - rho2 I_P|| < 1e-10
  attr(B0, "rho2")      <- rho2
  attr(B0, "sigma_eps") <- sqrt(1 - rho2)
  attr(B0, "rowortho_err") <- err
  B0
}

# --- Corrected estimator on one replicate --------------------------
# Returns the point estimate on the corrected scale plus everything the
# SE routine needs.  b_hat = sqrt(M) * fit$Bz  (reuses package solver).
corrected_fit <- function(dat, B0, lambda = 1) {
  K   <- dat$params$K
  M   <- nrow(dat$W)
  t0  <- proc.time()
  fit <- sgscatm(dat$W, dat$C, K = K, lambda = lambda, rotate = TRUE)
  t_fit <- (proc.time() - t0)[3L]

  Zstar   <- fit$Z                       # M x (K-1), Z*'Z* = I
  Ztilde  <- sqrt(M) * Zstar             # standardized scores
  b_hat   <- sqrt(M) * fit$Bz            # = (C'C)^{-1} C' Ztilde
  pa      <- procrustes_align(b_hat, B0) # align to truth (SVD of b_hat'B0)

  list(fit = fit, M = M, Ztilde = Ztilde, b_hat = b_hat,
       R_hat = pa$R, b_hat_al = pa$Bz_aligned, time_fit = t_fit)
}

# --- Dual fixed-dimension influence-function SE --------------------
# H = [W_tilde  sqrt(lambda) Q_C]  (M x (N+P))
# G = H'H ;  eigen(G/M) -> (V_dual, rho).  Empirical influence per doc,
# vectorized over documents.  SE(vec b_hat_al) = sqrt(diag(Sigma/M)).
dual_influence_se <- function(cf, dat, B0, lambda = 1) {
  fit <- cf$fit
  M   <- cf$M
  Wt  <- fit$W_tilde                     # M x N centred DTM (package-stored)
  Cc  <- fit$C_centred                   # M x P centred covariates
  N   <- ncol(Wt);  P <- ncol(Cc);  Km1 <- fit$K - 1L
  NP  <- N + P

  Q_C <- qr.Q(qr(Cc))                    # M x P
  H   <- cbind(Wt, sqrt(lambda) * Q_C)   # M x (N+P)   (never M x M densely)
  G   <- crossprod(H)                    # (N+P) x (N+P)
  eg  <- eigen(G / M, symmetric = TRUE)
  rho <- eg$values                       # length N+P (decreasing)
  Vd  <- eg$vectors                      # (N+P) x (N+P)
  rho_top <- rho[seq_len(Km1)]

  # Raw dual spectral scores (unit norm): u_k = H v_k / sqrt(M rho_k)
  pos   <- pmax(rho_top, .Machine$double.eps)
  Zraw  <- (H %*% Vd[, seq_len(Km1), drop = FALSE]) %*%
           diag(1 / sqrt(M * pos), Km1)  # M x (K-1)

  # Rotation mapping raw dual basis -> package rotated basis (varimax)
  svT   <- svd(crossprod(Zraw, fit$Z))
  Tmat  <- svT$u %*% t(svT$v)            # (K-1) x (K-1) orthogonal
  R_tot <- Tmat %*% cf$R_hat             # raw basis -> aligned-to-B0

  # Raw-basis point estimate + residuals (term 1 is OLS influence)
  Ct   <- crossprod(Cc)
  B_raw    <- sqrt(M) * safe_solve(Ct, crossprod(Cc, Zraw))   # P x (K-1)
  Ztil_raw <- sqrt(M) * Zraw                                  # M x (K-1)
  Resid    <- Ztil_raw - Cc %*% B_raw                          # M x (K-1)

  SigC     <- Ct / M
  SigC_inv <- safe_solve(SigC)                                 # P x P

  # Precomputations (all O(M N (N+P)) at worst, never M x M dense)
  VW   <- Vd[seq_len(N), , drop = FALSE]      # N x (N+P) : W-block of eigvecs
  A    <- Wt %*% VW                            # M x (N+P): a_il = w~_i' vW_l
  Smat <- crossprod(A) / M                     # (N+P)x(N+P): vW_l'SigW vW_k
  Gcw  <- crossprod(Cc, A) / M                 # P x (N+P): (E[c w']) vW_l

  psi_raw <- array(0, dim = c(M, P, Km1))
  for (k in seq_len(Km1)) {
    rk <- rho[k]
    d  <- 1 / (rk - rho)                       # length N+P
    d[k] <- 0                                  # l != k
    if (P > 0L) d[(N + 1L):NP] <- 0            # l <= N (drop covariate dirs)
    d[!is.finite(d)] <- 0
    d[abs(rk - rho) < 1e-10] <- 0              # guard degenerate eigengap

    Ad   <- sweep(A, 2L, d, "*")               # M x (N+P)
    P1   <- Ad %*% t(Gcw)                       # M x P
    T1   <- A[, k] * P1                          # row-scale by a_ik
    cvec <- Gcw %*% (Smat[, k] * d)             # P-vector
    T2   <- matrix(cvec, M, P, byrow = TRUE)
    term2 <- (T1 - T2) / sqrt(pos[k])          # M x P

    term1 <- Cc * Resid[, k]                    # row-scale Cc by residual
    psi_raw[, , k] <- (term1 + term2) %*% SigC_inv
  }

  # Rotate influence to the aligned basis: psi_al_i = psi_raw_i %*% R_tot
  psi_al <- array(matrix(psi_raw, M * P, Km1) %*% R_tot, dim = c(M, P, Km1))

  # SE(vec b_hat_al) = sqrt( (1/M^2) sum_i psi_al_i^2 )
  se_mat <- sqrt(apply(psi_al^2, c(2L, 3L), sum) / M^2)

  # eigengap diagnostic (signal vs bulk)
  gap <- if (length(rho) > Km1)
    (rho_top[Km1] - rho[Km1 + 1L]) / max(abs(rho_top)) else NA_real_

  list(se = se_mat, rho_top = rho_top, eig_gap = gap)
}

# --- Independent leave-one-out jackknife SE (documents) -------------
# Expensive (M refits): call selectively.  Each LOO estimate is aligned
# to the full-sample aligned point estimate to remove rotation drift.
loo_jackknife_se <- function(dat, B_target_al, lambda = 1) {
  K <- dat$params$K; M <- nrow(dat$W); P <- ncol(dat$C); Km1 <- K - 1L
  acc  <- array(NA_real_, dim = c(M, P, Km1))
  for (i in seq_len(M)) {
    Wi <- dat$W[-i, , drop = FALSE]; Ci <- dat$C[-i, , drop = FALSE]
    fi <- tryCatch(sgscatm(Wi, Ci, K = K, lambda = lambda, rotate = TRUE),
                   error = function(e) NULL)
    if (is.null(fi)) next
    b_i <- sqrt(M - 1L) * fi$Bz
    b_i <- procrustes_align(b_i, B_target_al)$Bz_aligned
    acc[i, , ] <- b_i
  }
  ok    <- apply(acc, 1L, function(x) all(is.finite(x)))
  A_ok  <- acc[ok, , , drop = FALSE]
  n     <- sum(ok)
  mbar  <- apply(A_ok, c(2L, 3L), mean)
  ss    <- apply(sweep(A_ok, c(2L, 3L), mbar, "-")^2, c(2L, 3L), sum)
  se    <- sqrt((n - 1) / n * ss)
  list(se = se, n = n)
}
