#' ===================================================================
#'  Feasibility round — F2: split-document jackknife bias correction
#'  + sandwich covariance (verified in the audit)
#' ===================================================================
#'
#'  Per document, counts are split by cell-wise binomial thinning
#'  (a_ij ~ Bin(w_ij, 1/2), b = w - a).  Each half's z_hat is refit by
#'  damped GN with Phi held FIXED at the full-data refined value; the
#'  halves' OLS coefficients enter
#'      B_jk = 2 B_full - (B_A + B_B) / 2 ,
#'  which cancels the O(1/L) incidental-parameter bias to first order
#'  (halves have doubled bias 2c/L).
#' ===================================================================

#' Sandwich covariance for multivariate OLS (copied verbatim from
#' replication/audit_block1_stm/01_audit_block1.R::aud_sandwich — that
#' file is a run-script whose sourcing would launch the audit, so the
#' 15 lines are duplicated here with this note; unit-tested there
#' against direct simulation, rel. Frobenius error 0.040).
fs_sandwich <- function(Z, C, B_hat) {
  M <- nrow(C); P <- ncol(C); Kp <- ncol(Z)
  R    <- Z - C %*% B_hat
  XtXi <- solve(crossprod(C))
  Meat <- matrix(0, P * Kp, P * Kp)
  for (k in seq_len(Kp)) for (kp in seq_len(Kp)) {
    blk <- crossprod(C, C * (R[, k] * R[, kp]))
    Meat[(k - 1L) * P + seq_len(P), (kp - 1L) * P + seq_len(P)] <- blk
  }
  Bread <- kronecker(diag(Kp), XtXi)
  (M / (M - P)) * (Bread %*% Meat %*% Bread)
}

#' Cell-wise binomial thinning split of a count matrix
fs_thin_split <- function(W, seed) {
  set.seed(seed)
  A <- matrix(rbinom(length(W), as.vector(W), 0.5), nrow(W), ncol(W))
  list(A = A, B = W - A)
}

#' Refit per-document z on a half's frequencies, Phi fixed, warm start
fs_zfit_half <- function(Wpart, Phi, V, Z_start, n_gn = 8L) {
  rs <- rowSums(Wpart)
  Wf <- Wpart / pmax(rs, 1)
  bc_z_step(Z_start, Phi, Wf, V, lambda = 0, CB = NULL,
            nu = rep(1e-6, nrow(Wpart)), n_gn = n_gn)$Z
}

#' Split-document jackknife around a full fit (Z_full, Phi_full)
fs_jackknife_B <- function(W, C, V, Z_full, Phi_full, seed) {
  th <- fs_thin_split(W, seed)
  Z_A <- fs_zfit_half(th$A, Phi_full, V, Z_full)
  Z_B <- fs_zfit_half(th$B, Phi_full, V, Z_full)
  B_full <- bc_b_step(Z_full, C)
  B_A <- bc_b_step(Z_A, C)
  B_B <- bc_b_step(Z_B, C)
  list(B_full = B_full, B_A = B_A, B_B = B_B,
       B_jk = 2 * B_full - (B_A + B_B) / 2)
}

#' Entrywise coverage of Bz0 by CIs centred at B (Procrustes-aligned),
#' with the sandwich rotated consistently; optional entrywise variance
#' addition (unaligned frame) for the jackknife-inflated variant.
fs_coverage_entry <- function(B, Bz0, Sigma, P, Kp, var_add = NULL) {
  pa <- procrustes_align(B, Bz0)
  Sig <- Sigma
  if (!is.null(var_add)) Sig <- Sig + diag(as.vector(var_add))
  RkI <- kronecker(t(pa$R), diag(P))
  Sig_rot <- RkI %*% Sig %*% t(RkI)
  se <- matrix(sqrt(pmax(diag(Sig_rot), 0)), P, Kp)
  z_q <- qnorm(0.975)
  covers <- (Bz0 >= pa$Bz_aligned - z_q * se) &
            (Bz0 <= pa$Bz_aligned + z_q * se)
  list(covers = covers, se = se, mse = pa$mse,
       std_err = (pa$Bz_aligned - Bz0) / se)
}

#' Row-norm (alignment-free) coverage via the delta method
fs_coverage_rownorm <- function(B, Bz0, Sigma, P, Kp, var_add = NULL) {
  Sig <- Sigma
  if (!is.null(var_add)) Sig <- Sig + diag(as.vector(var_add))
  s_hat <- rowSums(B^2)
  s_se <- vapply(seq_len(P), function(j) {
    idx <- j + (seq_len(Kp) - 1L) * P
    sqrt(max(4 * B[j, ] %*% Sig[idx, idx] %*% B[j, ], 0))
  }, numeric(1))
  s_true <- rowSums(Bz0^2)
  z_q <- qnorm(0.975)
  (s_true >= s_hat - z_q * s_se) & (s_true <= s_hat + z_q * s_se)
}
