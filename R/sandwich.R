#' Heteroskedasticity-robust sandwich covariance for the chain path coefficients
#'
#' Implements the Lemma-17 sandwich covariance of \eqn{\mathrm{vec}(\hat B_z)}
#' for the multivariate regression of the refined ILR scores \eqn{\hat Z} on the
#' covariates \eqn{C}. With refined residual rows
#' \eqn{\hat r_i = \hat z_i - \hat B_z^\top c_i},
#' \deqn{\hat\Sigma_B = \tfrac{M}{M-P}\,(I_{K-1}\otimes(C^\top C)^{-1})
#'   \Big[\textstyle\sum_i (\hat r_i\hat r_i^\top)\otimes(c_ic_i^\top)\Big]
#'   (I_{K-1}\otimes(C^\top C)^{-1}).}
#' Promoted verbatim from `replication/feasibility/03_jackknife.R::fs_sandwich`
#' (unit-tested against direct simulation).
#'
#' The covariance estimates the sampling variance of \eqn{\hat B_z} around its
#' own mean, which includes the \eqn{\sqrt M/L} incidental-parameter drift and
#' the \eqn{\sqrt M\,\delta_\Phi} orientation drift; therefore Wald coverage of
#' the population \eqn{B_{z,0}} is nominal only when both drifts vanish
#' (Corollary 18), while the variance itself stays calibrated regardless.
#'
#' @param object An `"sgscatm_chain"` object from [sgscatm_chain()].
#' @param ... Unused.
#' @return A \eqn{P(K-1)\times P(K-1)} covariance matrix (column-major
#'   `vec(Bz)` order: covariate fastest within component).
#' @seealso [sgscatm_chain()], [ilr_se()], [sgscatm_vcov()]
#' @export
vcov.sgscatm_chain <- function(object, ...) {
  .sg_sandwich(object$Z, object$C_centred, object$Bz)
}

#' Full-chain document-bootstrap standard errors for the chain path coefficients
#'
#' The calibrated **feasible** SE for `sgscatm_chain()`. Unlike the plain Lemma-17
#' sandwich ([vcov.sgscatm_chain()]) — which is calibrated only when the anchor
#' orientation error \eqn{\delta_\Phi} is negligible — this resamples the M
#' documents with replacement and **re-runs the entire chain** (anchors,
#' orientation, frozen-\eqn{\Phi} refinement) on each resample, so the anchor
#' re-estimation variance that the sandwich (and the split-document jackknife)
#' omit is captured. Each bootstrap \eqn{\hat B_z^\ast} is permutation+sign
#' aligned to the **point estimate** (not the truth) before the per-entry SD is
#' taken.
#'
#' @param fit An `"sgscatm_chain"` point-estimate object.
#' @param W,C The original count DTM and covariate matrix the fit was built from.
#' @param B Integer bootstrap replicates. Default 200.
#' @param conf Confidence level. Default 0.95.
#' @param seed Integer or NULL.
#' @param max_sweeps,dz_cap Passed through to the per-replicate
#'   `sgscatm_chain(refine = "frozen_phi", ...)`.
#' @param verbose Logical; progress every 50 replicates.
#' @return A list with `se` (P x (K-1)), `ci_lower`, `ci_upper`, `B` (successful
#'   replicates), and `boot` (P x (K-1) x B aligned replicates).
#' @seealso [sgscatm_chain()], [vcov.sgscatm_chain()], [perm_sign_align()]
#' @export
chain_boot_se <- function(fit, W, C, B = 200L, conf = 0.95, seed = NULL,
                          max_sweeps = 100L, dz_cap = 1, verbose = FALSE) {
  stopifnot(inherits(fit, "sgscatm_chain"))
  W <- as.matrix(W); C <- as.matrix(C)
  M <- nrow(W); K <- fit$K; V <- fit$V
  Bref <- fit$Bz; P <- nrow(Bref); Km1 <- ncol(Bref)
  if (!is.null(seed)) set.seed(seed)
  boot <- array(NA_real_, dim = c(P, Km1, B))
  for (b in seq_len(B)) {
    if (verbose && b %% 50L == 0L) message(sprintf("boot %d/%d", b, B))
    idx <- sample.int(M, M, replace = TRUE)
    fb <- tryCatch(
      sgscatm_chain(W[idx, , drop = FALSE], C[idx, , drop = FALSE], K = K,
                    refine = "frozen_phi", max_sweeps = max_sweeps,
                    dz_cap = dz_cap, V = V),
      error = function(e) NULL)
    if (!is.null(fb)) boot[, , b] <- perm_sign_align(fb$Bz, Bref, V)$B
  }
  ok <- apply(boot, 3L, function(x) all(is.finite(x)))
  if (sum(ok) < 2L) stop("Too many bootstrap chain fits failed.")
  bok <- boot[, , ok, drop = FALSE]
  se     <- apply(bok, c(1L, 2L), stats::sd)
  se_mad <- apply(bok, c(1L, 2L), stats::mad)   # robust scale (1.4826 x MAD)
  z   <- stats::qnorm(1 - (1 - conf) / 2)
  list(se = se, se_mad = se_mad,
       ci_lower = Bref - z * se, ci_upper = Bref + z * se,
       B = sum(ok), boot = bok)
}

#' Standard errors and Wald intervals for the chain path coefficients
#'
#' @param object An `"sgscatm_chain"` object.
#' @param conf Confidence level for the Wald intervals. Default 0.95.
#' @return A list with `se` (P x (K-1)), `ci_lower`, `ci_upper`, and the full
#'   `vcov`.
#' @export
chain_se <- function(object, conf = 0.95) {
  stopifnot(inherits(object, "sgscatm_chain"))
  Sig <- vcov.sgscatm_chain(object)
  P <- nrow(object$Bz); Km1 <- ncol(object$Bz)
  se <- matrix(sqrt(pmax(diag(Sig), 0)), P, Km1)
  z <- stats::qnorm(1 - (1 - conf) / 2)
  list(se = se, ci_lower = object$Bz - z * se, ci_upper = object$Bz + z * se,
       vcov = Sig)
}

# Lemma-17 sandwich (internal), verbatim from fs_sandwich.
.sg_sandwich <- function(Z, C, B_hat) {
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

# All K! permutations of 1..K (verbatim from fs_all_perms)
.sg_all_perms <- function(K) {
  if (K == 1L) return(matrix(1L, 1L, 1L))
  sub <- .sg_all_perms(K - 1L)
  do.call(rbind, lapply(seq_len(K), function(pos) {
    left  <- if (pos > 1L) sub[, seq_len(pos - 1L), drop = FALSE] else
      matrix(0L, nrow(sub), 0L)
    right <- if (pos <= K - 1L) sub[, pos:(K - 1L), drop = FALSE] else
      matrix(0L, nrow(sub), 0L)
    cbind(left, rep(K, nrow(sub)), right)
  }))
}

#' Permutation + sign alignment of an estimate to a reference
#'
#' Aligns `B_hat` (P x (K-1)) to `B_ref` by the topic-relabeling action only:
#' a K-topic permutation acts on ILR coordinates as the orthogonal
#' \eqn{Q_P = V^\top P V}. Returns the relabeled estimate minimising entrywise
#' MSE and the acting orthogonal map `Q`. **No** continuous Procrustes and **no**
#' general-linear map to the truth is used — these would absorb a degree of
#' freedom the feasible estimator does not have. Enumerates K! (intended for the
#' small K of the certification designs).
#'
#' @param B_hat,B_ref P x (K-1) matrices.
#' @param V K x (K-1) ILR contrast matrix.
#' @return List with `B` (aligned), `Q` (the (K-1)x(K-1) orthogonal map), `mse`.
#' @export
perm_sign_align <- function(B_hat, B_ref, V) {
  K <- nrow(V)
  perms <- .sg_all_perms(K)
  best <- Inf; bestQ <- diag(ncol(V)); bestB <- B_hat
  for (r in seq_len(nrow(perms))) {
    P <- diag(K)[perms[r, ], , drop = FALSE]
    Qp <- crossprod(V, P %*% V)                # (K-1)x(K-1) orthogonal
    cand <- B_hat %*% t(Qp)
    m <- mean((cand - B_ref)^2)
    if (m < best) { best <- m; bestQ <- t(Qp); bestB <- cand }
  }
  list(B = bestB, Q = bestQ, mse = best)
}

#' Permutation+sign entrywise Wald coverage of a reference (simulation only)
#'
#' Aligns `B_hat` to `B_ref` (= truth) by permutation+sign, rotates the
#' sandwich `Sigma` consistently, and returns the entrywise 95% Wald coverage
#' indicator and SEs. This is the honest coverage metric for the chain (no
#' Procrustes).
#' @param B_hat,B_ref P x (K-1). @param Sigma P(K-1) square sandwich.
#' @param V K x (K-1). @param conf confidence level.
#' @return list(covers, se, mse).
#' @export
perm_sign_coverage <- function(B_hat, B_ref, Sigma, V, conf = 0.95) {
  al <- perm_sign_align(B_hat, B_ref, V)
  P <- nrow(B_hat); Km1 <- ncol(B_hat)
  RkI <- kronecker(al$Q, diag(P))             # vec(B Q) = (Q' x I) vec B; Q=al$Q
  Sig_rot <- RkI %*% Sigma %*% t(RkI)
  se <- matrix(sqrt(pmax(diag(Sig_rot), 0)), P, Km1)
  z <- stats::qnorm(1 - (1 - conf) / 2)
  covers <- (B_ref >= al$B - z * se) & (B_ref <= al$B + z * se)
  list(covers = covers, se = se, mse = al$mse)
}
