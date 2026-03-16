#' Standard Errors for ILR Path Coefficients
#'
#' Estimates standard errors and confidence intervals for the ILR path
#' coefficient matrix \eqn{\hat{\mathbf{B}}_z} via nonparametric bootstrap
#' resampling of documents.
#'
#' The bootstrap SE is robust to near-degenerate eigenspectra and does not
#' require the eigengap assumption. For large corpora (`M > 5000`) consider
#' reducing `B` for speed.
#'
#' @param fit An `"sgscatm"` object returned by [sgscatm()].
#' @param W Numeric M x N document-term matrix used to fit `fit`.
#' @param C Numeric M x P covariate matrix used to fit `fit`.
#' @param B Integer. Number of bootstrap replicates. Default 200.
#' @param conf Numeric in (0,1). Confidence level for intervals. Default 0.95.
#' @param seed Integer or NULL. Random seed for reproducibility.
#' @param verbose Logical. Print progress every 50 replicates. Default FALSE.
#'
#' @return A list with:
#'   \describe{
#'     \item{se}{P x (K-1) matrix of bootstrap standard errors for `Bz`.}
#'     \item{ci_lower}{P x (K-1) matrix of lower confidence bounds
#'       (percentile method).}
#'     \item{ci_upper}{P x (K-1) matrix of upper confidence bounds.}
#'     \item{Bz_boot}{P x (K-1) x B array of bootstrap replicates.}
#'     \item{B}{Number of replicates used.}
#'     \item{conf}{Confidence level used.}
#'   }
#' @seealso [sgscatm()], [ilr_se_analytical()]
#' @export
ilr_se <- function(fit, W, C, B = 200L, conf = 0.95,
                   seed = NULL, verbose = FALSE) {
  stopifnot(inherits(fit, "sgscatm"))
  W  <- as.matrix(W)
  C  <- as.matrix(C)
  M  <- nrow(W)
  K  <- fit$K
  Kp <- K - 1L
  P  <- ncol(C)
  stopifnot(nrow(C) == M)
  stopifnot(B >= 2L, conf > 0, conf < 1)

  if (!is.null(seed)) set.seed(seed)

  Bz_boot <- array(NA_real_, dim = c(P, Kp, B))

  for (b in seq_len(B)) {
    if (verbose && b %% 50L == 0L)
      message(sprintf("Bootstrap replicate %d / %d", b, B))

    idx <- sample.int(M, M, replace = TRUE)
    fb  <- tryCatch(
      sgscatm(W[idx, , drop = FALSE],
              C[idx, , drop = FALSE],
              K        = K,
              lambda   = fit$lambda,
              r        = nrow(fit$Psi) + P,   # same rank budget
              V        = fit$V,
              scale_W  = fit$scale_W,
              rotate   = fit$rotate),
      error = function(e) NULL
    )
    if (!is.null(fb)) Bz_boot[, , b] <- fb$Bz
  }

  # drop failed replicates
  ok <- apply(Bz_boot, 3L, function(x) all(is.finite(x)))
  if (sum(ok) < 2L) stop("Too many bootstrap replicates failed.")
  Bz_ok <- Bz_boot[, , ok, drop = FALSE]

  se_mat    <- apply(Bz_ok, c(1L, 2L), sd)
  alpha     <- (1 - conf) / 2
  ci_lower  <- apply(Bz_ok, c(1L, 2L), quantile, probs = alpha)
  ci_upper  <- apply(Bz_ok, c(1L, 2L), quantile, probs = 1 - alpha)

  list(
    se        = se_mat,
    ci_lower  = ci_lower,
    ci_upper  = ci_upper,
    Bz_boot   = Bz_ok,
    B         = sum(ok),
    conf      = conf
  )
}


#' Analytical First-Order SE for ILR Path Coefficients (Experimental)
#'
#' Computes the first-order asymptotic covariance of \eqn{\hat{\mathbf{B}}_z}
#' using the eigenvector perturbation formula from the paper. This function
#' is **experimental**: it requires a well-separated eigenspectrum
#' (Assumption 4 in the paper) and uses the stored eigenvectors from
#' the fit object.
#'
#' When eigenvalues are nearly equal the formula is ill-conditioned; a warning
#' is issued and results should be treated with caution. Prefer [ilr_se()] for
#' routine use.
#'
#' @param fit An `"sgscatm"` object returned by [sgscatm()].
#' @param tol_gap Numeric. Relative eigengap threshold. Default 0.05.
#'
#' @return A list with:
#'   \describe{
#'     \item{se}{P x (K-1) matrix of analytical standard errors.}
#'     \item{Sigma_Bz}{P(K-1) x P(K-1) asymptotic covariance matrix.}
#'   }
#' @export
ilr_se_analytical <- function(fit, tol_gap = 0.05) {
  stopifnot(inherits(fit, "sgscatm"))
  if (is.null(fit$U_all))
    stop("fit was created with an old version of sgscatm(); please refit.")

  K   <- fit$K
  Kp  <- K - 1L
  C   <- fit$C_centred          # M x P
  M   <- nrow(C);  P <- ncol(C)

  # Sample eigenvalues/eigenvectors scaled to M^{-1}S_z
  sigma_top <- fit$eigenvalues / M        # K-1 (population-scale)
  U_all     <- fit$U_all                  # M x L
  sigma_all <- fit$eigenvalues_all / M    # L

  # Check eigengap
  if (length(sigma_all) > Kp) {
    gap     <- sigma_top[Kp] - sigma_all[Kp + 1L]
    rel_gap <- gap / max(abs(sigma_top))
    if (rel_gap < tol_gap)
      warning(sprintf(
        "Small relative eigengap (%.4f) at position %d/%d. ",
        rel_gap, Kp, Kp + 1L
      ), "Analytical SEs may be unreliable.")
  }

  W_tilde  <- fit$W_tilde
  Sigma_C  <- crossprod(C) / M
  SigC_inv <- solve(Sigma_C)

  # Precompute P_C * each column of U_all
  CtU_all <- crossprod(C, U_all)                       # P x L
  PC_Uall <- C %*% (SigC_inv %*% CtU_all)             # M x L

  # W~W~' u_l = S_z u_l - lambda*P_C u_l
  # In population-scaled form: (1/M)(W~W~' u_l) = sigma_all[l]*u_l - (lambda/M)*PC_Uall
  # For the influence score we use the unscaled form divided by M at the end
  WW_Uall <- sweep(U_all, 2L, fit$eigenvalues_all, "*") -
    fit$lambda * PC_Uall                               # M x L  (S_z U_all)

  L     <- ncol(U_all)
  Gamma <- matrix(0, P * Kp, P * Kp)
  eps   <- .Machine$double.eps * max(abs(sigma_all)) * 1e4

  for (k in seq_len(Kp)) {
    for (kp in seq_len(Kp)) {
      G_kkp <- matrix(0, P, P)
      for (l in seq_len(L)) {
        d_l <- sigma_top[k] - sigma_all[l]
        if (abs(d_l) < eps) next
        s_lk <- .iscore(U_all[, l], WW_Uall[, k], PC_Uall[, k],
                        fit$lambda, M)
        for (lp in seq_len(L)) {
          d_lp <- sigma_top[kp] - sigma_all[lp]
          if (abs(d_lp) < eps) next
          s_lp <- .iscore(U_all[, lp], WW_Uall[, kp], PC_Uall[, kp],
                          fit$lambda, M)
          G_kkp <- G_kkp +
            (1 / (d_l * d_lp)) *
            crossprod(C * s_lk, C * s_lp) / M
        }
      }
      rk  <- (k  - 1L) * P + seq_len(P)
      rkp <- (kp - 1L) * P + seq_len(P)
      Gamma[rk, rkp] <- G_kkp
    }
  }

  IKp            <- diag(Kp)
  SigC_inv_block <- kronecker(IKp, SigC_inv)
  Sigma_Bz       <- SigC_inv_block %*% Gamma %*% SigC_inv_block
  se_vec         <- sqrt(pmax(diag(Sigma_Bz), 0) / M)

  list(se       = matrix(se_vec, P, Kp),
       Sigma_Bz = Sigma_Bz)
}

# Internal: per-document influence score (M-vector), centred.
# u_l[i] * (S_z u_k)[i] / M,  minus mean.
.iscore <- function(u_l, Sz_uk, PC_uk, lambda, M) {
  # Sz_uk = WW_uk + lambda*PC_uk already passed as S_z u_k (precomputed)
  # We divide by M to match population-scaled formula
  raw <- u_l * Sz_uk / M
  raw - mean(raw)
}
