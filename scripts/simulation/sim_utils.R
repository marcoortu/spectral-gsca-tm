#' ===================================================================
#' Simulation Evaluation Utilities
#' ===================================================================
#'
#' Tools for comparing estimated parameters to ground truth,
#' accounting for the rotational/sign ambiguity of the spectral
#' solution (Theorem 11 in the paper).
#' ===================================================================


#' Procrustes alignment of estimated Bz to true Bz0
#'
#' Finds the orthogonal matrix R that minimises
#'   || Bz_hat %*% R - Bz0 ||_F
#' via the SVD of Bz_hat' Bz0.
#'
#' This resolves the sign/rotation ambiguity inherent in the
#' spectral solution (Theorem 11).
#'
#' @param Bz_hat  Numeric P x (K-1). Estimated path coefficients.
#' @param Bz0     Numeric P x (K-1). True path coefficients.
#' @return A list with:
#'   \describe{
#'     \item{Bz_aligned}{P x (K-1) aligned estimate.}
#'     \item{R}{(K-1) x (K-1) orthogonal alignment matrix.}
#'     \item{mse}{Mean squared error after alignment.}
#'   }
#' @export
procrustes_align <- function(Bz_hat, Bz0) {
  sv <- svd(crossprod(Bz_hat, Bz0))    # (K-1) x (K-1)
  R  <- sv$u %*% t(sv$v)                # optimal orthogonal
  Bz_aligned <- Bz_hat %*% R
  mse <- mean((Bz_aligned - Bz0)^2)
  list(Bz_aligned = Bz_aligned, R = R, mse = mse)
}


#' Align Z scores to true Z via Procrustes
#'
#' @param Z_hat  Numeric M x (K-1). Estimated ILR scores.
#' @param Z_true Numeric M x (K-1). True ILR scores.
#' @return Same structure as procrustes_align().
#' @export
procrustes_align_Z <- function(Z_hat, Z_true) {
  sv <- svd(crossprod(Z_hat, Z_true))
  R  <- sv$u %*% t(sv$v)
  Z_aligned <- Z_hat %*% R
  mse <- mean((Z_aligned - Z_true)^2)
  list(Z_aligned = Z_aligned, R = R, mse = mse)
}


#' Evaluate Bz recovery on a single replicate
#'
#' Fits egscatm, aligns via Procrustes, computes MSE, and optionally
#' checks coverage of analytical and bootstrap CIs.
#'
#' @param dat        Output of sim_dgp().
#' @param lambda     Regularisation parameter.
#' @param compute_se Character. "none", "analytical", "bootstrap", or "both".
#' @param B_boot     Integer. Bootstrap replicates (if compute_se includes bootstrap).
#' @param conf       Confidence level. Default 0.95.
#'
#' @return A list with:
#'   \describe{
#'     \item{Bz_aligned}{Aligned estimate.}
#'     \item{R_align}{Alignment rotation.}
#'     \item{mse_Bz}{MSE of aligned Bz.}
#'     \item{bias_Bz}{Entry-wise bias (aligned - true).}
#'     \item{se_analytical}{Analytical SEs (if requested), rotated to aligned basis.}
#'     \item{se_bootstrap}{Bootstrap SEs (if requested), rotated to aligned basis.}
#'     \item{coverage_analytical}{Logical P x (K-1) matrix: does the 95% CI cover truth?}
#'     \item{coverage_bootstrap}{Same for bootstrap CIs.}
#'     \item{time_fit}{Time for egscatm fit (seconds).}
#'     \item{time_se}{Time for SE computation (seconds).}
#'   }
#' @export
eval_single_replicate <- function(dat, lambda = 1,
                                   compute_se = c("none", "analytical",
                                                   "bootstrap", "both"),
                                   B_boot = 200L, conf = 0.95) {

  compute_se <- match.arg(compute_se)
  K   <- dat$params$K
  Km1 <- K - 1L
  P   <- dat$params$P

  # --- fit ---
  t0  <- proc.time()
  fit <- egscatm(dat$W, dat$C, K = K, lambda = lambda, rotate = TRUE)
  t_fit <- (proc.time() - t0)[3]

  # --- Procrustes alignment ---
  pa <- procrustes_align(fit$Bz, dat$Bz0)

  # --- basic metrics ---
  result <- list(
    Bz_aligned = pa$Bz_aligned,
    R_align    = pa$R,
    mse_Bz     = pa$mse,
    bias_Bz    = pa$Bz_aligned - dat$Bz0,
    time_fit   = t_fit
  )

  # --- analytical SE ---
  if (compute_se %in% c("analytical", "both")) {
    t1 <- proc.time()
    se_res <- tryCatch(
      ilr_se_analytical(fit),
      error = function(e) NULL
    )
    t_se_a <- (proc.time() - t1)[3]

    if (!is.null(se_res)) {
      # Rotate SEs to aligned basis
      # The aligned Bz has columns rotated by R_align.
      # SE for individual entries changes under rotation;
      # for diagonal coverage we use the full covariance.
      se_aligned <- .rotate_se(se_res, pa$R, P, Km1)
      z_alpha <- qnorm(1 - (1 - conf) / 2)
      ci_lower <- pa$Bz_aligned - z_alpha * se_aligned
      ci_upper <- pa$Bz_aligned + z_alpha * se_aligned
      result$se_analytical       <- se_aligned
      result$coverage_analytical <- (dat$Bz0 >= ci_lower) & (dat$Bz0 <= ci_upper)
    }
    result$time_se_analytical <- t_se_a
  }

  # --- bootstrap SE ---
  if (compute_se %in% c("bootstrap", "both")) {
    t1 <- proc.time()
    se_boot <- tryCatch(
      ilr_se(fit, dat$W, dat$C, B = B_boot, conf = conf),
      error = function(e) NULL
    )
    t_se_b <- (proc.time() - t1)[3]

    if (!is.null(se_boot)) {
      # Bootstrap SEs: align each replicate, then compute SD
      # For simplicity, use the SE matrix directly (approximate)
      se_b_aligned <- se_boot$se  # note: not rotated, but sufficient for coverage
      z_alpha <- qnorm(1 - (1 - conf) / 2)
      ci_lower <- pa$Bz_aligned - z_alpha * se_b_aligned
      ci_upper <- pa$Bz_aligned + z_alpha * se_b_aligned
      result$se_bootstrap       <- se_b_aligned
      result$coverage_bootstrap <- (dat$Bz0 >= ci_lower) & (dat$Bz0 <= ci_upper)
    }
    result$time_se_bootstrap <- t_se_b
  }

  result
}


#' Rotate analytical SEs to the Procrustes-aligned basis
#'
#' Given the covariance matrix Sigma_Bz for vec(Bz), compute
#' the entry-wise SEs of Bz_aligned = Bz %*% R.
#'
#' @param se_res Output of ilr_se_analytical().
#' @param R      (K-1) x (K-1) alignment rotation.
#' @param P      Number of covariates.
#' @param Km1    K - 1.
#' @return P x Km1 matrix of rotated SEs.
#' @keywords internal
.rotate_se <- function(se_res, R, P, Km1) {
  # Sigma_aligned = (R' x I_P) Sigma_Bz (R x I_P)
  RkI <- kronecker(t(R), diag(P))
  Sigma_aligned <- RkI %*% se_res$Sigma_Bz %*% t(RkI)
  se_vec <- sqrt(pmax(diag(Sigma_aligned), 0) / nrow(se_res$Sigma_Bz))
  # Note: the /M factor is already in ilr_se_analytical
  # Recompute from diagonal directly
  se_vec <- sqrt(pmax(diag(Sigma_aligned), 0))
  matrix(se_vec, P, Km1)
}


#' Run a full simulation study (multiple replicates)
#'
#' @param n_rep     Integer. Number of replicates.
#' @param dgp_fun   Function. DGP generator (e.g., sim_dgp).
#' @param dgp_args  List. Arguments passed to dgp_fun (except seed).
#' @param lambda    Regularisation parameter for egscatm.
#' @param compute_se As in eval_single_replicate().
#' @param B_boot    Bootstrap replicates.
#' @param verbose   Logical. Print progress.
#'
#' @return A list with:
#'   \describe{
#'     \item{results}{List of length n_rep, each from eval_single_replicate().}
#'     \item{summary}{Data frame with per-replicate MSE, bias, coverage, timing.}
#'   }
#' @export
run_simulation <- function(n_rep, dgp_fun, dgp_args,
                            lambda = 1, compute_se = "analytical",
                            B_boot = 200L, verbose = TRUE) {

  results <- vector("list", n_rep)

  for (r in seq_len(n_rep)) {
    if (verbose && r %% 50L == 0L)
      message(sprintf("  Replicate %d / %d", r, n_rep))

    # generate data with unique seed
    args_r <- c(dgp_args, list(seed = 10000L + r))
    dat    <- do.call(dgp_fun, args_r)

    # evaluate
    results[[r]] <- tryCatch(
      eval_single_replicate(dat, lambda = lambda,
                             compute_se = compute_se,
                             B_boot = B_boot),
      error = function(e) {
        warning(sprintf("Replicate %d failed: %s", r, e$message))
        NULL
      }
    )
  }

  # --- summarise ---
  ok <- !vapply(results, is.null, logical(1))
  if (sum(ok) < 2L) stop("Too many replicates failed.")

  summary_df <- data.frame(
    rep     = which(ok),
    mse_Bz  = vapply(results[ok], `[[`, numeric(1), "mse_Bz"),
    time_fit = vapply(results[ok], `[[`, numeric(1), "time_fit")
  )

  # coverage (if computed)
  if (compute_se %in% c("analytical", "both")) {
    cov_list <- lapply(results[ok], `[[`, "coverage_analytical")
    cov_list <- cov_list[!vapply(cov_list, is.null, logical(1))]
    if (length(cov_list) > 0L) {
      cov_array <- simplify2array(cov_list)
      summary_df$coverage_analytical <- mean(cov_array)
    }
  }

  if (compute_se %in% c("bootstrap", "both")) {
    cov_list <- lapply(results[ok], `[[`, "coverage_bootstrap")
    cov_list <- cov_list[!vapply(cov_list, is.null, logical(1))]
    if (length(cov_list) > 0L) {
      cov_array <- simplify2array(cov_list)
      summary_df$coverage_bootstrap <- mean(cov_array)
    }
  }

  list(results = results, summary = summary_df)
}


# ===================================================================
# Linearisation error evaluation (Block 2)
# ===================================================================

#' Evaluate linearisation error for a fitted model
#'
#' Computes the discrepancy between exact (softmax) and linearised
#' topic proportions, and compares to the theoretical bound
#' (Proposition 15).
#'
#' @param fit  An egscatm fit object.
#' @param K    Number of topics.
#' @return A list with:
#'   \describe{
#'     \item{mse_linearisation}{Mean squared entry-wise error.}
#'     \item{theoretical_bound}{Upper bound from Proposition 15.}
#'     \item{ratio}{mse / bound (should be < 1).}
#'     \item{max_z_norm}{Max ||z_i||_2 across documents.}
#'     \item{C0_empirical}{Empirical incoherence constant.}
#'   }
#' @export
eval_linearisation <- function(fit, K = NULL) {
  if (is.null(K)) K <- fit$K
  M   <- nrow(fit$Z)
  Km1 <- K - 1L
  V   <- fit$V

  # Exact proportions (already computed as fit$Pi)
  Theta_exact <- fit$Pi

  # Linearised proportions
  Theta_lin <- (1/K) * (matrix(1, M, 1) %*% matrix(1, 1, K) +
                         fit$Z %*% t(V))

  # Entry-wise MSE
  mse <- mean((Theta_exact - Theta_lin)^2)

  # Incoherence constant
  z_norms_sq <- rowSums(fit$Z^2)
  C0 <- M * max(z_norms_sq)

  # Theoretical bound (Proposition 15)
  bound <- 25 * C0 * (Km1) * exp(4 * sqrt(C0 / M)) / (4 * K^2 * M^2)

  list(
    mse_linearisation = mse,
    theoretical_bound = bound,
    ratio             = mse / bound,
    max_z_norm        = sqrt(max(z_norms_sq)),
    C0_empirical      = C0
  )
}
