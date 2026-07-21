## =====================================================================
## estimators.R  —  EXPERIMENTAL one-step-efficient estimators of B_z0
## =====================================================================
##
## Reuses (does NOT reimplement / does NOT modify) the package:
##   sgscatm(), sgscatm_vcov(), ilr_se(), ilr_contrast(), and
##   replication/simulation/sim_utils.R::procrustes_align().
##
## Conventions are fixed in replication/one_step/preregistration.md §1.
## Every symbol here matches sgscatm_fit.R / sim_dgp.R.
##
## Estimators implemented (all return B_hat_z (P x (K-1)), an HC0 SE
## matrix, the recovered ILR scores H (M x (K-1)), and the raw regression
## vcov of vec(B) in column-major order):
##   (A) baseline_std  : current standardized-eigenvector-score estimator
##   (B) proj          : loading-projection scores (natural scale)
##   (C) onestep_uw    : one UNWEIGHTED Gauss-Newton closure step from proj
##   (D) onestep_mw    : one MULTINOMIAL-WEIGHTED Gauss-Newton step from proj
## =====================================================================


## ---- small linear-algebra helpers -----------------------------------

# Symmetric (pseudo-)inverse via eigen-thresholding.
.ginv_sym <- function(S, tol = 1e-8) {
  eg <- eigen(S, symmetric = TRUE)
  d  <- eg$values
  keep <- d > tol * max(abs(d))
  if (!any(keep)) return(matrix(0, nrow(S), ncol(S)))
  eg$vectors[, keep, drop = FALSE] %*%
    (t(eg$vectors[, keep, drop = FALSE]) / d[keep])
}

# Guarded solve of a small SPD system, optional ridge.
.solve_ridge <- function(A, b, ridge = 0) {
  if (ridge > 0) A <- A + diag(ridge, nrow(A))
  ch <- tryCatch(chol(A), error = function(e) NULL)
  if (!is.null(ch)) return(chol2inv(ch) %*% b)
  .ginv_sym(A) %*% b
}


## ---- HC0 sandwich on the final multivariate ILR regression ----------
##
## Model:  H = C B + E,  B = (C'C)^{-1} C' H,  e_i = H_i - B' c_i.
## Score g_i = e_i (x) c_i (column-major), bread = -(I_{K-1} (x) C'C),
## meat = sum_i (e_i e_i') (x) (c_i c_i'). Returns vcov of vec(B) in
## column-major order (block k stacks the P coefficients of column k),
## matching kronecker(t(R), diag(P)) rotation used by sgscatm_vcov().
## Optional per-document weights `w` (length M) give WLS:
##   B = (C'WC)^{-1} C'W H,  with the HC sandwich carrying the weights
##   (meat block (k,l) = C' diag(w_i^2 e_ik e_il) C).
.hc_ilr_sandwich <- function(C, H, w = NULL) {
  M <- nrow(C); P <- ncol(C); Km1 <- ncol(H)
  if (is.null(w)) w <- rep(1, M)
  Cw     <- C * w                            # rows scaled by w_i
  CtWC   <- crossprod(C, Cw)                  # C' W C
  CtWCinv <- solve(CtWC)
  B <- CtWCinv %*% crossprod(Cw, H)          # P x (K-1)
  E <- H - C %*% B                            # M x (K-1)

  # meat = sum_i w_i^2 (e_i e_i') (x) (c_i c_i')
  w2 <- w^2
  meat <- matrix(0, P * Km1, P * Km1)
  for (k in seq_len(Km1)) for (l in seq_len(Km1)) {
    blk <- crossprod(C, C * (w2 * E[, k] * E[, l]))
    rk <- (k - 1L) * P + seq_len(P); rl <- (l - 1L) * P + seq_len(P)
    meat[rk, rl] <- blk
  }
  bread_inv <- kronecker(diag(Km1), CtWCinv)
  V <- bread_inv %*% meat %*% bread_inv     # HC0 vcov of vec(B), column-major
  se <- matrix(sqrt(pmax(diag(V), 0)), P, Km1)
  list(B = B, vcov = V, se = se, resid = E)
}


## ---- reconstruct the (row-normalised) term-frequency matrix ---------
##
## sgscatm(scale_W=TRUE) row-normalises W then centres:
##   W_tilde = F - 1 w_bar',  so  F = W_tilde + 1 w_bar'.
.reconstruct_F <- function(fit) {
  sweep(fit$W_tilde, 2L, fit$w_bar, "+")    # M x N term frequencies
}


## ---- proportion estimate by topic-word right-inverse projection -----
##
## pi_i = Beta' theta_i  is EXACTLY linear in theta_i; with a NATURAL-SCALE
## topic-word anchor Phi_anchor (K x N, rows sum to 1) the ordinary right
## inverse
##   theta_hat_i = (Phi Phi')^{-1} Phi f_i
## returns theta_i exactly when f_i = Phi' theta_i, at natural scale, for
## any document length (multinomial noise is mean-zero, O(L^{-1/2})).
##
## IMPORTANT (deviation logged in report.md / prereg addendum): the fit's
## OWN Phi (fit$Phi) has the correct row-space but is scale-inflated by the
## unit-norm score normalisation (Z*'Z*=I), so projecting on it reproduces
## the collapsed baseline scale. The loading-projection therefore requires
## an EXTERNAL natural-scale topic-word anchor. In simulation this is the
## generative Beta (a consistent topic-word estimate); G6 degrades it to
## characterise sensitivity to anchor error. If Phi_anchor is NULL we fall
## back to fit$Phi (used ONLY to show the collapse in the G1 diagnosis).
.proportions_hat <- function(fit, Phi_anchor = NULL) {
  Fmat <- .reconstruct_F(fit)               # M x N
  Phi  <- if (is.null(Phi_anchor)) fit$Phi else Phi_anchor   # K x N
  PPt  <- tcrossprod(Phi)                    # K x K
  # theta_hat = F Phi' (Phi Phi')^{-1}
  Fmat %*% t(Phi) %*% .ginv_sym(PPt)         # M x K
}


## =====================================================================
## Estimator (A): baseline standardized-eigenvector-score estimator
## =====================================================================
est_baseline_std <- function(fit) {
  H  <- fit$Z                                # M x (K-1), Z*'Z* = I (unit norm)
  sw <- .hc_ilr_sandwich(fit$C_centred, H)
  list(name = "baseline_std", B = fit$Bz, se = sw$se, vcov = sw$vcov, H = H)
}


## =====================================================================
## Estimator (B): loading-projection scores (natural scale)
## =====================================================================
##  H_proj = K * theta_hat %*% V   (anchored at the centroid p = 1/K,
##  since V'1 = 0 makes any constant proportion offset irrelevant).
est_proj <- function(fit, theta_hat = NULL, Phi_anchor = NULL) {
  if (is.null(theta_hat)) theta_hat <- .proportions_hat(fit, Phi_anchor)
  V <- fit$V; K <- fit$K
  H <- K * theta_hat %*% V                   # M x (K-1)  natural-scale scores
  sw <- .hc_ilr_sandwich(fit$C_centred, H)
  list(name = "proj", B = sw$B, se = sw$se, vcov = sw$vcov,
       H = H, theta_hat = theta_hat)
}


## =====================================================================
## Estimators (C)/(D): one Gauss-Newton closure step from the proj start
## =====================================================================
##  Re-linearise p_i = softmax(V eta) about eta_hat_i (NOT the centroid):
##    s_i = softmax(V eta_hat_i),  J_i = (diag(s_i) - s_i s_i') V,
##    r_i = theta_hat_i - s_i,
##    uw:  delta_i = (J_i'J_i + ridge)^{-1} J_i' r_i
##    mw:  Omega_i = (1/L_i)(diag(th) - th th');  Wi = ginv(Omega_i),
##         delta_i = (J_i' Wi J_i + ridge)^{-1} J_i' Wi r_i
##  Exactly ONE step; barrier-robust via ridge + pseudo-inverse.
.softmax <- function(x) { x <- x - max(x); e <- exp(x); e / sum(e) }

## final_weight: pooled-regression weighting of the per-document scores.
##   "none"      -> OLS  (onestep_uw / onestep_mw)
##   "L"         -> WLS with w_i = L_i (multinomial precision proxy;
##                  downweights short/noisy documents). This is where the
##                  multinomial weight actually bites — the per-document GN
##                  weight is a scalar that cancels (see report G5c).
est_onestep <- function(fit, L, weighted = TRUE, proj = NULL,
                        Phi_anchor = NULL, ridge_rel = 1e-6,
                        final_weight = c("none", "L")) {
  final_weight <- match.arg(final_weight)
  if (is.null(proj)) proj <- est_proj(fit, Phi_anchor = Phi_anchor)
  V <- fit$V; K <- fit$K; Km1 <- K - 1L
  M <- nrow(proj$H)
  th <- proj$theta_hat
  Hnew <- matrix(0, M, Km1)
  if (length(L) == 1L) L <- rep(L, M)

  for (i in seq_len(M)) {
    eta <- proj$H[i, ]
    s   <- .softmax(as.numeric(V %*% eta))            # K
    J   <- (diag(s) - tcrossprod(s)) %*% V            # K x (K-1)
    r   <- th[i, ] - s                                 # K
    JtJ <- crossprod(J)
    ridge <- ridge_rel * (sum(diag(JtJ)) / Km1 + 1e-12)
    if (weighted) {
      thc <- pmax(th[i, ], 1e-8); thc <- thc / sum(thc)  # clamp for the working cov only
      Omega <- (diag(thc) - tcrossprod(thc)) / L[i]
      Wi    <- .ginv_sym(Omega)
      A <- crossprod(J, Wi) %*% J
      bb <- crossprod(J, Wi) %*% r
    } else {
      A  <- JtJ
      bb <- crossprod(J, r)
    }
    delta <- .solve_ridge(A, bb, ridge = ridge)
    Hnew[i, ] <- eta + as.numeric(delta)
  }

  wt <- if (final_weight == "L") L / mean(L) else NULL
  sw <- .hc_ilr_sandwich(fit$C_centred, Hnew, w = wt)
  nm <- if (final_weight == "L") "onestep_wls"
        else if (weighted) "onestep_mw" else "onestep_uw"
  list(name = nm, B = sw$B, se = sw$se, vcov = sw$vcov, H = Hnew)
}


## =====================================================================
## Driver: fit once, run all four estimators
## =====================================================================
sg_all_estimators <- function(W, C, K, L, lambda = 1, V = NULL,
                              rotate = TRUE, Phi_anchor = NULL,
                              which = c("baseline_std", "proj",
                              "onestep_uw", "onestep_mw")) {
  fit <- sgscatm(W, C, K = K, lambda = lambda, V = V, rotate = rotate)
  out <- list(fit = fit)
  proj <- NULL
  if ("baseline_std" %in% which) out$baseline_std <- est_baseline_std(fit)
  if (any(c("proj","onestep_uw","onestep_mw","onestep_wls") %in% which)) {
    proj <- est_proj(fit, Phi_anchor = Phi_anchor)
    if ("proj" %in% which) out$proj <- proj
  }
  if ("onestep_uw" %in% which)
    out$onestep_uw <- est_onestep(fit, L, weighted = FALSE, proj = proj)
  if ("onestep_mw" %in% which)
    out$onestep_mw <- est_onestep(fit, L, weighted = TRUE, proj = proj)
  if ("onestep_wls" %in% which)
    out$onestep_wls <- est_onestep(fit, L, weighted = FALSE, proj = proj,
                                   final_weight = "L")
  out
}


## =====================================================================
## Evaluation: Procrustes-align to B_z0, rotate the sandwich vcov, and
## report direct + rotation-aligned RMSE, per-entry SE, coverage.
## =====================================================================
##  Direct RMSE uses B_z0's own scale; rotation-aligned RMSE minimises
##  over the orthogonal group (procrustes_align). If the two ~coincide,
##  the estimator recovered SCALE, not merely rotation.
eval_vs_Bz0 <- function(est, Bz0, conf = 0.95) {
  P <- nrow(Bz0); Km1 <- ncol(Bz0)
  pa <- procrustes_align(est$B, Bz0)          # from sim_utils.R
  Ba <- pa$Bz_aligned; R <- pa$R

  # rotate HC vcov to the aligned basis: (R' (x) I_P) V (R (x) I_P)
  Kk <- kronecker(t(R), diag(P))
  Vr <- Kk %*% est$vcov %*% t(Kk)
  # The Procrustes-aligned estimator has NO variance along the rotational
  # tangent (Bz0 is identified only up to rotation); project it out, exactly
  # as sgscatm_vcov(identified=TRUE) does, so the SE matches the aligned SD.
  Pr <- getFromNamespace(".rot_tangent_projector", "sgscatm")(Ba)
  Vr <- Pr %*% Vr %*% Pr
  se_al <- matrix(sqrt(pmax(diag(Vr), 0)), P, Km1)

  rmse_dir <- sqrt(mean((Ba - Bz0)^2))        # aligned = direct here (Ba on Bz0 scale)
  # "rotation-only-aligned" reference: best orthogonal fit of the RAW estimate,
  # already what procrustes gives; the scale check compares this to the
  # unrotated direct error of the raw estimate on the Bz0 scale.
  rmse_raw <- sqrt(mean((est$B - Bz0)^2))

  z <- qnorm(1 - (1 - conf) / 2)
  cover <- (Bz0 >= Ba - z * se_al) & (Bz0 <= Ba + z * se_al)

  list(B_aligned = Ba, R = R, se_aligned = se_al,
       rmse = rmse_dir,               # RMSE(B_z0) after Procrustes (primary)
       rmse_raw = rmse_raw,           # RMSE of the raw (unrotated) estimate
       bias = Ba - Bz0,
       coverage = cover,
       norm_Bz0 = sqrt(mean(Bz0^2)))
}
