#' Fit the three-stage anchored k-step spectral estimator (Theorem-16 chain)
#'
#' The manuscript's deliverable estimator (Algorithm 1): a spectral pilot, an
#' anchored **general-linear** orientation, and a fixed budget of **least-squares**
#' Gauss-Newton refinement sweeps with a B-stationarity stop. Its inference is the
#' heteroskedasticity-robust sandwich of Lemma 17 ([vcov.sgscatm_chain()]).
#'
#' This differs from [sgscatm()], which returns only the raw spectral pilot whose
#' path-coefficient block `Bz = solve(C'C, C'Z*)` is the degenerate, scale-collapsed
#' regression of Corollary 11. Use `sgscatm_chain()` for coefficient estimation and
#' inference; the raw pilot is retained here as `pilot_Bz` / `pilot_Z` for the
#' subspace-recovery / collapse diagnostics.
#'
#' @details
#' Stages:
#' \enumerate{
#'   \item **Pilot** — [sgscatm()] internals with `rotate = FALSE`: unit-norm ILR
#'     score subspace \eqn{Z^\ast}.
#'   \item **Anchored orientation** — anchor-word recovery of \eqn{\hat\Phi}
#'     (Section 3.6(i), [.sg_anchor_pipeline]) with an optional pooled
#'     multinomial-EM polish (the single deviance step); per-document Gauss-Newton
#'     read-out on \eqn{\lVert w_i - \hat\Phi^\top f(z)\rVert^2} (trust-region
#'     capped); **general-linear** alignment of the pilot to the read-out
#'     ([.sg_gl_align]).
#'   \item **k-step LS refinement** — damped per-document Gauss-Newton on the exact
#'     \eqn{\lVert w_i - \hat\Phi^\top f(z)\rVert^2} at \eqn{\lambda=0} from the
#'     oriented pilot, with the B-stationarity stop
#'     (`max|dBz|/rms(Bz) < rule_tol` for `patience` consecutive sweeps, cap
#'     `max_sweeps`).
#' }
#' `refine = "joint"` (default) re-estimates \eqn{\Phi} (unconstrained ridge LS)
#' each sweep — the empirically-best "V4" variant of `replication/feasibility/`,
#' stable and monotonically improving, and the estimator the manuscript's reported
#' feasible numbers come from (Prop-19 alternating mode). `refine = "frozen_phi"`
#' holds \eqn{\hat\Phi} fixed — the strict regime for which the Lemma-17 sandwich is
#' derived — but **empirically confirms Proposition 21**: the fixed-\eqn{\Phi}
#' criterion optimum is inconsistent, so B̂z is best after ~1 sweep and drifts away
#' if refined further; if used, apply the sandwich at an early-stopped iterate
#' (`max_sweeps` small). See `replication/certification/CHANGES.md`.
#'
#' @param W Numeric M x N non-negative document-term **count** matrix.
#' @param C Numeric M x P covariate matrix (column-centred internally).
#' @param K Integer >= 2, number of topics.
#' @param lambda Pilot regularisation (Stage 1 only; the refinement uses
#'   \eqn{\lambda=0}). Default 1.
#' @param refine `"frozen_phi"` (default, theory) or `"joint"` (V4).
#' @param polish_em Integer; pooled multinomial-EM polish iterations for the
#'   anchored \eqn{\hat\Phi}. Default 3.
#' @param readout_gn,refine_gn Integer GN iterations per read-out / per refinement
#'   sweep. Defaults 10 and 2.
#' @param max_sweeps,rule_tol,patience B-stationarity stop parameters
#'   (defaults 100, 1e-3, 2).
#' @param dz_cap Trust-region cap on the per-document GN step norm. Default 1.
#' @param V Optional K x (K-1) ILR contrast matrix.
#' @param verbose Logical.
#' @return An object of class `"sgscatm_chain"`: `Bz` (P x (K-1), oriented,
#'   generative-scale), `Z` (refined scores), `Phi` (anchored, K x N), `Theta`
#'   (M x K), `C_centred`, `W_tilde`, `V`, `K`, `pilot_Bz`, `pilot_Z`, the pilot
#'   `eigenvalues`, refinement diagnostics (`sweeps`, `rule_stop`, `monotone_ok`,
#'   `n_fail`), `anchor_Phi`, and `call`.
#' @seealso [vcov.sgscatm_chain()], [chain_se()], [sgscatm()]
#' @export
sgscatm_chain <- function(W, C, K, lambda = 1,
                          refine = c("joint", "frozen_phi"),
                          polish_em = 3L, readout_gn = 10L, refine_gn = 2L,
                          max_sweeps = 100L, rule_tol = 1e-3, patience = 2L,
                          dz_cap = 1, V = NULL, verbose = FALSE) {
  cl <- match.call()
  refine <- match.arg(refine)
  W <- as.matrix(W); C <- as.matrix(C)
  stopifnot(is.numeric(W), all(W >= 0), is.numeric(C), K >= 2L)
  K <- as.integer(K)
  M <- nrow(W); N <- ncol(W); P <- ncol(C)
  if (nrow(C) != M) stop("W and C must have the same number of rows.")

  if (is.null(V)) V <- ilr_contrast(K)
  Wf <- W / pmax(rowSums(W), 1)                     # row frequencies

  # --- Stage 1: spectral pilot (unit-norm subspace) ---
  if (verbose) message("Stage 1: spectral pilot")
  fit <- sgscatm(W, C, K = K, lambda = lambda, V = V, scale_W = TRUE,
                 rotate = FALSE)
  Cc <- fit$C_centred                               # centred covariates

  # --- Stage 2: anchored orientation ---
  if (verbose) message("Stage 2: anchors + read-out + general-linear alignment")
  ap  <- .sg_anchor_pipeline(W, K)
  Phi <- ap$Phi
  if (polish_em > 0L) {
    Th0 <- .sg_theta(.sg_readout_gn(Phi, Wf, V, n_gn = readout_gn,
                                    dz_cap = dz_cap)$Z, V)
    Phi <- .sg_phi_em_polish(Phi, W, Th0, n_it = polish_em)
  }
  ro  <- .sg_readout_gn(Phi, Wf, V, n_gn = readout_gn, dz_cap = dz_cap)
  gl  <- .sg_gl_align(fit$Z, ro$Z)                  # general-linear (not Procrustes)
  Z0  <- gl$Z

  # --- Stage 3: k-step LS refinement (frozen Phi or joint) ---
  if (verbose) message(sprintf("Stage 3: k-step refinement (%s)", refine))
  Phi_start <- if (refine == "frozen_phi") Phi else .sg_phi_step(Z0, Wf, V)
  rf <- .sg_refine(Z0, Phi_start, Wf, Cc, V, mode = refine,
                   max_sweeps = max_sweeps, rule_tol = rule_tol,
                   patience = patience, n_gn = refine_gn, dz_cap = dz_cap)

  Bz <- .sg_b_step(rf$Z, Cc)
  dimnames(Bz) <- list(colnames(C), NULL)
  Theta <- .sg_theta(rf$Z, V)

  structure(list(
    Bz = Bz, Z = rf$Z, Phi = rf$Phi, Theta = Theta,
    C_centred = Cc, W_tilde = fit$W_tilde, Wf = Wf, V = V, K = K, N = N, M = M,
    pilot_Bz = fit$Bz, pilot_Z = fit$Z, eigenvalues = fit$eigenvalues,
    anchor_Phi = ap$Phi, readout_Z = ro$Z, gl_A = gl$A,
    refine = refine, sweeps = rf$sweeps, rule_stop = rf$rule_stop,
    monotone_ok = rf$monotone_ok, n_fail = rf$n_fail, lambda = lambda,
    call = cl
  ), class = "sgscatm_chain")
}

#' @export
print.sgscatm_chain <- function(x, ...) {
  cat(sprintf("sgscatm_chain fit: M=%d docs, N=%d terms, K=%d topics, P=%d covariates\n",
              x$M, x$N, x$K, nrow(x$Bz)))
  cat(sprintf("  refinement: %s | sweeps=%d (rule stop=%s) | monotone=%s | GN failures=%d\n",
              x$refine, x$sweeps,
              ifelse(is.na(x$rule_stop), "cap", x$rule_stop),
              x$monotone_ok, x$n_fail))
  cat("  Bz (path coefficients, generative scale):\n")
  print(round(x$Bz, 4))
  invisible(x)
}

#' @export
coef.sgscatm_chain <- function(object, ...) object$Bz
