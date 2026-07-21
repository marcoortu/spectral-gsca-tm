#' Anchor-word recovery of the topic-term matrix (internal)
#'
#' Promoted verbatim (up to renaming) from the validated
#' `replication/feasibility/01_anchors.R`. Implements Section 3.6(i): the
#' factorial word co-occurrence matrix with the multinomial diagonal
#' correction, successive-projection (Gillis-Vavasis / Arora et al.) anchor
#' selection on the row-normalised co-occurrence rows, constrained
#' least-squares recovery of the word-topic posterior on the simplex, and Bayes
#' inversion to a row-stochastic topic-term matrix. The unit test
#' `.sg_test_anchors()` (gate in `tests/testthat/test-anchors.R`) checks recovery
#' of a known Phi up to topic permutation.
#' @keywords internal
#' @name sgscatm-anchors-internal
NULL

# Duchi et al. (2008) Euclidean projection of each row onto the simplex
.sg_proj_simplex_rows <- function(X) {
  for (i in seq_len(nrow(X))) {
    v <- X[i, ]
    u <- sort(v, decreasing = TRUE)
    cs <- cumsum(u)
    rho <- max(which(u - (cs - 1) / seq_along(u) > 0))
    tau <- (cs[rho] - 1) / rho
    X[i, ] <- pmax(v - tau, 0)
  }
  X
}

# Exact simplex-constrained least squares by KKT enumeration over supports
# (K small). Solves min ||T - C Basis||_F^2 with rows of C on the simplex.
.sg_simplex_ls <- function(Tmat, Basis, ridge = 1e-10) {
  K <- nrow(Basis)
  G <- tcrossprod(Basis) + ridge * diag(K)
  H <- Tmat %*% t(Basis)
  n <- nrow(Tmat)
  supports <- lapply(seq_len(2^K - 1L), function(m)
    which(bitwAnd(m, 2^(seq_len(K) - 1L)) > 0L))
  kkt <- lapply(supports, function(S) {
    k <- length(S)
    A <- rbind(cbind(2 * G[S, S, drop = FALSE], 1), c(rep(1, k), 0))
    tryCatch(solve(A), error = function(e) NULL)
  })
  Cm <- matrix(0, n, K)
  for (i in seq_len(n)) {
    h <- H[i, ]
    best_f <- Inf; best_c <- rep(1 / K, K)
    for (s in seq_along(supports)) {
      if (is.null(kkt[[s]])) next
      S <- supports[[s]]
      sol <- kkt[[s]] %*% c(2 * h[S], 1)
      cS <- sol[seq_along(S)]
      if (any(cS < -1e-12)) next
      cS <- pmax(cS, 0); cS <- cS / sum(cS)
      cfull <- rep(0, K); cfull[S] <- cS
      f <- sum(cfull * (G %*% cfull)) - 2 * sum(cfull * h)
      if (f < best_f) { best_f <- f; best_c <- cfull }
    }
    Cm[i, ] <- best_c
  }
  Cm
}

# Word co-occurrence with the multinomial diagonal correction: unbiased for
# E[(Phi' theta)(Phi' theta)'] entrywise.
.sg_build_Q <- function(W) {
  L <- rowSums(W)
  Ws <- W / sqrt(pmax(L * (L - 1), 1))
  Q <- crossprod(Ws)
  diag(Q) <- diag(Q) - colSums(W / pmax(L * (L - 1), 1))
  Q / nrow(W)
}

# Successive projection on the row-normalised rows of Q (candidates restricted
# to words with document frequency >= min_docfreq_frac).
.sg_select_anchors <- function(Qbar, W, K, min_docfreq_frac = 0.05) {
  cand <- which(colSums(W > 0) >= min_docfreq_frac * nrow(W))
  Xc <- Qbar[cand, , drop = FALSE]
  anchors <- integer(K)
  ctr <- colMeans(Xc)
  d2 <- rowSums((Xc - matrix(ctr, nrow(Xc), ncol(Xc), byrow = TRUE))^2)
  pick <- which.max(d2)
  anchors[1L] <- cand[pick]
  p1 <- Xc[pick, ]
  Xd <- Xc - matrix(p1, nrow(Xc), ncol(Xc), byrow = TRUE)
  basis <- NULL
  for (k in 2L:K) {
    R <- Xd
    if (!is.null(basis))
      R <- R - (R %*% t(basis)) %*% basis
    d2 <- rowSums(R^2)
    pick <- which.max(d2)
    anchors[k] <- cand[pick]
    u <- R[pick, ] / sqrt(d2[pick])
    basis <- rbind(basis, u)
  }
  list(anchors = anchors, Qbar = Qbar)
}

#' Anchor-word recovery of a row-stochastic topic-term matrix
#'
#' Recovers \eqn{\hat{\boldsymbol\Phi}} (K x N, rows on the simplex) from a
#' non-negative document-term count matrix via the anchor-word pipeline of
#' Section 3.6(i). Truth-free.
#'
#' @param W Numeric M x N non-negative document-term count matrix.
#' @param K Integer number of topics.
#' @param min_docfreq_frac Minimum document-frequency share for anchor
#'   candidates. Default 0.05.
#' @param reselect Logical; run one exclusivity re-selection round. Default TRUE.
#' @param post_thresh Numeric; zero word-topic posterior entries below this and
#'   renormalise (near-separable denoising). Default 0.10.
#' @return A list with `Phi` (K x N row-stochastic), `anchors`, `anchors2`, `Q`.
#' @keywords internal
.sg_anchor_pipeline <- function(W, K, min_docfreq_frac = 0.05,
                                reselect = TRUE, post_thresh = 0.10) {
  Q  <- .sg_build_Q(W)
  eg <- eigen(Q, symmetric = TRUE)
  Vk <- eg$vectors[, seq_len(K), drop = FALSE]
  Qbar <- (Q / pmax(rowSums(Q), 1e-12)) %*% Vk %*% t(Vk)
  p_w  <- colSums(W) / sum(W)
  cand <- which(colSums(W > 0) >= min_docfreq_frac * nrow(W))

  recover <- function(anchors) {
    A <- Qbar[anchors, , drop = FALSE]
    Cmat <- .sg_simplex_ls(Qbar, A)
    if (post_thresh > 0) {
      Cmat[Cmat < post_thresh] <- 0
      Cmat <- Cmat / pmax(rowSums(Cmat), 1e-12)
    }
    Phi <- t(Cmat * p_w)
    list(Phi = Phi / pmax(rowSums(Phi), 1e-12), Cmat = Cmat)
  }

  sel <- .sg_select_anchors(Qbar, W, K, min_docfreq_frac)
  rec <- recover(sel$anchors)
  anchors2 <- NULL
  if (reselect) {
    Cc <- rec$Cmat[cand, , drop = FALSE]
    anchors2 <- integer(K)
    taken <- rep(FALSE, length(cand))
    for (k in order(-apply(Cc, 2L, max))) {
      sc <- Cc[, k]; sc[taken] <- -Inf
      pick <- which.max(sc)
      anchors2[k] <- cand[pick]
      taken[pick] <- TRUE
    }
    rec <- recover(anchors2)
  }
  list(Phi = rec$Phi, anchors = sel$anchors, anchors2 = anchors2, Q = Q)
}

# Pooled multinomial-EM polish of Phi at fixed Theta (the single deviance step;
# Phi-only, immune to the incidental-parameter mechanism). n_it E/M sweeps.
.sg_phi_em_polish <- function(Phi, W, Theta, n_it = 3L, eps = 1e-8) {
  Wc <- W                                   # counts
  for (it in seq_len(n_it)) {
    # E-step token responsibilities are implicit in the pooled M-step:
    #   phi_kj propto sum_i theta_ik w_ij / (Theta phi)_ij
    denom <- Theta %*% Phi                   # M x N expected freq
    denom[denom < eps] <- eps
    ratio <- Wc / denom                      # M x N
    Phi_new <- Phi * crossprod(Theta, ratio) # K x N (unnormalised)
    Phi_new <- Phi_new + eps
    Phi <- Phi_new / rowSums(Phi_new)
  }
  Phi
}
