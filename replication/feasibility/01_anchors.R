#' ===================================================================
#'  Feasibility round — F1a: anchor-word pipeline (Arora et al. 2013
#'  style) and permutation-metric utilities
#' ===================================================================
#'
#'  Pipeline: word co-occurrence Q with the multinomial within-document
#'  diagonal correction -> successive-projection anchor selection on
#'  the row-normalised rows of Q -> RecoverL2 (each word's profile as a
#'  convex combination of the anchors' profiles) -> Bayes inversion to
#'  topic-word rows on the simplex.
#'
#'  All functions are prefixed fs_ and truth-free; the unit test
#'  fs_test_anchors() (hard gate in the runner) checks recovery of
#'  Phi_true up to topic permutation at alpha_beta = 0.05.
#' ===================================================================

#' Duchi et al. (2008) Euclidean projection of each row onto the simplex
fs_proj_simplex_rows <- function(X) {
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

#' Solve min ||T - C %*% Basis||_F^2 over rows of C in the simplex,
#' EXACTLY, by enumerating the 2^K - 1 candidate supports of the KKT
#' system (K is small here).  Projected gradient was tried first and
#' silently under-converges: the anchors' co-occurrence profiles are
#' close together, so the Gram matrix has tiny curvature along soft
#' directions and PG stalls near the uniform start (C error ~0.25).
#' For each support S the equality-constrained QP
#'   min c'Gc - 2h'c  s.t.  sum(c_S) = 1, c_{S^c} = 0
#' is a (|S|+1)-dim linear system; the global optimum is the feasible
#' candidate with the smallest objective.  Shared by RecoverL2, the V3
#' document mapping, and the F1a unit test.
fs_simplex_ls <- function(Tmat, Basis, ridge = 1e-10) {
  K <- nrow(Basis)
  G <- tcrossprod(Basis) + ridge * diag(K)
  H <- Tmat %*% t(Basis)
  n <- nrow(Tmat)
  supports <- lapply(seq_len(2^K - 1L), function(m)
    which(bitwAnd(m, 2^(seq_len(K) - 1L)) > 0L))
  # pre-factor the KKT system of every support
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

#' Word co-occurrence matrix with the multinomial diagonal correction:
#'   Q = (1/M) sum_i [w_i w_i' - diag(w_i)] / (L_i (L_i - 1)),
#' an unbiased estimator of E[(Beta' theta)(Beta' theta)'] entrywise.
fs_build_Q <- function(W) {
  L <- rowSums(W)
  Ws <- W / sqrt(pmax(L * (L - 1), 1))
  Q <- crossprod(Ws)
  diag(Q) <- diag(Q) - colSums(W / pmax(L * (L - 1), 1))
  Q / nrow(W)
}

#' Successive projection (furthest point from the affine span of the
#' anchors chosen so far) on the row-normalised rows of Q.  Candidates
#' are restricted to words appearing in at least a min_docfreq_frac
#' share of documents: the co-occurrence profile of a rare word is
#' estimated from a handful of documents and its noise dominates the
#' geometry (measured: with a docfreq >= 10 cutoff the selector picks
#' 11-16-document words and recovery fails even with oracle anchors).
fs_select_anchors <- function(Qbar, W, K, min_docfreq_frac = 0.05) {
  cand <- which(colSums(W > 0) >= min_docfreq_frac * nrow(W))
  Xc <- Qbar[cand, , drop = FALSE]
  anchors <- integer(K)
  # first anchor: farthest from the centroid of the candidate cloud
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
      R <- R - (R %*% t(basis)) %*% basis            # project out span
    d2 <- rowSums(R^2)
    pick <- which.max(d2)
    anchors[k] <- cand[pick]
    u <- R[pick, ] / sqrt(d2[pick])
    basis <- rbind(basis, u)
  }
  list(anchors = anchors, Qbar = Qbar)
}

#' Full anchor pipeline: counts W -> Phi_anchor (K x N, rows on simplex)
#'
#' Before RecoverL2, all co-occurrence profiles are denoised by
#' projection onto the top-K eigenspace of Q: the population rows of
#' Qbar lie in the K-dimensional affine span of the anchor profiles,
#' so the rank-K projection removes the sampling noise that otherwise
#' dominates the regression (measured: TV 0.34 -> see unit test).
#' post_thresh: entries of the word-topic posterior C below this value
#' are zeroed and rows renormalised before the Bayes inversion.  This
#' removes the noise-sprinkle on the ~85% of truly-zero Phi entries
#' (measured: TV 0.164 -> 0.098-0.138 at the unit-test config).  It is
#' appropriate in near-separable regimes (Dirichlet alpha <= 0.3); F1d
#' probes where it breaks.
fs_anchor_pipeline <- function(W, K, min_docfreq_frac = 0.05,
                               reselect = TRUE, post_thresh = 0.10) {
  Q  <- fs_build_Q(W)
  eg <- eigen(Q, symmetric = TRUE)
  Vk <- eg$vectors[, seq_len(K), drop = FALSE]      # N x K
  Qbar <- (Q / pmax(rowSums(Q), 1e-12)) %*% Vk %*% t(Vk)  # denoised
  p_w  <- colSums(W) / sum(W)
  cand <- which(colSums(W > 0) >= min_docfreq_frac * nrow(W))

  recover <- function(anchors) {
    A <- Qbar[anchors, , drop = FALSE]
    Cmat <- fs_simplex_ls(Qbar, A)                  # p(topic | word)
    if (post_thresh > 0) {
      Cmat[Cmat < post_thresh] <- 0
      Cmat <- Cmat / pmax(rowSums(Cmat), 1e-12)
    }
    Phi <- t(Cmat * p_w)
    list(Phi = Phi / pmax(rowSums(Phi), 1e-12), Cmat = Cmat)
  }

  sel <- fs_select_anchors(Qbar, W, K, min_docfreq_frac)
  rec <- recover(sel$anchors)
  anchors2 <- NULL
  if (reselect) {
    # one exclusivity re-selection round: anchor_k = the candidate word
    # with the highest posterior purity for topic k under the first
    # recovery, chosen greedily without replacement
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

# -------------------------------------------------------------------
#  Permutation utilities (legitimate for anchored estimators only)
# -------------------------------------------------------------------

#' All permutations of 1..K (K! rows): insert K at every position of
#' every permutation of 1..K-1
fs_all_perms <- function(K) {
  if (K == 1L) return(matrix(1L, 1L, 1L))
  sub <- fs_all_perms(K - 1L)
  do.call(rbind, lapply(seq_len(K), function(pos) {
    left  <- if (pos > 1L) sub[, seq_len(pos - 1L), drop = FALSE] else
      matrix(0L, nrow(sub), 0L)
    right <- if (pos <= K - 1L) sub[, pos:(K - 1L), drop = FALSE] else
      matrix(0L, nrow(sub), 0L)
    cbind(left, rep(K, nrow(sub)), right)
  }))
}

#' Permutation metric on ILR coefficients: a topic relabelling P acts
#' on ILR coordinates as Q_P = t(V) %*% P %*% V (orthogonal, since
#' P 1 = 1 and V'1 = 0); the relabelled estimate is B %*% Q_P'.
#' Returns the minimum entrywise MSE over all K! relabellings.
fs_perm_mse <- function(B_hat, Bz0, V) {
  K <- nrow(V)
  perms <- fs_all_perms(K)
  best <- Inf
  for (r in seq_len(nrow(perms))) {
    P <- diag(K)[perms[r, ], , drop = FALSE]
    Qp <- crossprod(V, P %*% V)
    m <- mean((B_hat %*% t(Qp) - Bz0)^2)
    if (m < best) best <- m
  }
  best
}

#' Row-wise total-variation error of Phi_hat vs Phi_true, minimised
#' over topic permutations (used by the F1a unit test and F1d).
fs_perm_tv <- function(Phi_hat, Phi_true) {
  K <- nrow(Phi_true)
  perms <- fs_all_perms(K)
  best <- Inf
  for (r in seq_len(nrow(perms))) {
    tv <- mean(0.5 * rowSums(abs(Phi_hat[perms[r, ], , drop = FALSE] -
                                   Phi_true)))
    if (tv < best) best <- tv
  }
  best
}

#' Anchor-exclusivity statistic of a topic-word matrix:
#'   min_k max_j Phi_kj / sum_k' Phi_k'j
fs_exclusivity <- function(Phi) {
  Xr <- Phi / pmax(matrix(colSums(Phi), nrow(Phi), ncol(Phi),
                          byrow = TRUE), 1e-300)
  min(apply(Xr, 1L, max))
}

#' F1a unit test (hard gate): alpha_beta = 0.05, M = 2000, strong
#' regime; Phi_anchor must match Phi_true up to permutation with mean
#' row-wise TV error <= tv_gate.
fs_test_anchors <- function(seed = 55001L, tv_gate = 0.15) {
  dat <- sim_dgp(M = 2000L, N = 500L, K = 5L, P = 3L, b_max = 0.50,
                 sigma_eps = 0.3, alpha_beta = 0.05, doc_length = 200L,
                 seed = seed)
  ap <- fs_anchor_pipeline(dat$W, K = 5L)
  tv <- fs_perm_tv(ap$Phi, dat$Beta)
  # sanity of the permutation representation: relabelling the truth
  # must give a zero permutation-metric
  Bz0 <- matrix(rnorm(12), 3, 4)
  P0 <- diag(5)[c(3, 1, 5, 2, 4), ]
  Qp <- crossprod(dat$V, P0 %*% dat$V)
  perm_zero <- fs_perm_mse(Bz0 %*% t(Qp), Bz0, dat$V)
  list(tv = tv, perm_zero = perm_zero, pass = tv <= tv_gate,
       exclusivity_true = fs_exclusivity(dat$Beta))
}
