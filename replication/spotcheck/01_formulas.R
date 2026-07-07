#' ===================================================================
#'  Spot-check — Section 4 objects: bias field b(z) and gradient
#'  direction G0, plus the certified per-document GN solver
#' ===================================================================
#'
#'  All formulas were re-derived independently before implementation
#'  (expansion of the LS normal equations psi(z) = J(z)'(w - p(z)) = 0;
#'  the third-derivative term drops at O(1/L) because multinomial third
#'  moments are O(L^-2)); the transcription in the brief matches the
#'  re-derivation.  SC0 unit-tests the derivative building blocks
#'  against numDeriv and the optimized implementation against a naive
#'  one.  Functions are prefixed sc_.
#' ===================================================================

sc_softmax <- function(x) { x <- x - max(x); e <- exp(x); e / sum(e) }

#' All Section-4 objects at a single z (optimized: no N x N products —
#' Sigma_e = diag(p) - p p' is applied analytically; the sum_j H2f_j
#' contractions are collapsed over topics, K terms instead of N).
sc_objects <- function(z, Phi0, V, want_G0 = TRUE) {
  K <- nrow(Phi0); N <- ncol(Phi0); Km1 <- length(z)
  th <- sc_softmax(as.vector(V %*% z))
  Om <- diag(th) - tcrossprod(th)
  p  <- as.vector(crossprod(Phi0, th))
  J  <- crossprod(Phi0, Om %*% V)               # N x (K-1)
  H  <- crossprod(J)
  Hi <- solve(H)
  Jp <- as.vector(crossprod(J, p))              # J'p
  # T0 = J' Sigma_e  (K-1) x N  = J'diag(p) - (J'p) p'
  T0 <- t(J * p) - tcrossprod(Jp, p)
  JSJ <- crossprod(J, J * p) - tcrossprod(Jp)   # J' Sigma_e J
  W  <- Hi %*% JSJ %*% Hi

  # Hessians of theta_k
  VOV <- crossprod(V, Om %*% V)
  H2th <- lapply(seq_len(K), function(k) {
    vk <- as.vector(crossprod(V, (seq_len(K) == k) - th))
    th[k] * (tcrossprod(vk) - VOV)
  })

  # b1 = Hi sum_j H2f_j T1[, j],  T1 = Hi T0; collapse over topics:
  # sum_j H2f_j x_j = sum_k H2th_k (X %*% Phi0[k, ])
  T1 <- Hi %*% T0
  s1 <- Reduce(`+`, lapply(seq_len(K), function(k)
    H2th[[k]] %*% (T1 %*% Phi0[k, ])))
  b1 <- Hi %*% s1

  # a2 = sum_j tr(H2f_j W) J[j, ] + 2 sum_j H2f_j W J[j, ]
  JtPhi <- lapply(seq_len(K), function(k) as.vector(crossprod(J, Phi0[k, ])))
  a2 <- Reduce(`+`, lapply(seq_len(K), function(k)
    sum(diag(H2th[[k]] %*% W)) * JtPhi[[k]] +
      2 * H2th[[k]] %*% W %*% JtPhi[[k]]))
  b <- b1 - 0.5 * Hi %*% a2

  G0 <- NULL
  if (want_G0) {
    OV <- Om %*% V
    G0 <- OV %*% Hi %*% T0 - tcrossprod(th, J %*% b) - OV %*% W %*% t(J)
  }
  list(theta = th, p = p, J = J, H = H, Hi = Hi, W = W,
       b1 = as.vector(b1), a2 = as.vector(a2), b = as.vector(b),
       G0 = G0, H2th = H2th)
}

#' Naive reference implementation (dense Sigma_e, explicit j-sums) —
#' used only by the SC0 unit test to validate the optimized version.
sc_objects_naive <- function(z, Phi0, V) {
  K <- nrow(Phi0); N <- ncol(Phi0)
  th <- sc_softmax(as.vector(V %*% z))
  Om <- diag(th) - tcrossprod(th)
  p  <- as.vector(crossprod(Phi0, th))
  Se <- diag(p) - tcrossprod(p)
  J  <- crossprod(Phi0, Om %*% V)
  H  <- crossprod(J); Hi <- solve(H)
  VOV <- crossprod(V, Om %*% V)
  H2th <- lapply(seq_len(K), function(k) {
    vk <- as.vector(crossprod(V, (seq_len(K) == k) - th))
    th[k] * (tcrossprod(vk) - VOV)
  })
  H2f <- lapply(seq_len(N), function(j)
    Reduce(`+`, lapply(seq_len(K), function(k) Phi0[k, j] * H2th[[k]])))
  b1 <- Hi %*% Reduce(`+`, lapply(seq_len(N), function(j)
    H2f[[j]] %*% (Hi %*% crossprod(J, Se[, j]))))
  W <- Hi %*% crossprod(J, Se %*% J) %*% Hi
  a2 <- Reduce(`+`, lapply(seq_len(N), function(j)
    sum(diag(H2f[[j]] %*% W)) * J[j, ] + 2 * H2f[[j]] %*% W %*% J[j, ]))
  b <- b1 - 0.5 * Hi %*% a2
  G0 <- Om %*% V %*% Hi %*% crossprod(J, Se) -
    tcrossprod(th, J %*% b) - Om %*% V %*% W %*% t(J)
  list(b1 = as.vector(b1), a2 = as.vector(a2), b = as.vector(b), G0 = G0)
}

# -------------------------------------------------------------------
#  Certified batch GN solver (reuses bc_z_step; all replicate
#  documents for one (z, L) cell are solved as one batch)
# -------------------------------------------------------------------

#' Per-row LS gradient 2 J_i'(p_i - w_i), vectorised
sc_grad_batch <- function(Z, Phi, Wf, V) {
  Th <- bc_theta(Z, V)
  R  <- Th %*% Phi - Wf
  U  <- R %*% t(Phi)
  TU <- Th * U
  2 * (TU - Th * rowSums(TU)) %*% V
}

#' Solve the per-document LS problem for every row of Wf, starting at
#' z_true, iterating bc_z_step until ||grad||_inf < tol per row.
#' Returns Z, the per-row certificate, and the number of failures.
sc_gn_batch <- function(Wf, Phi, V, z_true, tol = 1e-10,
                        max_calls = 60L) {
  Rn <- nrow(Wf)
  Z  <- matrix(z_true, Rn, length(z_true), byrow = TRUE)
  nu <- rep(1e-6, Rn)
  for (it in seq_len(max_calls)) {
    g <- sc_grad_batch(Z, Phi, Wf, V)
    gm <- apply(abs(g), 1L, max)
    if (max(gm) < tol) break
    zs <- bc_z_step(Z, Phi, Wf, V, lambda = 0, CB = NULL, nu = nu,
                    n_gn = 2L)
    Z <- zs$Z; nu <- zs$nu
  }
  gm <- apply(abs(sc_grad_batch(Z, Phi, Wf, V)), 1L, max)
  list(Z = Z, gmax = gm, ok = gm < tol, n_fail = sum(gm >= tol))
}

# -------------------------------------------------------------------
#  SC0 unit tests (hard gates; K = 3, N = 20)
# -------------------------------------------------------------------
sc_verify <- function(seed = 77000L) {
  set.seed(seed)
  K <- 3L; N <- 20L
  V <- ilr_contrast(K)
  Phi0 <- .rdirichlet_matrix(K, N, 0.3)
  z <- rnorm(K - 1L, 0, 0.5)

  # (i) H2theta_k vs numDeriv::hessian of theta_k(z)
  ob <- sc_objects(z, Phi0, V)
  e1 <- max(vapply(seq_len(K), function(k) {
    Hn <- numDeriv::hessian(function(zz)
      sc_softmax(as.vector(V %*% zz))[k], z)
    max(abs(ob$H2th[[k]] - Hn)) / max(abs(Hn))
  }, numeric(1)))
  if (e1 > 1e-8) stop(sprintf("SC0(i) H2theta FAILED: rel err %.2e", e1))

  # (ii) H2f_j vs numDeriv::hessian of f_j(z) = p_j(z)
  e2 <- max(vapply(c(1L, 7L, N), function(j) {
    Hn <- numDeriv::hessian(function(zz)
      sum(Phi0[, j] * sc_softmax(as.vector(V %*% zz))), z)
    Ha <- Reduce(`+`, lapply(seq_len(K), function(k)
      Phi0[k, j] * ob$H2th[[k]]))
    max(abs(Ha - Hn)) / max(abs(Hn))
  }, numeric(1)))
  if (e2 > 1e-6) stop(sprintf("SC0(ii) H2f FAILED: rel err %.2e", e2))

  # (ii') optimized vs naive implementation of b1, a2, b, G0
  nv <- sc_objects_naive(z, Phi0, V)
  e3 <- max(abs(c(ob$b1 - nv$b1, ob$a2 - nv$a2, ob$b - nv$b)),
            max(abs(ob$G0 - nv$G0))) /
    max(abs(c(nv$b, 1)))
  if (e3 > 1e-10)
    stop(sprintf("SC0(ii') optimized-vs-naive FAILED: %.2e", e3))

  # (iii) certified GN from the truth
  set.seed(seed + 1L)
  n <- t(rmultinom(200L, 100L, ob$p))
  gn <- sc_gn_batch(n / 100, Phi0, V, z, tol = 1e-10)
  if (gn$n_fail > 0L)
    stop(sprintf("SC0(iii) GN certificate FAILED: %d/200 rows", gn$n_fail))

  list(H2theta_relerr = e1, H2f_relerr = e2, opt_vs_naive = e3,
       gn_gmax = max(gn$gmax))
}
