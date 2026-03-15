#' Fit ILR-EGSCA Structural Topic Model
#'
#' Estimates document-topic proportions, topic-term loadings, and covariate
#' path coefficients by solving a covariate-augmented eigenvalue problem in
#' ILR coordinates.
#'
#' @param W Numeric matrix M x N. Document-term matrix (non-negative).
#' @param C Numeric matrix M x P. Covariate matrix. Will be column-centred
#'   internally if not already centred.
#' @param K Integer >= 2. Number of topics.
#' @param lambda Numeric >= 0. Regularisation weight balancing reconstruction
#'   fidelity and covariate structure. Default 1.
#' @param r Integer. Rank of the truncated SVD of the centred document-term
#'   matrix. Defaults to min(M, N, 100).
#' @param V Numeric K x (K-1) ILR contrast matrix. If NULL, computed via
#'   [ilr_contrast()].
#' @param scale_W Logical. If TRUE (default), row-normalise W to term
#'   frequencies before fitting.
#' @param rotate Logical. If TRUE (default), apply varimax rotation to the ILR
#'   loading matrix `Psi` after estimation (step 8 of the algorithm). Rotation
#'   is cost-free (Proposition on rotation invariance): it does not change the
#'   value of the profiled objective but improves topic interpretability by
#'   inducing simple structure in the factor pattern. Ignored when K = 2
#'   (trivial single direction).
#'
#' @return An object of class `"egscatm"` with components:
#'   \describe{
#'     \item{Pi}{M x K matrix of document-topic proportions (rows sum to 1).}
#'     \item{Phi}{K x N matrix of topic-term loadings.}
#'     \item{Bz}{P x (K-1) matrix of ILR path coefficients.}
#'     \item{Z}{M x (K-1) matrix of ILR topic scores.}
#'     \item{Psi}{(K-1) x N matrix of ILR topic-term deviations.}
#'     \item{eigenvalues}{Top K-1 eigenvalues of the augmented similarity matrix.}
#'     \item{V}{The K x (K-1) ILR contrast matrix used.}
#'     \item{K}{Number of topics.}
#'     \item{lambda}{Regularisation parameter used.}
#'     \item{rotate}{Logical. Whether varimax rotation was applied.}
#'     \item{R_star}{(K-1) x (K-1) varimax rotation matrix, or NULL if
#'       `rotate = FALSE`.}
#'     \item{call}{The matched call.}
#'   }
#' @export
egscatm <- function(W, C, K, lambda = 1, r = NULL, V = NULL, scale_W = TRUE,
                    rotate = TRUE) {
  cl <- match.call()

  # --- input checks ---
  W <- as.matrix(W)
  C <- as.matrix(C)
  stopifnot(is.numeric(W), all(W >= 0))
  stopifnot(is.numeric(C))
  stopifnot(is.numeric(K), length(K) == 1L, K >= 2L)
  stopifnot(is.numeric(lambda), length(lambda) == 1L, lambda >= 0)
  K <- as.integer(K)

  M <- nrow(W); N <- ncol(W); P <- ncol(C)
  if (nrow(C) != M) stop("W and C must have the same number of rows.")
  if (K > min(M, N)) stop("K must be <= min(nrow(W), ncol(W)).")

  # --- optional row-normalisation ---
  if (scale_W) {
    rs <- rowSums(W)
    rs[rs == 0] <- 1
    W <- W / rs
  }

  # --- column-centre C ---
  C <- scale(C, center = TRUE, scale = FALSE)

  # --- ILR contrast matrix ---
  if (is.null(V)) V <- ilr_contrast(K)
  stopifnot(nrow(V) == K, ncol(V) == K - 1L)

  # --- Step 1: centre W ---
  w_bar  <- colMeans(W)                        # length N
  W_tilde <- sweep(W, 2L, w_bar, "-")         # M x N

  # --- Step 2: truncated SVD ---
  if (is.null(r)) r <- min(M, N, 100L)
  r <- min(r, M - 1L, N, K - 1L + ncol(C))    # safety cap
  svd_res <- .trunc_svd(W_tilde, r)
  U_r  <- svd_res$u                            # M x r
  sig_r <- svd_res$d                           # length r

  # --- Step 3: QR of C ---
  qr_C <- qr(C)
  Q_C  <- qr.Q(qr_C)                          # M x P

  # --- Step 4: augmented factor H ---
  H <- cbind(U_r * rep(sig_r, each = M),       # M x r  (U_r %*% diag(sig_r))
             sqrt(lambda) * Q_C)               # M x P

  # --- Step 5: eigendecompose H'H ---
  HtH <- crossprod(H)                          # (r+P) x (r+P)
  eig  <- eigen(HtH, symmetric = TRUE)

  # top K-1 eigenvectors/values
  idx   <- seq_len(K - 1L)
  E_top <- eig$vectors[, idx, drop = FALSE]    # (r+P) x (K-1)
  s_top <- eig$values[idx]                     # K-1
  s_top <- pmax(s_top, 0)                      # numerical safety

  # all eigenvalues / eigenvectors of H'H (for SE computation)
  # U_H  = H E diag(s^{-1/2}): columns are eigenvectors of S_z = H H'
  n_eig  <- length(eig$values)
  s_all  <- pmax(eig$values, 0)
  # non-zero eigenvalues only
  nz     <- which(s_all > .Machine$double.eps * max(s_all) * 1e4)

  # --- Step 6: ILR topic scores Z* ---
  Z_star <- H %*% E_top %*% diag(1 / sqrt(s_top), K - 1L)  # M x (K-1)

  # --- Step 7: ILR topic-term deviations ---
  Psi_hat <- K * crossprod(Z_star, W_tilde)    # (K-1) x N

  # --- Step 8: varimax rotation (optional) ---
  # Rotation is cost-free: it leaves tr(Z'S_z Z) unchanged (Proposition on
  # rotation invariance) while inducing simple structure in Psi, which
  # improves topic exclusivity. Oblique rotations are not supported.
  R_star <- NULL
  if (isTRUE(rotate) && K > 2L) {
    vx     <- varimax(t(Psi_hat), normalize = FALSE)  # varimax on N x (K-1)
    R_star <- vx$rotmat                               # (K-1) x (K-1) orthogonal
    Z_star  <- Z_star  %*% R_star                     # M x (K-1)
    Psi_hat <- t(t(Psi_hat) %*% R_star)              # R*^T Psi: (K-1) x N
  }

  # --- Step 9: topic-term matrix ---
  Phi_hat <- V %*% Psi_hat +
    matrix(1, K, 1) %*% matrix(w_bar, 1, N)  # K x N

  # --- Step 10: path coefficients ---
  Bz_hat <- solve(crossprod(C), crossprod(C, Z_star))  # P x (K-1)

  # --- Step 11: topic proportions via softmax(V z_i*) ---
  scores <- Z_star %*% t(V)                   # M x K  (each row = V z_i*)
  Pi_hat <- .softmax_rows(scores)             # M x K

  # All M-dim eigenvectors of S_z stored for SE computation:
  # U_all[, j] = H %*% eig$vectors[, j] / sqrt(s_all[j])  for j in nz
  U_all <- H %*% eig$vectors[, nz, drop = FALSE] *
    rep(1 / sqrt(s_all[nz]), each = M)        # M x |nz|

  # --- assemble output ---
  structure(
    list(
      Pi          = Pi_hat,
      Phi         = Phi_hat,
      Bz          = Bz_hat,
      Z           = Z_star,
      Psi         = Psi_hat,
      eigenvalues = s_top,
      eigenvalues_all = s_all[nz],   # for SE: all non-zero eigenvalues of S_z
      U_all       = U_all,           # for SE: all non-zero eigenvectors of S_z
      W_tilde     = W_tilde,         # centred DTM (stored for SE)
      w_bar       = w_bar,
      C_centred   = C,
      V           = V,
      K           = K,
      lambda      = lambda,
      scale_W     = scale_W,
      rotate      = isTRUE(rotate) && K > 2L,
      R_star      = R_star,
      call        = cl
    ),
    class = "egscatm"
  )
}

# ---------- internal helpers ----------

.trunc_svd <- function(X, r, n_iter = 2L, oversampling = 10L) {
  # Truncated SVD keeping top-r left singular vectors and singular values.
  #
  # For small matrices (M*N < 1e6) or r close to full rank: uses exact LAPACK
  # svd() to preserve numerical exactness (needed for unit tests).
  #
  # For large matrices: uses the randomized algorithm of
  # Halko, Martinsson & Tropp (2011), "Finding structure with randomness",
  # SIAM Review 53(2):217-288.  With n_iter=2 power iterations and
  # oversampling=10 the approximation error is O(sigma_{r+1}) * machine_eps
  # in practice — negligible for topic modelling.
  r <- min(r, nrow(X) - 1L, ncol(X))

  if (nrow(X) * ncol(X) < 1e6 || r >= min(nrow(X), ncol(X)) / 2L) {
    # exact path (small matrices / near-full-rank)
    sv <- svd(X, nu = r, nv = 0L)
    return(list(u = sv$u[, seq_len(r), drop = FALSE],
                d = sv$d[seq_len(r)]))
  }

  # --- randomized path -------------------------------------------------
  # Stage A: build a rank-(r+oversampling) orthonormal basis Q for range(X)
  r_eff  <- min(r + oversampling, nrow(X) - 1L, ncol(X))
  Omega  <- matrix(rnorm(ncol(X) * r_eff), ncol(X), r_eff)  # random sketch
  Y      <- X %*% Omega                                       # M x r_eff
  # Power iterations amplify the top-r signal vs noise
  for (i in seq_len(n_iter)) {
    Y <- X %*% crossprod(X, Y)                               # (X X')^i * Y
  }
  Q <- qr.Q(qr(Y))                                           # M x r_eff, orthonormal

  # Stage B: small SVD in the projected subspace
  B    <- crossprod(Q, X)                                     # r_eff x N
  sv_B <- svd(B, nu = r, nv = 0L)
  list(u = Q %*% sv_B$u[, seq_len(r), drop = FALSE],
       d = sv_B$d[seq_len(r)])
}

.softmax_rows <- function(X) {
  # numerically stable row-wise softmax
  X <- X - apply(X, 1L, max)
  E <- exp(X)
  E / rowSums(E)
}
