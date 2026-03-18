#' ===================================================================
#' Data Generating Process for ILR-Spectral-GSCA Simulation Study
#' ===================================================================
#'
#' Generates synthetic corpora from a structural topic model consistent
#' with the spectral-gsca framework:
#'
#'   1. Covariates  C  ~ N(0, Sigma_C), then centred
#'   2. ILR scores  z_i = B_z0' c_i + eps_i,  eps ~ N(0, sigma_eps^2 I)
#'   3. Proportions  theta_i = closure(exp(V z_i))
#'   4. Word counts  w_i ~ Multinomial(L_i, theta_i %*% Beta)
#'
#' The function returns everything needed for estimation and evaluation:
#' the observed data (W, C), the ground truth (Bz0, Z_true, Theta_true,
#' Beta), and the contrast matrix V.
#'
# ===================================================================


#' Generate a synthetic corpus from the structural topic model
#'
#' @param M       Integer. Number of documents.
#' @param N       Integer. Vocabulary size.
#' @param K       Integer. Number of topics (>= 2).
#' @param P       Integer. Number of covariates.
#' @param Bz0     Numeric P x (K-1) matrix. True path coefficients.
#'                If NULL, generated with entries ~ Uniform(-b_max, b_max).
#' @param b_max   Numeric. Max abs value for random Bz0 entries. Default 0.5.
#' @param sigma_eps Numeric > 0. Std dev of ILR residual noise. Default 0.3.
#'                Controls how much topic proportions deviate from the
#'                structural prediction. Small = covariates explain most;
#'                large = substantial residual variation.
#' @param Sigma_C  Numeric P x P PD matrix. Covariate covariance.
#'                If NULL, uses identity I_P.
#' @param Beta    Numeric K x N matrix. Topic-word distributions (rows
#'                sum to 1). If NULL, generated from Dirichlet(alpha_beta).
#' @param alpha_beta Numeric > 0. Dirichlet concentration for Beta
#'                generation. Small = sparse/distinct topics;
#'                large = diffuse/similar topics. Default 0.1.
#' @param doc_length Integer or function. If integer, all documents have
#'                this many words. If function, called as f(M) to produce
#'                an M-vector of document lengths. Default 200.
#' @param V       Numeric K x (K-1) ILR contrast matrix. If NULL,
#'                Helmert contrast via ilr_contrast().
#' @param seed    Integer or NULL. Random seed.
#'
#' @return A list with:
#'   \describe{
#'     \item{W}{M x N document-term matrix (integer counts).}
#'     \item{C}{M x P column-centred covariate matrix.}
#'     \item{Bz0}{P x (K-1) true path coefficient matrix.}
#'     \item{Z_true}{M x (K-1) true ILR scores.}
#'     \item{Theta_true}{M x K true topic proportions.}
#'     \item{Beta}{K x N true topic-word distributions.}
#'     \item{V}{K x (K-1) ILR contrast matrix.}
#'     \item{doc_lengths}{M-vector of document lengths.}
#'     \item{params}{List of all generation parameters for reproducibility.}
#'   }
#' @export
sim_dgp <- function(M, N, K, P,
                    Bz0       = NULL,
                    b_max     = 0.5,
                    sigma_eps = 0.3,
                    Sigma_C   = NULL,
                    Beta      = NULL,
                    alpha_beta = 0.1,
                    doc_length = 200L,
                    V          = NULL,
                    seed       = NULL) {

  # --- reproducibility ---
  if (!is.null(seed)) set.seed(seed)

  # --- dimensions ---
  stopifnot(K >= 2L, P >= 1L, M >= K, N >= K)
  Km1 <- K - 1L

  # --- ILR contrast matrix ---
  if (is.null(V)) V <- ilr_contrast(K)

  # --- covariate covariance ---
  if (is.null(Sigma_C)) {
    Sigma_C <- diag(P)
  } else {
    stopifnot(nrow(Sigma_C) == P, ncol(Sigma_C) == P)
  }

  # --- true path coefficients ---
  if (is.null(Bz0)) {
    Bz0 <- matrix(runif(P * Km1, -b_max, b_max), P, Km1)
  } else {
    stopifnot(nrow(Bz0) == P, ncol(Bz0) == Km1)
  }

  # --- generate covariates ---
  # C_raw ~ N(0, Sigma_C), then column-centre
  if (P == 1L) {
    C_raw <- matrix(rnorm(M, 0, sqrt(Sigma_C[1, 1])), M, 1L)
  } else {
    # Cholesky factorisation for multivariate normal
    L_chol <- chol(Sigma_C)  # upper triangular: Sigma_C = L'L
    C_raw  <- matrix(rnorm(M * P), M, P) %*% L_chol
  }
  C <- scale(C_raw, center = TRUE, scale = FALSE)

  # --- generate ILR scores ---
  # z_i = Bz0' c_i + eps_i
  Eta   <- C %*% Bz0                                    # M x (K-1)
  Eps   <- matrix(rnorm(M * Km1, 0, sigma_eps), M, Km1) # M x (K-1)
  Z_true <- Eta + Eps                                    # M x (K-1)

  # --- topic proportions via ILR inverse ---
  scores     <- Z_true %*% t(V)                          # M x K
  scores_exp <- exp(scores - apply(scores, 1L, max))     # stable softmax
  Theta_true <- scores_exp / rowSums(scores_exp)          # M x K

  # --- topic-word distributions ---
  if (is.null(Beta)) {
    Beta <- .rdirichlet_matrix(K, N, alpha_beta)
  } else {
    stopifnot(nrow(Beta) == K, ncol(Beta) == N)
  }

  # --- document lengths ---
  if (is.function(doc_length)) {
    L <- doc_length(M)
  } else {
    L <- rep(as.integer(doc_length), M)
  }

  # --- generate word counts ---
  # For each document: sample from Multinomial(L_i, theta_i %*% Beta)
  word_probs <- Theta_true %*% Beta                       # M x N
  W <- matrix(0L, M, N)
  for (i in seq_len(M)) {
    W[i, ] <- rmultinom(1L, size = L[i], prob = word_probs[i, ])
  }

  # --- package output ---
  list(
    W          = W,
    C          = C,
    Bz0        = Bz0,
    Z_true     = Z_true,
    Theta_true = Theta_true,
    Beta       = Beta,
    V          = V,
    doc_lengths = L,
    params     = list(
      M = M, N = N, K = K, P = P,
      b_max = b_max, sigma_eps = sigma_eps,
      Sigma_C = Sigma_C, alpha_beta = alpha_beta,
      doc_length = doc_length, seed = seed
    )
  )
}


# ===================================================================
# Convenience wrappers for specific simulation scenarios
# ===================================================================

#' Generate DGP with correlated covariates
#'
#' Produces a covariate matrix with AR(1) correlation structure,
#' which is more realistic than independent covariates.
#'
#' @param M,N,K,P,... Passed to sim_dgp().
#' @param rho  Numeric in (-1,1). AR(1) correlation. Default 0.5.
#' @export
sim_dgp_correlated <- function(M, N, K, P, rho = 0.5, ...) {
  # AR(1) covariance: Sigma_C[i,j] = rho^|i-j|
  Sigma_C <- rho^abs(outer(seq_len(P), seq_len(P), "-"))
  sim_dgp(M = M, N = N, K = K, P = P, Sigma_C = Sigma_C, ...)
}


#' Generate DGP with varying document lengths
#'
#' Document lengths sampled from a negative binomial, mimicking
#' the heavy-tailed length distribution of real corpora.
#'
#' @param M,N,K,P,... Passed to sim_dgp().
#' @param mean_length Numeric. Mean document length. Default 200.
#' @param size_nb     Numeric. Size parameter of NegBin. Default 5.
#'                    Smaller = more overdispersion.
#' @export
sim_dgp_variable_length <- function(M, N, K, P,
                                     mean_length = 200, size_nb = 5,
                                     ...) {
  len_fun <- function(m) {
    pmax(rnbinom(m, size = size_nb, mu = mean_length), 10L)
  }
  sim_dgp(M = M, N = N, K = K, P = P, doc_length = len_fun, ...)
}


#' Generate DGP with controlled topic separation
#'
#' @param M,N,K,P,... Passed to sim_dgp().
#' @param separation Character. "high" (alpha=0.01), "medium" (0.1),
#'                   or "low" (1.0). Controls Dirichlet concentration
#'                   for Beta.
#' @export
sim_dgp_separation <- function(M, N, K, P,
                                separation = c("high", "medium", "low"),
                                ...) {
  separation <- match.arg(separation)
  alpha <- switch(separation,
                  high   = 0.01,
                  medium = 0.1,
                  low    = 1.0)
  sim_dgp(M = M, N = N, K = K, P = P, alpha_beta = alpha, ...)
}


# ===================================================================
# Internal helpers
# ===================================================================

#' Sample rows from a Dirichlet distribution
#' @param K Number of rows (topics).
#' @param N Number of columns (terms).
#' @param alpha Concentration parameter (scalar or N-vector).
#' @return K x N matrix, each row summing to 1.
#' @keywords internal
.rdirichlet_matrix <- function(K, N, alpha) {
  if (length(alpha) == 1L) alpha <- rep(alpha, N)
  mat <- matrix(rgamma(K * N, shape = rep(alpha, each = K)), K, N)
  mat / rowSums(mat)
}


# ===================================================================
# Validation: quick check that the DGP is self-consistent
# ===================================================================

#' Validate a DGP output
#'
#' Runs basic sanity checks on a sim_dgp() output.
#' Useful for debugging and unit tests.
#'
#' @param dat Output of sim_dgp().
#' @return Invisible TRUE if all checks pass.
#' @export
validate_dgp <- function(dat) {
  M <- dat$params$M
  N <- dat$params$N
  K <- dat$params$K
  P <- dat$params$P

  # dimensions
  stopifnot(nrow(dat$W) == M, ncol(dat$W) == N)
  stopifnot(nrow(dat$C) == M, ncol(dat$C) == P)
  stopifnot(nrow(dat$Bz0) == P, ncol(dat$Bz0) == K - 1L)
  stopifnot(nrow(dat$Z_true) == M, ncol(dat$Z_true) == K - 1L)
  stopifnot(nrow(dat$Theta_true) == M, ncol(dat$Theta_true) == K)
  stopifnot(nrow(dat$Beta) == K, ncol(dat$Beta) == N)
  stopifnot(nrow(dat$V) == K, ncol(dat$V) == K - 1L)

  # simplex constraints
  stopifnot(all(dat$Theta_true > 0))
  stopifnot(all(abs(rowSums(dat$Theta_true) - 1) < 1e-10))
  stopifnot(all(dat$Beta > 0))
  stopifnot(all(abs(rowSums(dat$Beta) - 1) < 1e-10))

  # non-negative counts
  stopifnot(all(dat$W >= 0L))

  # C is centred
  stopifnot(all(abs(colMeans(dat$C)) < 1e-10))

  # V properties
  VtV <- crossprod(dat$V)
  stopifnot(max(abs(VtV - diag(K - 1L))) < 1e-10)
  stopifnot(max(abs(crossprod(dat$V, rep(1, K)))) < 1e-10)

  message("All DGP checks passed.")
  invisible(TRUE)
}


# ===================================================================
# Quick demo (not run by default)
# ===================================================================

if (FALSE) {
  # Source the sgscatm package first
  # devtools::load_all("path/to/sgscatm")

  # Generate a small corpus
  dat <- sim_dgp(M = 1000, N = 500, K = 5, P = 3,
                 sigma_eps = 0.3, alpha_beta = 0.1,
                 doc_length = 200, seed = 42)
  validate_dgp(dat)

  # Fit sgscatm
  fit <- sgscatm(dat$W, dat$C, K = 5, lambda = 1)

  # Compare estimated vs true Bz
  cat("True Bz:\n"); print(round(dat$Bz0, 3))
  cat("Estimated Bz:\n"); print(round(fit$Bz, 3))

  # Note: Bz is identified only up to rotation (Theorem 3).
  # For comparison, align via Procrustes:
  #   R_align = svd(t(fit$Bz) %*% dat$Bz0)
  #   Bz_aligned = fit$Bz %*% R_align$u %*% t(R_align$v)
}
