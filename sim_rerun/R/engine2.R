# ===================================================================
#  engine2.R  —  clean-DGP estimator + packaged SE + jackknife (v2)
# ===================================================================
#  Reuses the package spectral solver sgscatm() and the NEW package
#  function sgscatm_vcov() (corrected 3-term influence SE).
# ===================================================================

# --- Build the clean well-separated DGP coefficient B0 -------------
build_clean_B0 <- function(K = K_TOPICS, P = P_COV, d = SCORE_D,
                           se2 = SIGMA_EPS2, seed = 20260704L) {
  stopifnot(P == K - 1L, length(d) == K - 1L, se2 < min(d))
  set.seed(seed)
  R0 <- qr.Q(qr(matrix(rnorm(P * P), P, P)))   # orthonormal 4x4 frame
  B0 <- R0 %*% diag(sqrt(d - se2))             # P x (K-1)
  eigcov <- eigen(crossprod(B0) + se2 * diag(P), only.values = TRUE)$values
  gaps <- abs(diff(sort(eigcov, decreasing = TRUE)))
  stopifnot(all(gaps > 0.1))                    # assert well-separated
  attr(B0, "R0") <- R0; attr(B0, "d") <- d
  attr(B0, "sigma_eps") <- sqrt(se2)
  attr(B0, "eig_cov_z") <- sort(eigcov, decreasing = TRUE)
  attr(B0, "min_gap") <- min(gaps)
  B0
}

gen_clean <- function(M, B0, seed) {
  sim_dgp(M = M, N = N_VOCAB, K = K_TOPICS, P = P_COV, Bz0 = B0,
          sigma_eps = attr(B0, "sigma_eps"), alpha_beta = ALPHA_BETA,
          doc_length = DOC_LEN, V = ilr_contrast(K_TOPICS), seed = seed)
}

# --- Standardized estimator + corrected SE from a fit --------------
# B_hat = sgscatm_vcov()$B (word-SVD standardized, lambda-inert); the SE
# is the packaged 3-term influence, rotationally identified via `rotation`.
corrected_fit2 <- function(dat, B0, lambda) {
  t0  <- proc.time()
  fit <- sgscatm(dat$W, dat$C, K = K_TOPICS, lambda = lambda, rotate = TRUE)
  t_fit <- (proc.time() - t0)[3L]
  v0  <- sgscatm_vcov(fit)                       # B_raw (standardized)
  list(fit = fit, B_raw = v0$B, time_fit = t_fit)
}

# Align to B_star and return aligned estimate + rotated SE.
aligned_se2 <- function(cf, B_star) {
  sv <- svd(crossprod(cf$B_raw, B_star))
  R  <- sv$u %*% t(sv$v)
  vr <- sgscatm_vcov(cf$fit, rotation = R)
  list(B_al = vr$B, se = vr$se, R = R)
}

# --- Pilot pseudo-true estimand B_star -----------------------------
pilot_Bstar <- function(B0, lambda, M_pilot = M_PILOT, seed = 7L) {
  dat <- gen_clean(M_pilot, B0, seed)
  fit <- sgscatm(dat$W, dat$C, K = K_TOPICS, lambda = lambda, rotate = TRUE)
  sgscatm_vcov(fit)$B                             # standardized estimand
}

# --- Delete-a-block jackknife over documents (independent SE) ------
loo_block_jack <- function(dat, B_star, lambda, G = JK_BLOCKS) {
  M <- nrow(dat$W); P <- ncol(dat$C); Km1 <- K_TOPICS - 1L
  blk <- rep(seq_len(G), length.out = M)
  Bj  <- array(NA_real_, dim = c(G, P, Km1))
  for (g in seq_len(G)) {
    keep <- blk != g
    fg <- tryCatch(sgscatm(dat$W[keep, , drop = FALSE],
                           dat$C[keep, , drop = FALSE],
                           K = K_TOPICS, lambda = lambda, rotate = TRUE),
                   error = function(e) NULL)
    if (is.null(fg)) next
    Bg <- sgscatm_vcov(fg)$B
    sv <- svd(crossprod(Bg, B_star))
    Bj[g, , ] <- Bg %*% (sv$u %*% t(sv$v))
  }
  ok  <- apply(Bj, 1L, function(x) all(is.finite(x)))
  Bok <- Bj[ok, , , drop = FALSE]; n <- sum(ok)
  mb  <- apply(Bok, c(2L, 3L), mean)
  ss  <- apply(sweep(Bok, c(2L, 3L), mb, "-")^2, c(2L, 3L), sum)
  sqrt((n - 1) / n * ss)                          # P x (K-1) jackknife SE
}
