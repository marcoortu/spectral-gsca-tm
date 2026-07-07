test_that("sgscatm_vcov returns a well-formed standardized-scale result", {
  set.seed(1)
  K <- 4L; P <- 3L; N <- 120L; M <- 400L
  Bz0 <- matrix(runif(P * (K - 1L), -0.4, 0.4), P, K - 1L)
  C  <- scale(matrix(rnorm(M * P), M, P), TRUE, FALSE)
  Z  <- C %*% Bz0 + matrix(rnorm(M * (K - 1L), 0, 0.4), M, K - 1L)
  V  <- ilr_contrast(K)
  Th <- ilr_to_proportions(Z, V)
  Beta <- matrix(rgamma(K * N, 0.2), K, N); Beta <- Beta / rowSums(Beta)
  WP <- Th %*% Beta
  W  <- t(vapply(seq_len(M), function(i) rmultinom(1, 200, WP[i, ]), numeric(N)))

  fit <- sgscatm(W, C, K = K, lambda = 0.1, rotate = TRUE)
  vc  <- sgscatm_vcov(fit)

  expect_type(vc, "list")
  expect_equal(dim(vc$se), c(P, K - 1L))
  expect_true(all(is.finite(vc$se)) && all(vc$se > 0))
  expect_equal(length(vc$rho), K - 1L)
  expect_identical(vc$scale, "standardized")
  # vcov is symmetric PSD-ish and square P(K-1)
  expect_equal(dim(vc$vcov), c(P * (K - 1L), P * (K - 1L)))
  expect_lt(max(abs(vc$vcov - t(vc$vcov))), 1e-8)
})

test_that("identity rotation leaves SEs unchanged", {
  set.seed(2)
  K <- 4L; P <- 3L; N <- 100L; M <- 350L
  C <- scale(matrix(rnorm(M * P), M, P), TRUE, FALSE)
  Z <- C %*% matrix(runif(P * (K - 1L), -.3, .3), P, K - 1L) +
    matrix(rnorm(M * (K - 1L), 0, .4), M, K - 1L)
  V <- ilr_contrast(K); Th <- ilr_to_proportions(Z, V)
  Beta <- matrix(rgamma(K * N, .2), K, N); Beta <- Beta / rowSums(Beta)
  WP <- Th %*% Beta
  W <- t(vapply(seq_len(M), function(i) rmultinom(1, 200, WP[i, ]), numeric(N)))
  fit <- sgscatm(W, C, K = K, lambda = .1)
  s0 <- sgscatm_vcov(fit)$se
  s1 <- sgscatm_vcov(fit, rotation = diag(K - 1L))$se
  expect_lt(max(abs(s0 - s1)), 1e-8)
})

test_that("standard errors shrink like M^{-1/2}", {
  skip_on_cran()
  set.seed(3)
  K <- 4L; P <- 3L; N <- 150L
  Bz0 <- matrix(runif(P * (K - 1L), -.4, .4), P, K - 1L)
  Beta <- matrix(rgamma(K * N, .2), K, N); Beta <- Beta / rowSums(Beta)
  V <- ilr_contrast(K)
  mean_se <- function(M) {
    C <- scale(matrix(rnorm(M * P), M, P), TRUE, FALSE)
    Z <- C %*% Bz0 + matrix(rnorm(M * (K - 1L), 0, .4), M, K - 1L)
    Th <- ilr_to_proportions(Z, V); WP <- Th %*% Beta
    W <- t(vapply(seq_len(M), function(i) rmultinom(1, 200, WP[i, ]), numeric(N)))
    mean(sgscatm_vcov(sgscatm(W, C, K = K, lambda = .1))$se)
  }
  s_small <- mean_se(500L); s_big <- mean_se(2000L)
  # 4x M should roughly halve the SE; allow generous tolerance
  expect_lt(s_big, s_small)
  expect_gt(s_big / s_small, 0.3)
  expect_lt(s_big / s_small, 0.85)
})

test_that("sgscatm(se = TRUE) attaches B_z_se without breaking defaults", {
  set.seed(4)
  K <- 4L; P <- 3L; N <- 100L; M <- 300L
  C <- scale(matrix(rnorm(M * P), M, P), TRUE, FALSE)
  Z <- C %*% matrix(runif(P * (K - 1L), -.3, .3), P, K - 1L) +
    matrix(rnorm(M * (K - 1L), 0, .4), M, K - 1L)
  V <- ilr_contrast(K); Th <- ilr_to_proportions(Z, V)
  Beta <- matrix(rgamma(K * N, .2), K, N); Beta <- Beta / rowSums(Beta)
  WP <- Th %*% Beta
  W <- t(vapply(seq_len(M), function(i) rmultinom(1, 200, WP[i, ]), numeric(N)))
  fit_no <- sgscatm(W, C, K = K)
  fit_se <- sgscatm(W, C, K = K, se = TRUE)
  expect_null(fit_no$B_z_se)
  expect_equal(dim(fit_se$B_z_se), c(P, K - 1L))
  expect_true(all(is.finite(fit_se$B_z_se)))
})
