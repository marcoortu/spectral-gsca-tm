test_that("Lemma-17 sandwich is symmetric, PSD, and matches an explicit loop", {
  set.seed(1)
  M <- 300L; P <- 3L; Kp <- 4L
  C <- scale(matrix(rnorm(M * P), M, P), TRUE, FALSE)
  B <- matrix(rnorm(P * Kp, 0, 0.3), P, Kp)
  Z <- C %*% B + matrix(rnorm(M * Kp, 0, 0.5), M, Kp)
  Sig <- sgscatm:::.sg_sandwich(Z, C, sgscatm:::.sg_b_step(Z, C))

  expect_equal(Sig, t(Sig), tolerance = 1e-10)
  expect_gte(min(eigen(Sig, symmetric = TRUE, only.values = TRUE)$values), -1e-8)

  # explicit reference: (M/(M-P)) (I x XtXi) [sum r_i r_i' x c_i c_i'] (I x XtXi)
  Bhat <- sgscatm:::.sg_b_step(Z, C)
  R <- Z - C %*% Bhat
  XtXi <- solve(crossprod(C))
  meat <- matrix(0, P * Kp, P * Kp)
  for (i in seq_len(M))
    meat <- meat + kronecker(tcrossprod(R[i, ]), tcrossprod(C[i, ]))
  bread <- kronecker(diag(Kp), XtXi)
  ref <- (M / (M - P)) * (bread %*% meat %*% bread)
  expect_equal(Sig, ref, tolerance = 1e-8)
})

test_that("sandwich SE scales as O(1/sqrt(M))", {
  set.seed(2)
  P <- 2L; Kp <- 2L; B <- matrix(c(0.4, -0.2, 0.1, 0.3), P, Kp)
  se_of <- function(M) {
    C <- scale(matrix(rnorm(M * P), M, P), TRUE, FALSE)
    Z <- C %*% B + matrix(rnorm(M * Kp, 0, 0.5), M, Kp)
    mean(sqrt(diag(sgscatm:::.sg_sandwich(Z, C, sgscatm:::.sg_b_step(Z, C)))))
  }
  r <- se_of(500L) / se_of(2000L)
  expect_gt(r, 1.6)   # ~ sqrt(4) = 2
  expect_lt(r, 2.5)
})
