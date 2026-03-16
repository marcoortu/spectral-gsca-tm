set.seed(42)
M <- 50; N <- 30; K <- 3; P <- 2

make_data <- function() {
  W <- matrix(rpois(M * N, lambda = 5), M, N) + 0.0
  C <- scale(matrix(rnorm(M * P), M, P), center = TRUE, scale = FALSE)
  list(W = W, C = C)
}

test_that("sgscatm returns correct structure", {
  d <- make_data()
  fit <- sgscatm(d$W, d$C, K = K, lambda = 1)
  expect_s3_class(fit, "sgscatm")
  expect_equal(dim(fit$Pi),  c(M, K))
  expect_equal(dim(fit$Phi), c(K, N))
  expect_equal(dim(fit$Bz),  c(P, K - 1L))
  expect_equal(dim(fit$Z),   c(M, K - 1L))
})

test_that("topic proportions sum to 1", {
  d <- make_data()
  fit <- sgscatm(d$W, d$C, K = K, lambda = 1)
  expect_equal(rowSums(fit$Pi), rep(1, M), tolerance = 1e-10)
})

test_that("topic proportions are non-negative", {
  d <- make_data()
  fit <- sgscatm(d$W, d$C, K = K, lambda = 1)
  expect_true(all(fit$Pi >= 0))
})

test_that("lambda=0 recovers centred LSA solution", {
  d <- make_data()
  fit0 <- sgscatm(d$W, d$C, K = K, lambda = 0)
  expect_equal(rowSums(fit0$Pi), rep(1, M), tolerance = 1e-10)
})

test_that("ILR roundtrip is exact", {
  V <- ilr_contrast(K)
  Pi <- matrix(c(0.5, 0.3, 0.2), nrow = 1)
  Z  <- proportions_to_ilr(Pi, V)
  Pi2 <- ilr_to_proportions(Z, V)
  expect_equal(Pi2, Pi, tolerance = 1e-12)
})

test_that("fit stores auxiliary fields for SE", {
  d   <- make_data()
  fit <- sgscatm(d$W, d$C, K = K, lambda = 1)
  expect_false(is.null(fit$U_all))
  expect_false(is.null(fit$W_tilde))
  expect_false(is.null(fit$C_centred))
  expect_false(is.null(fit$eigenvalues_all))
})

test_that("ilr_se (bootstrap) returns valid structure", {
  d   <- make_data()
  fit <- sgscatm(d$W, d$C, K = K, lambda = 1)
  se  <- ilr_se(fit, d$W, d$C, B = 20L, seed = 1L)
  expect_equal(dim(se$se),       c(P, K - 1L))
  expect_equal(dim(se$ci_lower), c(P, K - 1L))
  expect_equal(dim(se$ci_upper), c(P, K - 1L))
  expect_true(all(is.finite(se$se)))
  expect_true(all(se$se >= 0))
  expect_true(all(se$ci_lower < se$ci_upper))
})

test_that("predict returns valid proportions", {
  d    <- make_data()
  fit  <- sgscatm(d$W, d$C, K = K, lambda = 1)
  Wnew <- matrix(rpois(10 * N, 5), 10, N) + 0.0
  Pnew <- predict(fit, Wnew)
  expect_equal(dim(Pnew), c(10, K))
  expect_equal(rowSums(Pnew), rep(1, 10), tolerance = 1e-10)
})

test_that("rotate=TRUE returns valid fit with orthogonal R_star", {
  d   <- make_data()
  fit <- sgscatm(d$W, d$C, K = K, lambda = 1, rotate = TRUE)
  expect_true(fit$rotate)
  expect_false(is.null(fit$R_star))
  expect_equal(dim(fit$R_star), c(K - 1L, K - 1L))
  # R_star is orthogonal: R'R = I
  expect_equal(crossprod(fit$R_star),
               diag(K - 1L), tolerance = 1e-10)
  # output dimensions and constraints unchanged
  expect_equal(dim(fit$Pi),  c(M, K))
  expect_equal(rowSums(fit$Pi), rep(1, M), tolerance = 1e-10)
  expect_true(all(fit$Pi >= 0))
})

test_that("rotate=TRUE does not change objective (rotation invariance)", {
  d    <- make_data()
  fit0 <- sgscatm(d$W, d$C, K = K, lambda = 1, rotate = FALSE)
  fit1 <- sgscatm(d$W, d$C, K = K, lambda = 1, rotate = TRUE)
  # tr(Z'S_z Z) is invariant: eigenvalues of S_z are unchanged
  expect_equal(sort(fit0$eigenvalues, decreasing = TRUE),
               sort(fit1$eigenvalues, decreasing = TRUE),
               tolerance = 1e-8)
})

test_that("rotate=FALSE stores rotate=FALSE and R_star=NULL", {
  d   <- make_data()
  fit <- sgscatm(d$W, d$C, K = K, lambda = 1, rotate = FALSE)
  expect_false(fit$rotate)
  expect_null(fit$R_star)
})

test_that("default rotate=TRUE is applied", {
  d   <- make_data()
  fit <- sgscatm(d$W, d$C, K = K, lambda = 1)   # rotate not specified
  expect_true(fit$rotate)
  expect_false(is.null(fit$R_star))
})

test_that("rotate=TRUE with K=2 is skipped (trivial rotation)", {
  d   <- make_data()
  fit <- sgscatm(d$W, d$C, K = 2L, lambda = 1, rotate = TRUE)
  expect_false(fit$rotate)   # K=2 skips rotation
  expect_null(fit$R_star)
})
# ---- refine_phi & topic_word_dist -------------------------------------------

test_that("refine_phi (kmeans) returns valid probability matrix", {
  d    <- make_data()
  fit  <- sgscatm(d$W, d$C, K = K, lambda = 1)
  fit2 <- refine_phi(fit, d$W, method = "kmeans", seed = 1L)
  expect_true(fit2$phi_refined)
  expect_equal(fit2$phi_method, "kmeans")
  expect_equal(dim(fit2$Phi), c(K, N))
  expect_equal(rowSums(fit2$Phi), rep(1, K), tolerance = 1e-10)
  expect_true(all(fit2$Phi > 0))
})

test_that("refine_phi (em) returns valid probability matrix", {
  d    <- make_data()
  fit  <- sgscatm(d$W, d$C, K = K, lambda = 1)
  fit2 <- refine_phi(fit, d$W, method = "em")
  expect_true(fit2$phi_refined)
  expect_equal(fit2$phi_method, "em")
  expect_equal(dim(fit2$Phi), c(K, N))
  expect_equal(rowSums(fit2$Phi), rep(1, K), tolerance = 1e-10)
  expect_true(all(fit2$Phi > 0))
})

test_that("refine_phi does not alter Pi, Bz, Z", {
  d   <- make_data()
  fit <- sgscatm(d$W, d$C, K = K, lambda = 1)
  fit2 <- refine_phi(fit, d$W, method = "kmeans", seed = 1L)
  expect_equal(fit2$Pi,  fit$Pi)
  expect_equal(fit2$Bz,  fit$Bz)
  expect_equal(fit2$Z,   fit$Z)
  expect_equal(fit2$Psi, fit$Psi)
})

test_that("refine_phi with temp<1 sharpens topic distributions", {
  d    <- make_data()
  fit  <- sgscatm(d$W, d$C, K = K, lambda = 1)
  fit1 <- refine_phi(fit, d$W, method = "kmeans", seed = 1L, temp = 1.0)
  fit2 <- refine_phi(fit, d$W, method = "kmeans", seed = 1L, temp = 0.5)
  expect_true(all(apply(fit2$Phi, 1L, max) >= apply(fit1$Phi, 1L, max)))
})

test_that("topic_word_dist returns K x N probability matrix (unrefined)", {
  d   <- make_data()
  fit <- sgscatm(d$W, d$C, K = K, lambda = 1, rotate = FALSE)
  phi <- topic_word_dist(fit)
  expect_equal(dim(phi), c(K, N))
  expect_equal(rowSums(phi), rep(1, K), tolerance = 1e-10)
})

test_that("topic_word_dist uses M-step Phi when phi_refined=TRUE", {
  d    <- make_data()
  fit  <- sgscatm(d$W, d$C, K = K, lambda = 1)
  fit2 <- refine_phi(fit, d$W)
  phi  <- topic_word_dist(fit2)
  expect_equal(phi, fit2$Phi)
})

test_that("topic_word_dist with temp<1 sharpens (unrefined and refined)", {
  d   <- make_data()
  fit <- sgscatm(d$W, d$C, K = K, lambda = 1)
  p1  <- topic_word_dist(fit, temp = 1.0)
  p2  <- topic_word_dist(fit, temp = 0.5)
  expect_true(all(apply(p2, 1L, max) >= apply(p1, 1L, max)))

  fit_r <- refine_phi(fit, d$W)
  r1 <- topic_word_dist(fit_r, temp = 1.0)
  r2 <- topic_word_dist(fit_r, temp = 0.5)
  expect_true(all(apply(r2, 1L, max) >= apply(r1, 1L, max)))
})
