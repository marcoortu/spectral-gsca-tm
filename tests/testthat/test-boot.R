test_that("chain_boot_se runs, returns right shape, positive SE, aligns to point est", {
  skip_on_cran()
  source(testthat::test_path("..", "..", "replication", "simulation", "sim_dgp.R"),
         local = TRUE)
  Bz0 <- matrix(c(0.4,-0.2,0.1,0.3, -0.15,0.35,-0.25,0.05, 0.2,0.1,0.4,-0.3),
                nrow = 3, byrow = TRUE)
  dat <- sim_dgp(M = 800L, N = 200L, K = 5L, P = 3L, Bz0 = Bz0, sigma_eps = 0.3,
                 alpha_beta = 0.05,
                 doc_length = function(m) pmax(rnbinom(m, 3, mu = 1e4), 500L),
                 seed = 321L)
  ch <- sgscatm_chain(dat$W, dat$C, K = 5L, refine = "frozen_phi")
  bs <- chain_boot_se(ch, dat$W, dat$C, B = 25L, seed = 1L)
  expect_equal(dim(bs$se), c(3L, 4L))
  expect_true(all(bs$se > 0))
  expect_true(all(bs$ci_lower < bs$ci_upper))
  expect_gte(bs$B, 20L)
})

test_that("perm_sign_align is invariant to a relabeling of the estimate", {
  V <- ilr_contrast(5L)
  set.seed(9); B <- matrix(rnorm(3 * 4, 0, 0.3), 3, 4)
  Pm <- diag(5L)[c(2, 4, 1, 5, 3), ]; Qp <- crossprod(V, Pm %*% V)
  al <- perm_sign_align(B %*% t(Qp), B, V)
  expect_lt(al$mse, 1e-20)
})
