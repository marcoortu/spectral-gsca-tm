test_that("anchor pipeline recovers Phi up to topic permutation", {
  skip_on_cran()
  source(testthat::test_path("..", "..", "replication", "simulation", "sim_dgp.R"),
         local = TRUE)
  dat <- sim_dgp(M = 2000L, N = 500L, K = 5L, P = 3L, b_max = 0.5,
                 sigma_eps = 0.3, alpha_beta = 0.05, doc_length = 200L,
                 seed = 55001L)
  ap <- sgscatm:::.sg_anchor_pipeline(dat$W, 5L)
  # row-wise total-variation error minimised over topic permutations
  perms <- sgscatm:::.sg_all_perms(5L)
  tv <- min(apply(perms, 1L, function(p)
    mean(0.5 * rowSums(abs(ap$Phi[p, ] - dat$Beta)))))
  expect_lt(tv, 0.20)
  expect_equal(dim(ap$Phi), c(5L, 500L))
  expect_true(all(abs(rowSums(ap$Phi) - 1) < 1e-8))
})
