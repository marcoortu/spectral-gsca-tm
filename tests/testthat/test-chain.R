test_that("sgscatm_chain runs, collapses the raw pilot, and beats it", {
  skip_on_cran()
  source(testthat::test_path("..", "..", "replication", "simulation", "sim_dgp.R"),
         local = TRUE)
  Bz0 <- matrix(c(0.40,-0.20,0.10,0.30, -0.15,0.35,-0.25,0.05,
                  0.20,0.10,0.40,-0.30), nrow = 3, byrow = TRUE)
  dat <- sim_dgp(M = 1500L, N = 500L, K = 5L, P = 3L, Bz0 = Bz0,
                 sigma_eps = 0.3, alpha_beta = 0.05, doc_length = 5000L,
                 seed = 4242L)
  ch <- sgscatm_chain(dat$W, dat$C, K = 5L, refine = "joint", max_sweeps = 60L)

  expect_s3_class(ch, "sgscatm_chain")
  expect_equal(dim(ch$Bz), c(3L, 4L))
  expect_equal(dim(ch$Theta), c(1500L, 5L))
  expect_true(all(abs(rowSums(ch$Theta) - 1) < 1e-8))

  nB0 <- sqrt(sum(Bz0^2))
  raw_ratio <- sqrt(sum(ch$pilot_Bz^2)) / nB0        # Corollary 11 collapse
  expect_lt(raw_ratio, 0.2)

  chain_rmse <- sqrt(perm_sign_align(ch$Bz, Bz0, dat$V)$mse)
  pilot_rmse <- sqrt(perm_sign_align(ch$pilot_Bz, Bz0, dat$V)$mse)
  expect_lt(chain_rmse, pilot_rmse)                  # orientation/scale restored

  # sandwich SE well-defined and vcov S3 dispatch works
  Sig <- vcov(ch)
  expect_equal(dim(Sig), c(12L, 12L))
  expect_gte(min(eigen(Sig, symmetric = TRUE, only.values = TRUE)$values), -1e-8)
})

test_that("frozen_phi refinement confirms Proposition 21 (best early, then drifts)", {
  skip_on_cran()
  source(testthat::test_path("..", "..", "replication", "simulation", "sim_dgp.R"),
         local = TRUE)
  Bz0 <- matrix(c(0.4,-0.2,0.1,0.3, -0.15,0.35,-0.25,0.05, 0.2,0.1,0.4,-0.3),
                nrow = 3, byrow = TRUE)
  dat <- sim_dgp(M = 1000L, N = 500L, K = 5L, P = 3L, Bz0 = Bz0,
                 sigma_eps = 0.3, alpha_beta = 0.1, doc_length = 400L, seed = 11L)
  ch1 <- sgscatm_chain(dat$W, dat$C, K = 5L, refine = "frozen_phi", max_sweeps = 1L)
  ch9 <- sgscatm_chain(dat$W, dat$C, K = 5L, refine = "frozen_phi", max_sweeps = 15L)
  r1 <- sqrt(perm_sign_align(ch1$Bz, Bz0, dat$V)$mse)
  r9 <- sqrt(perm_sign_align(ch9$Bz, Bz0, dat$V)$mse)
  expect_lt(r1, r9)   # over-refining at frozen Phi degrades B (Prop 21)
})
