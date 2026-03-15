test_that("ilr_contrast produces valid V matrix", {
  for (K in 2:6) {
    V <- ilr_contrast(K)
    expect_equal(dim(V), c(K, K - 1L))
    expect_equal(crossprod(V), diag(K - 1L), tolerance = 1e-12)
    expect_equal(as.vector(colSums(V)), rep(0, K - 1L), tolerance = 1e-12)
  }
})
