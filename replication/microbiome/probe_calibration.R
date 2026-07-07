#!/usr/bin/env Rscript
# Probe: is sgscatm_vcov analytical SE calibrated against the empirical
# sampling SD across replicates, in the favorable regime?
# Estimand = across-replicate Procrustes-aligned mean (scale-consistent).
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
source("replication/simulation/sim_dgp.R")
source("replication/simulation/sim_utils.R")

R    <- 60L
M    <- 2000L; N <- 200L; K <- 5L; P <- 4L
bmax <- 0.25; lam <- 1
Km1  <- K - 1L

# fixed true coefficient matrix so estimand is stable across replicates
set.seed(100)
Bz0 <- matrix(runif(P * Km1, -bmax, bmax), P, Km1)

fits <- vector("list", R)
for (r in seq_len(R)) {
  dat <- sim_dgp(M = M, N = N, K = K, P = P, Bz0 = Bz0,
                 sigma_eps = 0.3, alpha_beta = 0.05,
                 doc_length = 200L, seed = 5000L + r)
  fit <- sgscatm(dat$W, dat$C, K = K, lambda = lam, rotate = FALSE)
  r_deloc <- M * max(rowSums(fit$Z^2)) / Km1
  fits[[r]] <- list(Bz = fit$Bz, fit = fit, r_deloc = r_deloc,
                    Z = fit$Z)
}

# --- standardized-scale point estimate per replicate (via vcov, no rotation) ---
Bstd <- lapply(fits, function(f) sgscatm_vcov(f$fit, identified = FALSE)$B)

# --- generalized Procrustes to a common reference ---
Ref <- Bstd[[1]]
align_to <- function(B, Ref) { sv <- svd(crossprod(B, Ref)); B %*% (sv$u %*% t(sv$v)) }
for (pass in 1:2) {
  Bal <- lapply(Bstd, align_to, Ref = Ref)
  Ref <- Reduce(`+`, Bal) / length(Bal)
}
Bal <- lapply(Bstd, align_to, Ref = Ref)   # final aligned estimates

# empirical sampling SD (entry-wise) across replicates
arr    <- array(unlist(Bal), dim = c(P, Km1, R))
emp_sd <- apply(arr, c(1, 2), sd)

# analytical SE per replicate, in the aligned basis, then averaged
se_list <- vector("list", R); cov_arr <- array(NA, dim = c(P, Km1, R))
for (r in seq_len(R)) {
  Rr <- { sv <- svd(crossprod(Bstd[[r]], Ref)); sv$u %*% t(sv$v) }
  vc <- sgscatm_vcov(fits[[r]]$fit, rotation = Rr, identified = TRUE)
  se_list[[r]] <- vc$se
  cov_arr[, , r] <- (Ref >= Bal[[r]] - 1.96 * vc$se) &
                    (Ref <= Bal[[r]] + 1.96 * vc$se)
}
se_mean <- Reduce(`+`, se_list) / R

cat(sprintf("cell: M=%d N=%d K=%d b_max=%.2f  reps=%d\n", M, N, K, bmax, R))
cat(sprintf("mean delocalization ratio r = %.3f\n",
            mean(vapply(fits, `[[`, numeric(1), "r_deloc"))))
cat("\nempirical SD (entry-wise):\n");    print(round(emp_sd, 4))
cat("\nmean analytical SE (entry-wise):\n"); print(round(se_mean, 4))
cat("\nSE / SD ratio (entry-wise):\n");    print(round(se_mean / emp_sd, 3))
cat(sprintf("\nmedian SE/SD ratio = %.3f\n", median(se_mean / emp_sd)))
cat(sprintf("mean 95%% coverage of aligned mean = %.3f\n", mean(cov_arr)))
