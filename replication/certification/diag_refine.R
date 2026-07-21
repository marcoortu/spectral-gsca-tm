#!/usr/bin/env Rscript
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
source("replication/simulation/sim_dgp.R")
set.seed(1)
dat <- sim_dgp(M = 1000L, N = 500L, K = 5L, P = 3L, b_max = 0.5,
               sigma_eps = 0.3, alpha_beta = 0.1, doc_length = 400L, seed = 11)
V <- dat$V; C <- scale(dat$C, TRUE, FALSE); Wf <- dat$W / rowSums(dat$W)
rmse <- function(B) sqrt(perm_sign_align(B, dat$Bz0, V)$mse)

fit <- sgscatm(dat$W, dat$C, K = 5L, lambda = 1, rotate = FALSE)
ap  <- sgscatm:::.sg_anchor_pipeline(dat$W, 5L)
ro  <- sgscatm:::.sg_readout_gn(ap$Phi, Wf, V, n_gn = 10L, dz_cap = 1)
cat(sprintf("read-out B RMSE            = %.4f\n", rmse(sgscatm:::.sg_b_step(ro$Z, C))))
gl  <- sgscatm:::.sg_gl_align(fit$Z, ro$Z)
cat(sprintf("oriented-start B(Z0) RMSE  = %.4f\n", rmse(sgscatm:::.sg_b_step(gl$Z, C))))

# frozen-phi: track B RMSE per sweep
Z <- gl$Z; nu <- rep(1e-6, nrow(Z)); Phi <- ap$Phi
cat("\nfrozen-phi refinement, B RMSE by sweep:\n")
for (s in 1:20) {
  zs <- sgscatm:::.sg_z_step(Z, Phi, Wf, V, lambda = 0, CB = NULL, nu = nu,
                             n_gn = 2L, dz_cap = 1)
  Z <- zs$Z; nu <- zs$nu
  B <- sgscatm:::.sg_b_step(Z, C)
  if (s %in% c(1,2,3,4,5,7,10,15,20))
    cat(sprintf("  sweep %2d: RMSE=%.4f  ||B||=%.2f  n_fail=%d\n",
                s, rmse(B), sqrt(sum(B^2)), zs$n_fail))
}

# joint (V4-style)
Zj <- gl$Z; nuj <- rep(1e-6, nrow(Zj)); Phij <- sgscatm:::.sg_phi_step(Zj, Wf, V)
cat("\njoint refinement, B RMSE by sweep:\n")
for (s in 1:20) {
  zs <- sgscatm:::.sg_z_step(Zj, Phij, Wf, V, lambda = 0, CB = NULL, nu = nuj,
                             n_gn = 2L, dz_cap = 1)
  Zj <- zs$Z; nuj <- zs$nu
  Phij <- sgscatm:::.sg_phi_step(Zj, Wf, V)
  B <- sgscatm:::.sg_b_step(Zj, C)
  if (s %in% c(1,2,3,4,5,7,10,15,20))
    cat(sprintf("  sweep %2d: RMSE=%.4f  ||B||=%.2f\n", s, rmse(B), sqrt(sum(B^2))))
}
