#!/usr/bin/env Rscript
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
source("replication/simulation/sim_dgp.R")
source("replication/simulation/sim_utils.R")

set.seed(1)
dat <- sim_dgp(M = 500, N = 200, K = 5, P = 4, b_max = 0.5,
               sigma_eps = 0.3, alpha_beta = 0.05, doc_length = 200L, seed = 7)
cat(sprintf("W: %d x %d, C: %d x %d\n", nrow(dat$W), ncol(dat$W), nrow(dat$C), ncol(dat$C)))

t0 <- proc.time()[3]
fit <- sgscatm(dat$W, dat$C, K = 5, lambda = 1, rotate = TRUE)
cat(sprintf("fit time: %.3fs\n", proc.time()[3] - t0))

# Procrustes align estimate to truth
pa <- procrustes_align(fit$Bz, dat$Bz0)
cat(sprintf("RMSE(Bz aligned) = %.4f\n", sqrt(pa$mse)))

# Delocalization ratio r = M * max_i ||z_i||^2 / (K-1)
zn2 <- rowSums(fit$Z^2)
M <- nrow(fit$Z); Km1 <- fit$K - 1L
r_deloc <- M * max(zn2) / Km1
cat(sprintf("delocalization ratio r = %.3f  (max||z||=%.3f)\n", r_deloc, sqrt(max(zn2))))

# Primary analytical SE via sgscatm_vcov, aligned basis
vc <- sgscatm_vcov(fit, rotation = pa$R, identified = TRUE)
cat("sgscatm_vcov SE (aligned):\n"); print(round(vc$se, 4))
z_alpha <- qnorm(0.975)
cover <- (dat$Bz0 >= pa$Bz_aligned - z_alpha*vc$se) &
         (dat$Bz0 <= pa$Bz_aligned + z_alpha*vc$se)
cat(sprintf("vcov coverage (this replicate, 20 entries): %.2f\n", mean(cover)))

# Legacy analytical SE for comparison
se_leg <- tryCatch(ilr_se_analytical(fit), error = function(e) NULL)
if (!is.null(se_leg)) {
  se_leg_al <- .rotate_se(se_leg, pa$R, ncol(dat$C), Km1)
  cat("legacy ilr_se_analytical SE (aligned):\n"); print(round(se_leg_al, 4))
  cat(sprintf("legacy/vcov SE ratio (median): %.2f\n",
              median(se_leg_al / vc$se, na.rm = TRUE)))
}
cat("SMOKE OK\n")
