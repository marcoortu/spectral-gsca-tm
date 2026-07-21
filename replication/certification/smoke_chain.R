#!/usr/bin/env Rscript
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
source("replication/simulation/sim_dgp.R")
set.seed(1)
# in-regime-ish: moderate L, K=5, P=3
dat <- sim_dgp(M = 1000L, N = 500L, K = 5L, P = 3L, b_max = 0.5,
               sigma_eps = 0.3, alpha_beta = 0.1, doc_length = 400L, seed = 11)
cat("truth ||Bz0||_F =", round(sqrt(sum(dat$Bz0^2)), 3), "\n")

t0 <- proc.time()[3]
ch <- sgscatm_chain(dat$W, dat$C, K = 5L, lambda = 1, refine = "frozen_phi",
                    verbose = TRUE)
cat(sprintf("chain fit time: %.2fs\n", proc.time()[3] - t0))
print(ch)

# perm+sign aligned RMSE vs truth (honest metric)
al <- perm_sign_align(ch$Bz, dat$Bz0, dat$V)
cat(sprintf("\nchain perm+sign RMSE(Bz) = %.4f   (||Bz0||=%.3f)\n",
            sqrt(al$mse), sqrt(sum(dat$Bz0^2))))
# raw pilot collapse
cat(sprintf("raw-pilot norm ratio ||Bpil||/||Bz0|| = %.3f (should be small; collapse)\n",
            sqrt(sum(ch$pilot_Bz^2)) / sqrt(sum(dat$Bz0^2))))

# sandwich SE + coverage of truth (perm+sign)
Sig <- vcov(ch)
cat(sprintf("sandwich PSD (min eig >= 0): %s ; symmetric: %s\n",
            min(eigen(Sig, only.values = TRUE)$values) > -1e-8,
            max(abs(Sig - t(Sig))) < 1e-10))
cv <- perm_sign_coverage(ch$Bz, dat$Bz0, Sig, dat$V)
cat(sprintf("entrywise 95%% coverage of Bz0 (this rep): %.2f\n", mean(cv$covers)))
cat("SE (mean) =", round(mean(cv$se), 4), "\n")

# joint (V4) variant for comparison
ch2 <- sgscatm_chain(dat$W, dat$C, K = 5L, lambda = 1, refine = "joint")
al2 <- perm_sign_align(ch2$Bz, dat$Bz0, dat$V)
cat(sprintf("\njoint(V4) perm+sign RMSE(Bz) = %.4f\n", sqrt(al2$mse)))
cat("SMOKE OK\n")
