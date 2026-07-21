## =====================================================================
## 01_estimand_check.R — reconcile B_z0 with Sigma_C^{-1} E(C' U_0)
## (JASA referee point). Prereg §3.
## =====================================================================
source(file.path(getwd(), "replication/one_step/_common.R"))
cat("=== 01 estimand reconciliation ===\n")

# Large single realisation to approximate population moments.
dat <- sim_dgp(M = 200000, N = 50, K = 4, P = 3, b_max = 0.5, sigma_eps = 0.3,
               alpha_beta = 0.1, doc_length = 50L, seed = 99)
C <- dat$C; M <- nrow(C)
SigC <- crossprod(C) / M
Zt   <- dat$Z_true                      # U_0 = true ILR scores

# (1) reviewer form with U_0 = true ILR scores:  Sigma_C^{-1} (1/M) C' Z_true
B_form <- solve(SigC, crossprod(C, Zt) / M)
cat("\n||Sigma_C^{-1} M^{-1} C'Z_true  -  B_z0|| / ||B_z0|| =",
    round(norm_F(B_form - dat$Bz0) / norm_F(dat$Bz0), 4), "\n")
cat("  (coincide exactly in population; residual is O(M^{-1/2}) sampling)\n")

# (2) E[C'U_0] is NOT zero though E[C]=0: show Cov(C,Z) = Sigma_C B_z0
CtZ <- crossprod(C, Zt) / M
cat("\n||M^{-1}C'Z_true - Sigma_C B_z0|| / ||.|| =",
    round(norm_F(CtZ - SigC %*% dat$Bz0) / norm_F(SigC %*% dat$Bz0), 4),
    " (E[C]=0 does NOT imply E[C'U_0]=0)\n")
cat("mean(C) per column:", round(colMeans(C), 5),
    " | ||E[C'U_0]||_F =", round(norm_F(CtZ), 4), "\n")

# (3) standardized eigen-score U~ gives B_z0 R S^{-1}, differs by scale S:
fit <- sgscatm(dat$W, dat$C, K = 4, lambda = 1, V = dat$V, rotate = TRUE)
Ztil <- sqrt(M) * fit$Z                                  # standardized scores
B_std <- solve(SigC, crossprod(C, Ztil) / M)
pa_std <- procrustes_align(B_std, dat$Bz0)
cat("\nstandardized-score form: RMSE(B_z0)/||B_z0|| =",
    round(sqrt(mean((pa_std$Bz_aligned - dat$Bz0)^2)) / norm_F(dat$Bz0), 3),
    " ||B_al||/||B_z0|| =", round(norm_F(pa_std$Bz_aligned)/norm_F(dat$Bz0), 3),
    "\n  (differs from B_z0 by the discarded scale S -- the G1 object)\n")

saveRDS(list(B_form = B_form, Bz0 = dat$Bz0, CtZ = CtZ, SigC = SigC),
        file.path(getwd(), "replication/one_step/out_01_estimand.rds"))
cat("\nCONCLUSION: with U_0 = generative ILR scores, Sigma_C^{-1}E(C'U_0)=B_z0\n",
    "exactly; with U_0 = standardized eigen-scores it returns B_z0 R S^{-1}.\n")
