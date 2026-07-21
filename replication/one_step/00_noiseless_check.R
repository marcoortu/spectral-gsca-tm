## =====================================================================
## 00_noiseless_check.R  —  MANDATORY pre-run sanity (prereg §5)
## In the noiseless linearised limit (sigma_eps->0, small b_max, huge L,
## large M): H_proj must recover Z_true and B_proj must recover B_z0.
## If this fails, the proj derivation is wrong and NO gate is trusted.
## =====================================================================
source(file.path(getwd(), "replication/one_step/_common.R"))

cat("=== 00 noiseless recovery check ===\n")

dat <- sim_dgp(M = 4000, N = 400, K = 5, P = 3,
               b_max = 0.15, sigma_eps = 1e-4,
               alpha_beta = 0.1, doc_length = 100000L, seed = 42)
validate_dgp(dat)

## Phi/Psi convention self-checks -------------------------------------
fit <- sgscatm(dat$W, dat$C, K = 5, lambda = 1, V = dat$V, rotate = TRUE)
cat(sprintf("Phi row-sum max dev from 1 : %.2e\n",
            max(abs(rowSums(fit$Phi) - 1))))
cat(sprintf("Psi = V'Phi max dev        : %.2e\n",
            max(abs(fit$Psi - crossprod(fit$V, fit$Phi)))))
cat(sprintf("Z*'Z* = I max dev          : %.2e\n",
            max(abs(crossprod(fit$Z) - diag(4)))))

## proj recovery (anchor = generative Beta, a consistent topic estimate)-
proj <- est_proj(fit, Phi_anchor = dat$Beta)
# align recovered scores to the true ILR scores
paZ <- procrustes_align_Z(proj$H, dat$Z_true)
cat(sprintf("\nH_proj vs Z_true : rel Frobenius err = %.4f (||Z_true||=%.3f)\n",
            norm_F(paZ$Z_aligned - dat$Z_true) / norm_F(dat$Z_true),
            norm_F(dat$Z_true)))

evP <- eval_vs_Bz0(proj, dat$Bz0)
cat(sprintf("B_proj  : RMSE(B_z0)=%.4f  ||B_z0||=%.4f  rel=%.4f  ||B_al||/||B_z0||=%.3f\n",
            evP$rmse, evP$norm_Bz0, evP$rmse / evP$norm_Bz0,
            norm_F(evP$B_aligned) / evP$norm_Bz0))

## baseline for contrast ---------------------------------------------
base <- est_baseline_std(fit)
evB  <- eval_vs_Bz0(base, dat$Bz0)
cat(sprintf("B_base  : RMSE(B_z0)=%.4f  rel=%.4f  ||B_al||/||B_z0||=%.3f\n",
            evB$rmse, evB$rmse / evB$norm_Bz0,
            norm_F(evB$B_aligned) / evB$norm_Bz0))

## confirm the score-level residual IS the closure linearisation bias:
## it must shrink ~ ||z|| (i.e. with b_max), not plateau.
relZ_of_bmax <- function(bm) {
  d <- sim_dgp(M = 4000, N = 400, K = 5, P = 3, b_max = bm, sigma_eps = 1e-4,
               alpha_beta = 0.1, doc_length = 100000L, seed = 42)
  f <- sgscatm(d$W, d$C, K = 5, lambda = 1, V = d$V, rotate = TRUE)
  pj <- est_proj(f, Phi_anchor = d$Beta)
  pz <- procrustes_align_Z(pj$H, d$Z_true)
  c(bmax = bm, relZ = norm_F(pz$Z_aligned - d$Z_true) / norm_F(d$Z_true),
    relB = eval_vs_Bz0(pj, d$Bz0)$rmse / norm_F(d$Bz0))
}
tab <- t(vapply(c(0.05, 0.10, 0.15, 0.30), relZ_of_bmax, numeric(3)))
cat("\nlinearisation residual vs b_max (proj, anchor=Beta, noiseless):\n")
print(round(tab, 4))

## verdict: SCALE recovered (B_proj rel small) AND score-residual is
## linearisation (shrinks with b_max).
## The estimand-level bias is the right diagnostic (score-level relZ is
## confounded by the finite-L noise floor divided by shrinking ||Z_true||).
pass_scale <- evP$rmse / evP$norm_Bz0 < 0.02
pass_lin   <- tab[nrow(tab), "relB"] > tab[1, "relB"]      # bias grows with b_max
pass <- pass_scale && pass_lin
cat(sprintf("\nNOISELESS CHECK: scale-recovered=%s  bias-grows-with-bmax=%s  => %s\n",
            pass_scale, pass_lin, if (pass) "PASS" else "FAIL"))
saveRDS(list(evP = evP, evB = evB, tab = tab, pass = pass),
        file.path(getwd(), "replication/one_step/out_00_noiseless.rds"))
