#!/usr/bin/env Rscript
# Simulation certification (reduced reps): G0 (subspace recovery + raw collapse
# + chain vs oracle), G7 (B-functional start-independence), G2c (in-regime
# large-L coverage of the chain + Lemma-17 sandwich). Metric alignment =
# permutation+sign only (never Procrustes).
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
source("replication/simulation/sim_dgp.R")
sg <- getNamespace("sgscatm")
K <- 5L; P <- 3L; N <- 500L
Bz0 <- matrix(c(0.40,-0.20,0.10,0.30, -0.15,0.35,-0.25,0.05,
                0.20,0.10,0.40,-0.30), nrow = P, byrow = TRUE)
rmse <- function(B, V) sqrt(perm_sign_align(B, Bz0, V)$mse)

principal_angle <- function(Z1, Z2) {
  q1 <- qr.Q(qr(scale(Z1, TRUE, FALSE))); q2 <- qr.Q(qr(scale(Z2, TRUE, FALSE)))
  acos(min(svd(crossprod(q1, q2))$d))
}
# oracle chain: GL-align pilot to TRUE scores, then joint refine (isolates
# refinement from anchor error)
oracle_B <- function(dat) {
  Wf <- dat$W / rowSums(dat$W); V <- dat$V; C <- scale(dat$C, TRUE, FALSE)
  fit <- sgscatm(dat$W, dat$C, K = K, lambda = 1, rotate = FALSE)
  gl <- sg$.sg_gl_align(fit$Z, scale(dat$Z_true, TRUE, FALSE))
  rf <- sg$.sg_refine(gl$Z, sg$.sg_phi_step(gl$Z, Wf, V), Wf, C, V,
                      mode = "joint", max_sweeps = 60L)
  sg$.sg_b_step(rf$Z, C)
}

## ---------- G0 : subspace recovery + collapse + chain vs oracle ----------
cat("== G0 : subspace recovery + raw collapse ==\n")
g0 <- list()
for (M in c(1000L, 2000L, 5000L)) {
  ang <- rr <- ch <- orc <- numeric(0)
  for (rep in 1:4) {
    dat <- sim_dgp(M = M, N = N, K = K, P = P, Bz0 = Bz0, sigma_eps = 0.3,
                   alpha_beta = 0.1, doc_length = 200L, seed = 60000L + M + rep)
    fit <- sgscatm(dat$W, dat$C, K = K, lambda = 1, rotate = FALSE)
    ang <- c(ang, principal_angle(fit$Z, dat$Z_true))
    rr  <- c(rr, sqrt(sum(fit$Bz^2)) / sqrt(sum(Bz0^2)))
    cf  <- sgscatm_chain(dat$W, dat$C, K = K, refine = "joint", max_sweeps = 60L)
    ch  <- c(ch, rmse(cf$Bz, dat$V))
    orc <- c(orc, rmse(oracle_B(dat), dat$V))
  }
  g0[[as.character(M)]] <- c(M = M, angle = mean(ang), raw_ratio = mean(rr),
                             chain_rmse = mean(ch), oracle_rmse = mean(orc))
  cat(sprintf("  M=%5d angle=%.2e raw_ratio=%.3f (M^-1/2=%.3f) chain=%.4f oracle=%.4f\n",
              M, mean(ang), mean(rr), 1/sqrt(M)*sqrt(1000), mean(ch), mean(orc)))
}
saveRDS(g0, "output/cert_G0.rds")

## ---------- G7 : B-functional start-independence ----------
cat("\n== G7 : start-independence (pilot-start vs truth-start joint refine) ==\n")
g7 <- numeric(0)
for (rep in 1:4) {
  dat <- sim_dgp(M = 2000L, N = N, K = K, P = P, Bz0 = Bz0, sigma_eps = 0.3,
                 alpha_beta = 0.1, doc_length = 1e4, seed = 61000L + rep)
  Wf <- dat$W / rowSums(dat$W); V <- dat$V; C <- scale(dat$C, TRUE, FALSE)
  cf <- sgscatm_chain(dat$W, dat$C, K = K, refine = "joint", max_sweeps = 60L)
  rt <- sg$.sg_refine(scale(dat$Z_true, TRUE, FALSE),
                      sg$.sg_phi_step(scale(dat$Z_true,TRUE,FALSE), Wf, V),
                      Wf, C, V, mode = "joint", max_sweeps = 60L)
  # align both to Bz0 frame, compare
  a1 <- perm_sign_align(cf$Bz, Bz0, V)$B
  a2 <- perm_sign_align(sg$.sg_b_step(rt$Z, C), Bz0, V)$B
  g7 <- c(g7, mean((a1 - a2)^2))
}
cat(sprintf("  mean B mse(pilot-start vs truth-start) = %.2e (gate <= 5e-4)\n", mean(g7)))
saveRDS(g7, "output/cert_G7.rds")

## ---------- G2c : in-regime large-L coverage ----------
cat("\n== G2c : in-regime coverage (large L), chain + Lemma-17 sandwich ==\n")
for (Lval in c(1e4)) {
  for (M in c(2000L)) {
    R <- 40L
    Bs <- vector("list", R); SEs <- vector("list", R); covs <- numeric(R)
    ratios <- numeric(R)
    for (rep in 1:R) {
      dat <- sim_dgp(M = M, N = N, K = K, P = P, Bz0 = Bz0, sigma_eps = 0.3,
                     alpha_beta = 0.05, doc_length = as.integer(Lval),
                     seed = 62000L + rep)
      cf <- sgscatm_chain(dat$W, dat$C, K = K, refine = "joint", max_sweeps = 60L)
      al <- perm_sign_align(cf$Bz, Bz0, dat$V)
      Bs[[rep]] <- al$B
      Sig <- vcov(cf)
      cv  <- perm_sign_coverage(cf$Bz, Bz0, Sig, dat$V)
      covs[rep] <- mean(cv$covers); SEs[[rep]] <- cv$se
    }
    arr <- array(unlist(Bs), c(P, K-1L, R))
    emp_sd <- apply(arr, c(1,2), sd)
    se_mean <- Reduce(`+`, SEs) / R
    cat(sprintf("  L=%.0e M=%d : coverage=%.3f  SE/SD median=%.3f (gate cov[.90,.97], SE/SD[.9,1.15])\n",
                Lval, M, mean(covs), median(se_mean / emp_sd)))
    saveRDS(list(cov = mean(covs), se_sd = median(se_mean/emp_sd),
                 emp_sd = emp_sd, se_mean = se_mean, covs = covs),
            sprintf("output/cert_G2c_L%.0e_M%d.rds", Lval, M))
  }
}
cat("\ncert_sim DONE\n")
