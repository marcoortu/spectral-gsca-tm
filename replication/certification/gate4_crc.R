#!/usr/bin/env Rscript
# Phase A — G4 CRC fork: on real CRC (Zeller 2014), does the plain Lemma-17
# sandwich agree with the full-chain bootstrap? Decides the paper's primary SE.
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
source("replication/simulation/sim_utils.R")   # (perm helpers live in package now)
set.seed(2026)
D <- readRDS("output/phase2_data.rds")
Xg <- as.matrix(D$Xg10); C <- D$C
Xg <- Xg / rowSums(Xg)
counts <- round(Xg * 1e6)                        # chain anchors need counts; large depth = in-regime
M <- nrow(counts); N <- ncol(counts); P <- ncol(C); K <- 5L
cat(sprintf("CRC genus: M=%d samples, N=%d genera, P=%d covariates, K=%d\n", M, N, P, K))

t0 <- proc.time()[3]
ch <- sgscatm_chain(counts, C, K = K, refine = "frozen_phi")
t_fit <- proc.time()[3] - t0
V <- ch$V; Km1 <- K - 1L
cat(sprintf("chain fit: %.2fs (sweeps=%d, rule_stop=%s, monotone=%s)\n",
            t_fit, ch$sweeps, ch$rule_stop, ch$monotone_ok))

# --- SE_sandwich (plain Lemma-17) ---
Sig_s <- vcov(ch)
SE_s  <- matrix(sqrt(pmax(diag(Sig_s), 0)), P, Km1)

# --- SE_boot (full-chain bootstrap) ---
t0 <- proc.time()[3]
bs <- chain_boot_se(ch, counts, C, B = 200L, seed = 99L)
t_boot <- proc.time()[3] - t0
SE_b <- bs$se
# empirical bootstrap covariance of vec(Bz) for the joint Wald
bootmat <- t(apply(bs$boot, 3L, as.vector))     # B x P(K-1)
Sig_b <- cov(bootmat)
cat(sprintf("bootstrap: %.1fs (B=%d ok)\n", t_boot, bs$B))

rownames(SE_s) <- rownames(SE_b) <- colnames(C)
ratio <- SE_s / SE_b
cat("\n=== B_hat (chain, generative scale) ===\n"); print(round(ch$Bz, 3))
cat("\n=== SE_sandwich ===\n"); print(round(SE_s, 4))
cat("\n=== SE_boot ===\n"); print(round(SE_b, 4))
cat("\n=== ratio SE_sandwich/SE_boot ===\n"); print(round(ratio, 3))
med <- median(ratio); pct_within <- mean(abs(log(ratio)) < log(1.25))
percov <- apply(ratio, 1L, median)
cat(sprintf("\nmedian ratio = %.3f ; %% entries within +/-25%% of parity = %.0f%%\n",
            med, 100*pct_within))
cat("per-covariate median ratio:\n"); print(round(percov, 3))

# --- delta_Phi proxy: anchor exclusivity (min_k max_j Phi_kj / colSums) ---
Phi <- ch$anchor_Phi
excl <- min(apply(sweep(Phi, 2L, pmax(colSums(Phi), 1e-300), "/"), 1L, max))
cat(sprintf("\ndelta_Phi proxy: anchor exclusivity = %.3f (near 1 = clean/exclusive anchors)\n", excl))

# --- joint Wald per covariate under both SEs ---
wald_p <- function(Sig) {
  vapply(seq_len(P), function(p) {
    idx <- (seq_len(Km1)-1L)*P + p; bp <- as.vector(ch$Bz[p, ]); Sp <- Sig[idx, idx, drop=FALSE]
    ei <- eigen(Sp, symmetric=TRUE); pos <- ei$values > 1e-10*max(ei$values)
    Si <- ei$vectors[,pos,drop=FALSE] %*% diag(1/ei$values[pos], sum(pos)) %*% t(ei$vectors[,pos,drop=FALSE])
    pchisq(as.numeric(t(bp)%*%Si%*%bp), df=sum(pos), lower.tail=FALSE)
  }, numeric(1))
}
comp <- data.frame(covariate = colnames(C),
                   wald_p_sandwich = wald_p(Sig_s),
                   wald_p_bootstrap = wald_p(Sig_b))
cat("\n=== joint Wald p per covariate, both SEs ===\n"); print(comp, row.names=FALSE, digits=3)

# --- fork decision ---
fork <- if (med >= 0.8 && med <= 1.25 && pct_within >= 0.60) "SANDWICH" else "BOOTSTRAP"
cat(sprintf("\n=== FORK: primary SE = %s (median %.2f, within25 %.0f%%) ===\n",
            fork, med, 100*pct_within))

saveRDS(list(Bz=ch$Bz, SE_s=SE_s, SE_b=SE_b, ratio=ratio, med=med,
             pct_within=pct_within, percov=percov, excl=excl, comp=comp,
             fork=fork, t_fit=t_fit, t_boot=t_boot, M=M, N=N),
        "output/gate4_crc.rds")
cat("\ngate4 DONE\n")
