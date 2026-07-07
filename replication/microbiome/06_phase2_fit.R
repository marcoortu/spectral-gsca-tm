#!/usr/bin/env Rscript
# ===================================================================
# Phase 2 fit — sgscatm on Zeller CRC genus compositions.
# Gates addressed here: G1 (delocalization), G4 (SE calibration), part of G6.
# ===================================================================
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
source("replication/simulation/sim_utils.R")   # procrustes_align
set.seed(2026)
D <- readRDS("output/phase2_data.rds")
C <- D$C; M <- nrow(C); P <- ncol(C)
cat(sprintf("M=%d samples, P=%d covariates\n", M, P))

# ---- delocalization ratio (two flavours) ----
# (a) fit-based, as pre-registered: r = M*max||z_i||^2/(K-1) on fit scores
r_fit_based <- function(X, K = 5L, lambda = 1) {
  fit <- sgscatm(as.matrix(X), C, K = K, lambda = lambda, rotate = FALSE)
  M * max(rowSums(fit$Z^2)) / (K - 1L)
}
# (b) composition-based natural delocalization on CLR coordinates:
#     r_nat = M*max_i ||clr(x_i)-mean||^2 / (N-1)  (extremity of compositions)
r_clr_based <- function(X) {
  X <- as.matrix(X); Xp <- X + 0.5 / ncol(X)      # pseudocount
  Xp <- Xp / rowSums(Xp)
  L  <- log(Xp); clr <- L - rowMeans(L)
  clr <- sweep(clr, 2, colMeans(clr))
  nrow(X) * max(rowSums(clr^2)) / (ncol(X) - 1L)
}

G1 <- data.frame(
  level = c("genus","genus","genus","species","species"),
  prev  = c(0.10,0.20,0.30,0.10,0.20),
  N     = c(ncol(D$Xg10),ncol(D$Xg20),ncol(D$Xg30),ncol(D$Xsp10),ncol(D$Xsp20)),
  r_fit = c(r_fit_based(D$Xg10), r_fit_based(D$Xg20), r_fit_based(D$Xg30),
            r_fit_based(D$Xsp10), r_fit_based(D$Xsp20)),
  r_clr = c(r_clr_based(D$Xg10), r_clr_based(D$Xg20), r_clr_based(D$Xg30),
            r_clr_based(D$Xsp10), r_clr_based(D$Xsp20))
)
cat("\n=== G1 delocalization ===\n"); print(G1, digits = 4)

# ---- K / lambda calibration (scree gap + covariate R^2) ----
covariate_r2 <- function(Z, C) { H <- C %*% solve(crossprod(C), t(C)); sum((H%*%Z)*Z)/sum(Z*Z) }
Xg <- as.matrix(D$Xg10)
calib <- list(); ii <- 0
for (K in 3:8) for (lam in c(0.5,1,2)) {
  f <- sgscatm(Xg, C, K = K, lambda = lam, rotate = TRUE)
  eg <- f$eigenvalues; gap <- if (length(eg)>=2) (eg[length(eg)-1]-eg[length(eg)])/eg[1] else NA
  ii <- ii+1; calib[[ii]] <- data.frame(K=K, lambda=lam, r2=covariate_r2(f$Z,C), eig_gap=gap)
}
calib <- do.call(rbind, calib)
cat("\n=== K/lambda calibration (lambda=1) ===\n")
print(calib[calib$lambda==1, ], row.names = FALSE, digits = 3)

Kstar <- 5L; lamstar <- 1
cat(sprintf("\nChosen K=%d, lambda=%.1f\n", Kstar, lamstar))

# ---- final fit + timing ----
t0 <- proc.time()[3]
fit <- sgscatm(Xg, C, K = Kstar, lambda = lamstar, rotate = TRUE)
t_fit <- proc.time()[3] - t0
t1 <- proc.time()[3]
vc <- sgscatm_vcov(fit, identified = TRUE)
t_se <- proc.time()[3] - t1
cat(sprintf("sgscatm fit: %.4fs, analytical SE: %.4fs, total: %.4fs\n",
            t_fit, t_se, t_fit + t_se))

Bz <- vc$B; SE <- vc$se                 # standardized, rotation-identified
rownames(Bz) <- rownames(SE) <- colnames(C)
Zscore <- Bz / SE
cat("\n=== Bz (standardized) ===\n"); print(round(Bz,3))
cat("\n=== Bz / SE (z-scores) ===\n"); print(round(Zscore,2))

# per-covariate joint Wald test (rotation-invariant magnitude)
vcm <- vc$vcov            # P(K-1) square, tangent-projected
Km1 <- Kstar - 1L
wald <- data.frame(covariate = colnames(C), stat = NA, df = NA, p = NA)
for (p in seq_len(P)) {
  idx <- (seq_len(Km1)-1L)*P + p          # entries for covariate p across comps
  bp  <- as.vector(Bz[p, ]); Sp <- vcm[idx, idx, drop = FALSE]
  ei  <- eigen(Sp, symmetric = TRUE); pos <- ei$values > 1e-10*max(ei$values)
  Sinv <- ei$vectors[,pos,drop=FALSE] %*% diag(1/ei$values[pos], sum(pos)) %*% t(ei$vectors[,pos,drop=FALSE])
  wald$stat[p] <- as.numeric(t(bp) %*% Sinv %*% bp)
  wald$df[p]   <- sum(pos)
  wald$p[p]    <- pchisq(wald$stat[p], df = sum(pos), lower.tail = FALSE)
}
cat("\n=== per-covariate joint Wald ===\n"); print(wald, row.names=FALSE, digits=3)

# ---- G4: document bootstrap SE, Procrustes+sign aligned to point estimate ----
B_boot <- 300L
Bref <- vc$B
align_to <- function(Bm, Ref){ sv <- svd(crossprod(Bm, Ref)); Bm %*% (sv$u %*% t(sv$v)) }
boot <- matrix(NA_real_, B_boot, P*Km1)
t2 <- proc.time()[3]
for (b in seq_len(B_boot)) {
  idx <- sample.int(M, M, replace = TRUE)
  fb <- tryCatch(sgscatm(Xg[idx,], C[idx,], K=Kstar, lambda=lamstar, rotate=FALSE),
                 error=function(e) NULL)
  if (!is.null(fb)) {
    Bb <- sgscatm_vcov(fb, identified=FALSE)$B
    boot[b,] <- as.vector(align_to(Bb, Bref))
  }
}
t_boot <- proc.time()[3] - t2
boot_sd <- matrix(apply(boot, 2, sd, na.rm=TRUE), P, Km1)
ratio <- SE / boot_sd
cat(sprintf("\n=== G4 SE calibration (bootstrap B=%d, %.1fs) ===\n", B_boot, t_boot))
cat("analytical SE:\n"); print(round(SE,4))
cat("bootstrap SD:\n"); print(round(boot_sd,4))
cat("analytical/bootstrap ratio:\n"); print(round(ratio,3))
cat(sprintf("median ratio = %.3f ; frac within +/-25%% = %.2f\n",
            median(ratio), mean(abs(log(ratio)) < log(1.25))))

# legacy analytical SE for contrast
se_leg <- tryCatch(ilr_se_analytical(fit)$se, error=function(e) NULL)
if (!is.null(se_leg))
  cat(sprintf("legacy ilr_se_analytical median SE = %.3f (vs vcov %.4f) -> ratio %.1f\n",
              median(se_leg), median(SE), median(se_leg)/median(SE)))

saveRDS(list(fit=fit, vc=vc, Bz=Bz, SE=SE, Zscore=Zscore, wald=wald,
             G1=G1, calib=calib, Kstar=Kstar, lamstar=lamstar,
             boot_sd=boot_sd, ratio=ratio, t_fit=t_fit, t_se=t_se, t_boot=t_boot,
             se_leg=se_leg),
        "output/phase2_fit.rds")
cat("\nSaved output/phase2_fit.rds\n")
