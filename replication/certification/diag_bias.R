#!/usr/bin/env Rscript
# Decisive A'/B' diagnostic. Part 1: oracle bias isolation (frozen k-curve vs
# joint@conv, L & M scaling). Part 2/3: feasible frozen@k* + jackknife SE.
# Metric: perm+sign align to Bz0 (never Procrustes). See PREREG_DIAG.md.
suppressPackageStartupMessages(devtools::load_all(".", quiet = TRUE))
source("replication/simulation/sim_dgp.R")
sg <- getNamespace("sgscatm")
K <- 5L; P <- 4L
set.seed(7)
Bz0 <- matrix(c(0.40,-0.20,0.10,0.30, -0.15,0.35,-0.25,0.05,
                0.20,0.10,0.40,-0.30, 0.05,-0.30,0.15,0.25), nrow = P, byrow = TRUE)
V <- ilr_contrast(K)
libfun <- function(mu) function(m) pmax(rnbinom(m, size = 3, mu = mu), 500L)

# split-document jackknife (Phi fixed) -> bias-corrected B and variance add
jackknife <- function(W, C, Vc, Z_full, Phi, seed) {
  set.seed(seed)
  A <- matrix(rbinom(length(W), as.vector(W), 0.5), nrow(W), ncol(W)); B <- W - A
  zfit <- function(Wp) { Wf <- Wp / pmax(rowSums(Wp), 1)
    sg$.sg_z_step(Z_full, Phi, Wf, Vc, 0, NULL, rep(1e-6, nrow(Wp)),
                  n_gn = 8L, dz_cap = 1)$Z }
  BA <- sg$.sg_b_step(zfit(A), C); BB <- sg$.sg_b_step(zfit(B), C)
  Bf <- sg$.sg_b_step(Z_full, C)
  list(B_jk = 2 * Bf - (BA + BB) / 2, var_add = (BA - BB)^2 / 4)
}

# aggregate per-rep aligned estimates + SEs into the metric row
agg <- function(Bal, SE, tag) {
  R <- length(Bal); arr <- array(unlist(Bal), c(P, K-1L, R))
  Bbar <- apply(arr, c(1,2), mean); emp_sd <- apply(arr, c(1,2), sd)
  se_mean <- Reduce(`+`, SE) / R
  true_cov <- mean(mapply(function(b, s) mean(abs(b - Bz0) <= 1.96 * s), Bal, SE))
  br_cov <- mean(mapply(function(b, s) mean(abs(b - Bbar) <= 1.96 * s), Bal, SE))
  data.frame(tag = tag, bias_norm = sqrt(sum((Bbar - Bz0)^2)),
             se_sd = median(se_mean / emp_sd), true_cov = true_cov,
             biasrem_cov = br_cov, rmse = sqrt(mean((arr - as.vector(Bz0))^2)))
}

## ================= Part 1: oracle bias isolation =================
oracle_cell <- function(M, Lmu, R) {
  ks <- c(1L,2L,3L,5L,10L)
  froz <- setNames(lapply(ks, function(k) list(Bal=list(), SE=list())), paste0("f",ks))
  jnt  <- list(Bal = list(), SE = list())
  for (r in seq_len(R)) {
    dat <- sim_dgp(M = M, N = 500L, K = K, P = P, Bz0 = Bz0, sigma_eps = 0.3,
                   alpha_beta = 0.05, doc_length = libfun(Lmu), seed = 70000L + r)
    Wf <- dat$W / rowSums(dat$W); C <- scale(dat$C, TRUE, FALSE)
    fit <- sgscatm(dat$W, dat$C, K = K, lambda = 1, rotate = FALSE)
    Z0  <- sg$.sg_gl_align(fit$Z, scale(dat$Z_true, TRUE, FALSE))$Z
    # frozen k-curve with TRUE Phi
    Z <- Z0; nu <- rep(1e-6, M)
    for (s in seq_len(max(ks))) {
      Z <- sg$.sg_z_step(Z, dat$Beta, Wf, V, 0, NULL, nu, n_gn = 2L, dz_cap = 1)$Z
      if (s %in% ks) {
        B <- sg$.sg_b_step(Z, C); cv <- perm_sign_coverage(B, Bz0, sg$.sg_sandwich(Z,C,B), V)
        key <- paste0("f", s)
        froz[[key]]$Bal <- c(froz[[key]]$Bal, list(perm_sign_align(B,Bz0,V)$B))
        froz[[key]]$SE  <- c(froz[[key]]$SE, list(cv$se))
      }
    }
    # joint@conv with TRUE Phi start
    rj <- sg$.sg_refine(Z0, sg$.sg_phi_step(Z0, Wf, V), Wf, C, V, mode="joint", max_sweeps=40L)
    Bj <- sg$.sg_b_step(rj$Z, C); cvj <- perm_sign_coverage(Bj, Bz0, sg$.sg_sandwich(rj$Z,C,Bj), V)
    jnt$Bal <- c(jnt$Bal, list(perm_sign_align(Bj,Bz0,V)$B)); jnt$SE <- c(jnt$SE, list(cvj$se))
  }
  rows <- lapply(names(froz), function(k) cbind(M=M, L=Lmu, agg(froz[[k]]$Bal, froz[[k]]$SE, k)))
  rbind(do.call(rbind, rows), cbind(M=M, L=Lmu, agg(jnt$Bal, jnt$SE, "joint")))
}

cat("== Part 1: oracle bias isolation ==\n")
p1 <- list()
for (Lmu in c(1e3, 1e4)) for (M in c(1000L, 2000L, 5000L)) {
  R <- if (Lmu == 1e4 && M == 2000L) 80L else 50L
  p1[[length(p1)+1L]] <- oracle_cell(M, Lmu, R)
  cat(sprintf("  done oracle L=%.0e M=%d\n", Lmu, M))
}
P1 <- do.call(rbind, p1); saveRDS(P1, "output/diag_P1_oracle.rds")
print(P1[, c("L","M","tag","bias_norm","se_sd","true_cov")], row.names=FALSE, digits=3)

## ================= Part 2/3: feasible frozen@k* + jackknife =================
KSTAR <- 1L
cat(sprintf("\n== Part 2/3: feasible frozen@k*=%d + jackknife (N=200, L=1e4) ==\n", KSTAR))
feasible_cell <- function(M, R) {
  bp <- bj <- list(); sp <- sj <- list()
  for (r in seq_len(R)) {
    dat <- sim_dgp(M = M, N = 200L, K = K, P = P, Bz0 = Bz0, sigma_eps = 0.3,
                   alpha_beta = 0.05, doc_length = libfun(1e4), seed = 71000L + r)
    ch <- sgscatm_chain(dat$W, dat$C, K = K, refine = "frozen_phi", max_sweeps = KSTAR)
    C <- ch$C_centred
    Sig <- vcov(ch)
    cvp <- perm_sign_coverage(ch$Bz, Bz0, Sig, V)
    jk <- jackknife(dat$W, C, V, ch$Z, ch$anchor_Phi, 71000L + r + 500L)
    Sig_infl <- Sig + diag(as.vector(jk$var_add))
    cvj <- perm_sign_coverage(jk$B_jk, Bz0, Sig_infl, V)
    bp <- c(bp, list(perm_sign_align(ch$Bz,Bz0,V)$B)); sp <- c(sp, list(cvp$se))
    bj <- c(bj, list(perm_sign_align(jk$B_jk,Bz0,V)$B)); sj <- c(sj, list(cvj$se))
  }
  rbind(cbind(M=M, SE="plain",     agg(bp, sp, "feas-frozen")),
        cbind(M=M, SE="jackknife", agg(bj, sj, "feas-frozen-jk")))
}
p2 <- list()
for (M in c(1000L, 2000L, 5000L)) {
  R <- if (M == 2000L) 80L else 50L
  p2[[length(p2)+1L]] <- feasible_cell(M, R); cat(sprintf("  done feasible M=%d\n", M))
}
P2 <- do.call(rbind, p2); saveRDS(P2, "output/diag_P2_feasible.rds")
print(P2[, c("M","SE","bias_norm","se_sd","true_cov","biasrem_cov")], row.names=FALSE, digits=3)
cat("\ndiag_bias DONE\n")
